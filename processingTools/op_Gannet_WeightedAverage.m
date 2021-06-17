function[MRS_struct] = op_Gannet_WeightedAverage(MRS_struct)
%op_Gannet_SpectralAlign(MRS_struct.out.SpecReg.data_align, water_flag, MRS_struct);
%adapted from Gannet 3.1

if MRS_struct.flags.averaged || MRS_struct.averages<2
    %DO NOTHING
    disp('WARNING: No averages found. Returning input without modification!');
    return;
end

if MRS_struct.dims.averages==0
    %DO NOTHING
    disp('WARNING: No averages found. Returning input without modification!');
    return;
end

fids = MRS_struct.fids;
water_flag = MRS_struct.out.RSpecReg.water_flag;

freqRange = MRS_struct.spectralwidth/(MRS_struct.txfrq*1e-6);
freq = (size(fids,1) + 1 - (1:size(fids,1))) / size(fids,1) * freqRange + 4.68 - freqRange/2;

D = zeros(MRS_struct.averages); %NtotalSpec/N = Nav

time = (0:(MRS_struct.sz(1)-1))'/MRS_struct.spectralwidth;
tMax = find(time <= 0.1,1,'last');

% Run similarity metric
for kk = 1:MRS_struct.averages
    for ll = 1:MRS_struct.averages
        tmp = sum((real(fids(1:tMax,kk)) - real(fids(1:tMax,ll))).^2) / 200;
        if tmp == 0
            D(kk,ll) = NaN;
        else
            D(kk,ll) = tmp;
        end
    end
end
d = nanmean(D);
w = 1./d.^2;
w = w/sum(w);
w = repmat(w, [size(fids,1) 1]);

if water_flag
    dataLim = ceil(size(fids,2)/3);
    MRS_struct.fids = sum(w(:,1:dataLim) .* fids(:,1:dataLim),2);
else
    MRS_struct.fids = sum(w .* fids,2);
end


MRS_struct.specs = fftshift(ifft(MRS_struct.fids,[],MRS_struct.dims.t),MRS_struct.dims.t);

% update dims
if MRS_struct.dims.t>MRS_struct.dims.averages
    dims.t=MRS_struct.dims.t-1;
else
    dims.t=MRS_struct.dims.t;
end
if MRS_struct.dims.coils>MRS_struct.dims.averages
    dims.coils=MRS_struct.dims.coils-1;
else
    dims.coils=MRS_struct.dims.coils;
end
dims.averages=0;
if MRS_struct.dims.subSpecs>MRS_struct.dims.averages
    dims.subSpecs=MRS_struct.dims.subSpecs-1;
else
    dims.subSpecs=MRS_struct.dims.subSpecs;
end
if MRS_struct.dims.extras>MRS_struct.dims.averages
    dims.extras=MRS_struct.dims.extras-1;
else
    dims.extras=MRS_struct.dims.extras;
end
MRS_struct.dims=dims;

% correct averages and sz
MRS_struct.sz = size(MRS_struct.   fids);
MRS_struct.averages=1;

% update flags
MRS_struct.flags.writtentostruct=1;
MRS_struct.flags.averaged=1;

end%func