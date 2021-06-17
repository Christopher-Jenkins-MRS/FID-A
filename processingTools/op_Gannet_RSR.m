function [MRS_struct] = op_Gannet_RSR(MRS_struct)
% function [MRS_struct] = op_Gannet_RSR(MRS_struct)
% Robust Spectral registration.Adapted from Gannet 3.1 (MM: August 2019)


%% Initialise and Pre-allocate memory
params = zeros(MRS_struct.averages,2); % size = NtotalSpec/NSubSpec=Nav, 2 (2 for f & phi)
MSE    = zeros(MRS_struct.averages,1); % size = NtotalSpec/NSubSpec=Nav, 1
w      = zeros(1,MRS_struct.averages); % size = 1 ,NtotalSpec/NSubSpec=Nav

MRS_struct.out.RSpecReg.freq(1,:)  = zeros(1,MRS_struct.sz(2)); % all zeros(1,NtotalSpec)
MRS_struct.out.RSpecReg.phase(1,:) = zeros(1,MRS_struct.sz(2));
MRS_struct.out.RSpecReg.MSE(1,:)   = zeros(1,MRS_struct.sz(2));

DataToAlign = complex(zeros(size(MRS_struct.fids))); % empty complex matrix size same size as fids

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Algorithm','levenberg-marquardt','Display','off');

% Automatic lipid/unstable residual water removal
freqRange = MRS_struct.spectralwidth/(MRS_struct.txfrq*1e-6); % double: size of window
freq = (MRS_struct.sz(1) + 1 - (1:MRS_struct.sz(1))) / MRS_struct.sz(1) * freqRange + 4.68 - freqRange/2; % Calculate frequency ppm


waterLim = freq <= 4.68 + 0.25 & freq >= 4.68 - 0.25; %Hard code frequency limits
lipidLim = freq <= 1.85 & freq >= 0;
noiseLim = freq <= 11 & freq >= 10;

%% Determine lipid and residual water contamination & filter if found

S = mean(real(fftshift(ifft(double(MRS_struct.fids),[],1),1)),2); %mean of all ffts
r = std(S(lipidLim)) / std(S(noiseLim)); % ~ SNR of lipid region
r_threshold = 40; %Threshold SNR hard coded for this region

spec = real(fftshift(ifft(double(MRS_struct.fids),[],1),1));

q = sum(abs(spec(waterLim,:))) * abs(freq(1) - freq(2)); %sum of magnitude water signal
q = q / max(q); %normalised to max
q = sum(q < 0.5) / length(q); %sum all q below 50% max / NtotalSpec
q_threshold = 0.1;

lipid_flag = 0;
water_flag = 0;

% If lipid and residual water are above threshold, then filter data.
if r > r_threshold || q > q_threshold
    if r > r_threshold
        lipid_flag = 1;
    end
    if q > q_threshold
        water_flag = 1;
    end
    
    spec = fftshift(ifft(double(MRS_struct.fids),[],1),1); % overwrite real spec with complex
    
    reverseStr = '';
    for jj = 1:size(MRS_struct.fids,2) % for each FId,
        if lipid_flag && ~water_flag
            msg = sprintf('\nLipid contamination detected. Applying lipid filter to transient: %d\n', jj);
        elseif ~lipid_flag && water_flag
            msg = sprintf('\nUnstable residual water detected. Applying residual water filter to transient: %d\n', jj);
        elseif lipid_flag && water_flag
            msg = sprintf('\nLipid contamination and unstable residual water detected. Applying lipid and residual water filters to transient: %d\n', jj);
        end
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        DataToAlign(:,jj) = op_Gannet_SignalFilter(spec(:,jj), lipid_flag, water_flag, MRS_struct); %apply filter depending on flags
    end
else
    DataToAlign = double(MRS_struct.fids); % if no flags, proceed
end

time = (0:(MRS_struct.sz(1)-1))'/MRS_struct.spectralwidth;

%% Use similarity index to determine weighting and order

% Use first n points of time-domain data, where n is the last point where abs(diff(mean(SNR))) > 0.5
signal = abs(DataToAlign);
noise = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
SNR = signal ./ repmat(noise, [size(DataToAlign,1) 1]);
SNR = abs(diff(mean(SNR,2)));
SNR = SNR(time <= 0.2);
tMax = find(SNR > 0.5,1,'last');
if isempty(tMax) || tMax < find(time <= 0.1,1,'last')
    tMax = find(time <= 0.1,1,'last');
end

% Flatten complex data for use in spectral registration
flatdata(:,1,:) = real(DataToAlign(1:tMax,:));
flatdata(:,2,:) = imag(DataToAlign(1:tMax,:));

D = zeros(size(flatdata,3)); % Ntotalspec x Ntotalspec
for jj = 1:size(flatdata,3) %Loop over Ntotalspec
    for kk = 1:size(flatdata,3) %Loop over NtotalSpec
        % tmp = sum of squared difference of real part of FID for
        % jj and kk. if jj==kk => D=NaN
        tmp = sum((real(DataToAlign(1:tMax,jj)) - real(DataToAlign(1:tMax,kk))).^2) / tMax;
        if tmp == 0
            D(jj,kk) = NaN;
        else
            D(jj,kk) = tmp;
        end
    end
end
d = nanmedian(D); %Median spectra
[~,alignOrd] = sort(d);

% Set initial reference transient based on similarity index
target = squeeze(flatdata(:,:,alignOrd(1)));
target = target(:);

% Scalar to normalize transients (reduces optimization time)
a = max(abs(target));

% Pre-allocate memory
m = zeros(length(target),size(flatdata,3));

% Starting values for optimization
f0 = Acquire_F0freq2(MRS_struct) * MRS_struct.txfrq*1e-6;
f0 = f0(alignOrd);
f0 = f0 - f0(1);
phi0 = zeros(size(f0));
x0 = [f0(:) phi0(:)];

%% Determine frequency and phase offsets by spectral registration
t = 0:(1/MRS_struct.spectralwidth):(length(target)/2-1)*(1/MRS_struct.spectralwidth);
reverseStr = '';

for jj = alignOrd
    

    msg = sprintf('\nRobust spectral registration - Iteration: %d', jj);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    transient = squeeze(flatdata(:,:,jj));
    fun = @(x) SpecReg(transient(:)/a, target/a, t, x);
    params(jj,:) = lsqnonlin(fun, x0(jj,:), [], [], lsqnonlinopts);

    f   = params(jj,1);
    phi = params(jj,2);
    m_c = complex(flatdata(:,1,jj), flatdata(:,2,jj));
    m_c = m_c .* exp(1i*pi*(t'*f*2+phi/180));
    m(:,jj) = [real(m_c); imag(m_c)];
    resid = target - m(:,jj);
    MSE(jj) = sum(resid.^2) / (length(resid) - 2); % mean square error

    % Update reference
    w(jj) = 0.5*corr(target, m(:,jj)).^2;
    target = (1 - w(jj))*target + w(jj)*m(:,jj);


end
fprintf('\n')
% Copy params to MRS_struct
%MRS_struct.out.RSpecReg.Alligned = 
MRS_struct.out.RSpecReg.freq  = params(:,1);
MRS_struct.out.RSpecReg.phase = params(:,2);
MRS_struct.out.RSpecReg.MSE = MSE;
MRS_struct.out.RSpecReg.water_flag = water_flag;
MRS_struct.out.RSpecReg.lipid_flag = lipid_flag;

% Apply freq & phase to FIDs
for jj = 1:size(flatdata,3)
    %MRS_struct.fids(:,jj) = MRS_struct.fids(:,jj) .* ...
    %    exp(1i*params(jj,1)*2*pi*time) * exp(1i*pi/180*params(jj,2));
    MRS_struct.fids(:,jj) = DataToAlign(:,jj) .* ...
        exp(1i*params(jj,1)*2*pi*time) * exp(1i*pi/180*params(jj,2));
end
MRS_struct.specs = fftshift(ifft(MRS_struct.fids,[],MRS_struct.dims.t),MRS_struct.dims.t);

end%func

function[F0freq2] = Acquire_F0freq2(Struct)
    % Find water frequency
    Struct = op_zeropad(Struct,10);
    Range = find(Struct.ppm<4.68+0.25 & Struct.ppm>4.68-0.25);
    [~,MaxLoc]=max(abs(Struct.specs(Range,:)));
    FR = Struct.ppm(Range);

    F0freq2 = FR(MaxLoc)-4.68;
end