function out = op_Gannet_SignalFilter(in, lipid_flag, water_flag, MRS_struct)
% adapted from Gannet 3.1
% Based on:
% Golotvin & Williams, 2000. Improved baseline recognition
%   and modeling of FT NMR spectra. J. Magn. Reson. 146, 122-125
% Cobas et al., 2006. A new general-purpose fully automatic
%   baseline-correction procedure for 1D and 2D NMR data. J. Magn. Reson.
%   183, 145-151

freqRange = MRS_struct.spectralwidth/(MRS_struct.txfrq*1e-6);


freq = (MRS_struct.sz(1) + 1 - (1:MRS_struct.sz(1))) / MRS_struct.sz(1) * freqRange + 4.68 - freqRange/2;

z = BaselineModeling(in, lipid_flag, water_flag, freqRange);

waterLim = freq <= 5.5 & freq >= 3.6;
noiseLim = freq <= 11 & freq >= 10;

y.real = real(in) - z.real;
y.imag = imag(in) - z.imag;
if water_flag % some residual water may remain so replace with noise
    %CWJ
    %rng('shuffle');
    noise.real = datasample(y.real(noiseLim), sum(waterLim), 'Replace', true);
    noise.imag = datasample(y.imag(noiseLim), sum(waterLim), 'Replace', true);
    y.real(waterLim) = noise.real;
    y.imag(waterLim) = noise.imag;
end
%cwj swapped to fft, in line with FID-A
out = fft(fftshift(complex(y.real, y.imag)));

end


function z = BaselineModeling(y, lipid_flag, water_flag, freqRange)

% Power spectrum of first-derivative of signal calculated by CWT
Wy = abs(cwt2(real(y), 10)).^2;

freq = (length(Wy) + 1 - (1:length(Wy))) / length(Wy) * freqRange + 4.68 - freqRange/2;

noiseLim = freq <= 11 & freq >= 10;

sigma = std(Wy(noiseLim));

w = 1:5;
k = 3;
baseline = zeros(length(Wy),1);
signal = zeros(length(Wy),1);

while 1
    if w(end) > length(Wy)
        break
    end
    h = max(Wy(w)) - min(Wy(w));
    if h < k*sigma
        baseline(w) = Wy(w);
    else
        signal(w) = Wy(w);
    end
    w = w + 1;
end

% Include lipids and water in baseline estimate, as appropriate
lipidLim = freq <= 1.85 & freq >= -2;
waterLim = find(freq <= 5.5 & freq >= 3.6);
if lipid_flag
    baseline(lipidLim) = Wy(lipidLim);
end
if water_flag
    baseline(waterLim) = Wy(waterLim);
end

z.real = real(y);
z.real(baseline == 0) = 0;

if lipid_flag
    z_lipids = whittaker(z.real(lipidLim), 2, 10);
end
if water_flag
    z_water = whittaker(z.real(waterLim), 2, 0.2);
end
z.real = whittaker(z.real, 2, 1e3);
if lipid_flag
    z.real(lipidLim) = z_lipids;
end
if water_flag
    z.real(waterLim) = z_water;
end

% cwj: Removed minus sign here due to difference in FID-A Fourier rep
z.imag = imag(hilbert(z.real));

end


function Wy = cwt2(y, a)

precis = 10;
coef = sqrt(2)^precis;
pas  = 1/2^precis;

lo    = [sqrt(2)*0.5 sqrt(2)*0.5];
hi    = lo .* [1 -1];
nbpts = (length(lo)-1)/pas+2;

psi = coef*upcoef2(lo,hi,precis);
psi = [0 psi 0];
x = linspace(0,(nbpts-1)*pas,nbpts);

step = x(2) - x(1);
wav = cumsum(psi)*step;
x = x - x(1);

y = y(:)';
j = 1+floor((0:a*x(end))/(a*step));
if length(j) == 1
    j = [1 1];
end
f = fliplr(wav(j));

Wy = -sqrt(a) * keepVec(diff(conv2(y,f,'full')),length(y));

    function out = upcoef2(lo,hi,precis)
        out = hi;
        for k = 2:precis
            out = conv2(dyadup2(out),lo,'full');
        end
        function out = dyadup2(in)
            z = zeros(1,length(in));
            out = [in; z];
            out(end) = [];
            out = out(:)';
        end
    end

    function out = keepVec(in,len)
        out = in;
        ok = len >= 0 && len < length(in);
        if ~ok
            return
        end
        d     = (length(in)-len)/2;
        first = 1+floor(d);
        last  = length(in)-ceil(d);
        out   = in(first:last);
    end

end


function z = whittaker(y, d, lambda)

% Code taken from: Eilers, 2003. A perfect smoother. Anal. Chem. 75,
% 3631-3636

y = y(:);

m = length(y);
E = speye(m);
D = diff(E,d);

w = double(y ~= 0);
W = spdiags(w,0,m,m);
C = chol(W + lambda*(D'*D));
z = C\(C'\(w.*y));

end



