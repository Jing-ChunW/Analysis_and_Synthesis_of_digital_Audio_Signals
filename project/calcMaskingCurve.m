% function curve = calcMaskingCurve(freqs, dBs, fNyq, M)
% freqs: list of tonal peak frequencies in Hz
% dBs: list of their amplitudes in dB (assuming 0 dB = full range)
% fNyq: Nyquist frequency
% M: FFT length divided by 2
function curve = calcMaskingCurve(freqs, dBs, fNyq, M)
J = length(freqs);
df = fNyq/M;
ff = df:df:fNyq;
bb = freq2bark(ff);
DELTA = 13.0;
spreadfuns = zeros(J,M); % storage place for individual spreading functions
z_masker = freq2bark(freqs);

for jj = 1:J
    spreadfuns(jj,:) = spreading(dBs(jj), DELTA, z_masker(jj), bb);
end
curve = max(spreadfuns,[],1);
return
