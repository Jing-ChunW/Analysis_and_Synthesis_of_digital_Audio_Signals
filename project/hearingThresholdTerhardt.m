function tiq = hearingThresholdTerhardt(NyquistRate, M)
  df = NyquistRate/M;  
  fsteps = df:df:NyquistRate;  
  fsteps = fsteps(:);
  f_kHz = fsteps/1000;
  tiq = 3.64* f_kHz.^(-0.8) -6.5*exp(-0.6*(f_kHz -3.3).^2) + 0.001* ...
	f_kHz.^4;
  
  %tiq will be column vector
  