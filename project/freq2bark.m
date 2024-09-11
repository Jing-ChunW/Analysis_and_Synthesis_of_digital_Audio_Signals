% freqs is a frequency vector in Hz, z is in Barks.

function z= freq2bark(freqs)
  f_kHz = freqs/1000;
  z = 13*atan(0.76*f_kHz) + 3.5*atan((f_kHz/7.5).^2);