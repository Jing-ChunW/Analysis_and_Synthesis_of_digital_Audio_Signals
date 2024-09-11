% function maskingCurve = spreading(dB, z_masker, bsteps)
% Alignment of levels: dB = 0 if amplitude is 1.
% 
% Reference: Bosi and Goldberg "introduction to digital audio
% coding and standards", pp.185
%
% Yi-Wen Liu, 9/10/2003
% PREVIOUS VERSION HAS A BUG! FIXED 2/3/2004

function maskingCurve = spreading(masker_dB, DELTA, z_masker, bsteps)
  bsteps = bsteps(:);
  applyNonlinearity =0; 
  if masker_dB > -56, applyNonlinearity = 1; 
  end
  dz = bsteps - z_masker;
  theta = (dz>=0);
  maskingCurve = masker_dB- DELTA + abs(dz).*(-27 + 0.37 * ...
		applyNonlinearity*(masker_dB+56)*theta);
  