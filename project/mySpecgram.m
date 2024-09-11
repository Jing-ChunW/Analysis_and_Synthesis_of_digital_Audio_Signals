%% function S = mySpecgram(x,win,Novlp,Nfft,param)
% Wrote this customized function for 2017 ASAS usages because the new
% default spectrogram() does not give me the plot in ways I like.
% 2017/12/21
% 
% x, win, Novlp, Nfft: same as definition in spectrogram()
% param is a structure that contains fields such as fs
function S = mySpecgram(x,win,Novlp,Nfft,param)
S = spectrogram(x,win,Novlp,Nfft);
Df = param.fs/ Nfft;
fsteps = 0:Df:param.fs/2;
Dt = Nfft/2/param.fs;
tsteps = 0:Dt:length(x)/param.fs;
cm = colormap('gray');
colormap(1-cm);
imagesc(tsteps,fsteps/1000,20*log10(abs(S)));
%imagesc(tsteps,fsteps/1000,abs(S));

set(gca,'ydir','normal');
%set(gca,'ytick',0:2:22);
ylabel('frequency (kHz)')
xlabel('time (sec)')

return