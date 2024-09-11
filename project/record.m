%% EE6641 HW: sinusoidal modeling, additive synthesis
% Adapted from a Stanford EE367 lab assignment in 2002.
% Revised Apr 10, 2014
% Revised Feb 21, 2017
% Touched Apr 21, 2020 for the year's HW5.
% Updated May 5, 2020 for S+N decomposition for HW7.
% Updated May 9, 2021 for correct handling of the noise part to allow
%       playback with a change of speed.
% Updated May 7, 2022 to use the new <mySinAnalysis.m> which replaces <MyFindPeaks2020.m> 
%       with improved readability
% Yi-Wen Liu

clear; close all;
DIR = './soundfiles/';
sw.WriteOut = 0;
%fname = 'peaches_16';
DIR_source = '../../english_children/english_words_sentences/11_M_native/port_mic/sentences/';
%DIR_source = '../../english_children/english_free_speech/files_in_one_part/03_F_native/';
%DIR_source = '../../hw3_107011133_¤ý¹t²E/HW3-sounds/';
%FILENAME_source = 'the_dog_is_in_front_of_the_horse';
FILENAME_source = 'the_fish_is_in_the_pond';
%FILENAME_source = 'port';
%fname = 'mymom_16';
fname = 'draw_16';
[x_source,fs1] = audioread([DIR_source FILENAME_source '.wav']);
%sound(x,fs);
fs = 16000;
x_source = resample(x_source,fs,fs1);
x = x_source(:,1);
x = shiftPitch(x, -3);
nx = length(x);

%% ANALYSIS PARAMETERS
frameRate =40;  % frames per second
M = floor(fs/frameRate);  
nFrames = floor(nx/M)*2-1;
R = floor(M/2);  % Note: exact COLA not required
N = 2^(1+floor(log2(5*M+1))); % FFT length, at least a factor of 5 zero-padding

maxPeaks = input('how many peaks to track? '); 
expandRatio = input('time expansion factor?'); % Default = 1.0
freqadjust = 0;
freqShift = 0;
fRatio = 2^(freqShift/12); 

%% VECTOR VARIABLES DECLARATION
fullRangedB = 108.5;    % Depending on actual recording condition. This value is the
                        % derived from calibration of the data in HW5-6.

amps = -fullRangedB * ones(maxPeaks,nFrames);
freqs = zeros(maxPeaks,nFrames);
phases = zeros(maxPeaks,nFrames);

%% ANALYSIS
%
w = blackman(M); % the first window is for analysis
w2 = hann(M+1); w2 = w2(1:end-1); w2 = w2(:); % this second window is for noise synthesis
df = fs/N;
ff = 0:df:(N-1)*df;
fNyq = fs/2;
normfac = sum(w)/2;
x_noise = zeros((nFrames+1)*R,1);
TIQ = hearingThresholdTerhardt(fNyq, N/2) - fullRangedB; % Threshold in Quiet
TIQ = TIQ(:)';
for m=1:nFrames
    tt = (m-1)*R+1:(m-1)*R+M;
    xw = w .* x(tt);
    Xw = fft(xw,N)/normfac; % so full range = 1.0 
    [amps(:,m),freqs(:,m),phases(:,m)] = mySinAnalysis(Xw(1:N/2),maxPeaks,fNyq, freqadjust); 
    %% removing spurious peaks below masking threshold
    numActualPeaks = sum(freqs(:,m)>0);
    freqsHz = freqs(1:numActualPeaks,m)/pi*fNyq;
    maskCurv = calcMaskingCurve(freqsHz, amps(1:numActualPeaks,m), fNyq, N/2);
    maskThres = max([maskCurv; TIQ],[],1);
    toKeep = zeros(numActualPeaks,1);
    for jj = 1:numActualPeaks
        ii = round(freqs(jj,m)/pi*fNyq/df);
        if amps(jj,m) < maskThres(ii)
            freqs(jj,m) = 0; % later this will be an indicator to abandon the sinusoid during synthesis.
        end
    end
    %% show spectrum and peaks
    figure(1);
    subplot(211);
    plot(tt/fs*1000, x(tt));
    xlabel('ms')
    set(gca,'xlim',[(m-1)*R+1 (m-1)*R+M]/fs*1000);
    
    subplot(212);
    Ymag = 20*log10(abs(Xw(1:N/2))); %
    plot(ff(1:N/2),Ymag); hold on;
    plot(freqs(:,m)/(2*pi)*fs, amps(:,m),'x');
    plot(ff(1:N/2)+df, maskThres,'g');
    
    hold off;
    set(gca,'xlim',[0 fNyq]);
    set(gca,'ylim',[-120 0]);
    xlabel('Hz');
    ylabel('dB');
    
    %% Sine + noise decomposition
    sinesum = zeros(M,1);
    for jj = 1:numActualPeaks
        ttsh = 0:M-1; ttsh = ttsh(:);
        if freqs(jj,m) > 0 % a legit peak
            sinesum = sinesum + 10^(amps(jj,m)/20)*cos(freqs(jj,m)*ttsh +phases(jj,m));
        end
    end
    x_noise(tt) = x_noise(tt) + (x(tt) - sinesum).*w2; % using the Hann for synthesis.
    
end

%sound(x_noise,fs);
fprintf('Now playing the noise part after decomposition. Press any key to continue.\n' );
pause();
clear sound;

%% LP coding of the noise part.
sw.emphasis = 0;
p = 8;

LPcoeffs = zeros(p+1,nFrames);
Nfreqs = 2^nextpow2(2*M-1)/2;
numFrames = floor(length(x_noise)/R);
newR = round(R*expandRatio);
noiseResynth = zeros(numFrames* newR, 1);
newM = 2*newR;
hannwin = hann(newM+1);
hannwin = hannwin(1:end-1); hannwin = hannwin(:);

if sw.emphasis == 1
    yemph = filter([1 -0.95],1,x_noise); 
                %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
else
    yemph = x_noise;
end
for kk = 1:nFrames
    fprintf('frame #%d\n',kk);
    ind = (kk-1)*R+1:(kk-1)*R+M;
    indSynth = (kk-1)*newR+1: (kk-1)*newR+newM;
    ywin = yemph(ind);
    %A = lpc(ywin,p); %% This is actually the more direct way to obtain the
    % LP coefficients. But, in this script, we used levinson() instead
    % because it gives us the "reflection coefficients". 
    % 
    Y = fft(ywin,2*Nfreqs);
    PSD = ifft(Y.*conj(Y));
    [A,errvar,K] = levinson(PSD,p); 
    LPcoeffs(:,kk) = A;
    excit_ori = filter(A,1,ywin);
    sigma_ori = sqrt(mean(excit_ori.^2));
    excit = sigma_ori * randn(p+newM,1);
    tmp = filter(1, A, excit);
    noiseResynth(indSynth) = noiseResynth(indSynth) + hannwin.* tmp(p+1:end);
end

%% A synthesizer using MyAdditivesynth2020.m
R = round(R* expandRatio);  % time expansion
y = zeros((nFrames+1)*R,1);
state = zeros(maxPeaks,3);  % [ampInitials, freqInitials, phaseInitials]
freqTrajectories = zeros(maxPeaks, (nFrames+1));
freqTrajTag = zeros(maxPeaks, (nFrames+1));
freqTrajectories(:,1) = 1:length(freqTrajectories(:,1));
for m=2:nFrames-1
    state(:,2) = fRatio * freqs(:,m-1);  % bug fixed May 2020. "m-1" used to be "m".
    tt = (m-1)*R+1: (m+1)*R;  
    [y_synth,phaseUpdate, freqUpdate, freqTag] = MyAdditivesynth2020(amps(:,m),fRatio*freqs(:,m),R,state,fs);
	y(tt)= y(tt) + y_synth;
    state(:,3) = phaseUpdate; % from previous round of synthesis.
    freqTrajectories(:, m) = freqUpdate;
    freqTrajTag(:, m) = freqTag;
end
y_sum = y + noiseResynth;
sound(y_sum,fs);

if sw.WriteOut
    audiowrite(sprintf('%sSine_%d.wav',fname,maxPeaks),y,fs);
    audiowrite(sprintf('%sNoise_%d.wav',fname,maxPeaks),x_noise,fs);
    audiowrite(sprintf('%sNsynth_%d.wav',fname,maxPeaks),noiseResynth,fs);
    audiowrite(sprintf('%sTotSynth_%d.wav',fname,maxPeaks),y+noiseResynth,fs);
end





%% Visualization
figure(10)
param.fs = fs;
S = mySpecgram(x,w,M*3/4,N,param); 
    % mySpecgram() is YWL's preferred way of seeing a spectrogram.
hold on;
yLIM = get(gca,'ylim');
ymax = yLIM(2);
plot((0:nFrames-1)*R/fs,freqs'/pi*ymax,'.','color','g');
set(gcf,'position',[360,80,800,600])

figure(9)
S = mySpecgram(x,w,M*3/4,N,param); 
hold on;
yLIM = get(gca,'ylim');
ymax = yLIM(2);
Ay = [];
for i=1:maxPeaks
    if (freqTrajTag(i, 2) == 1)
        Ay = [Ay, freqs(freqTrajectories(i, 2), 1)];
    end
end
Ay = Ay/pi*ymax;
Ax = ones(1, length(Ay));
Ax = Ax*R/fs;
for i=2:nFrames-1
    By = [];
    for j = 1:maxPeaks
        if (freqTrajTag(j, i) == 1)
            By = [By, freqs(j, i)];
        end
    end
    By = By/pi*ymax;
    Bx = ones(1, length(By));
    Bx = Bx*R/fs*i;
    X = [Ax;Bx];
    Y = [Ay;By];
    line(X, Y, 'color', 'g');
    hold on
    % A
    Ay = [];
    for j = 1:maxPeaks
        if (freqTrajTag(j, i+1) == 1)
            Ay = [Ay, freqs(freqTrajectories(j, i+1), i)];
        end
    end
    Ay = Ay/pi*ymax;
    Ax = ones(1, length(Ay));
    Ax = Ax*R/fs*i;
end
A_born = [];
num_born = zeros(1, nFrames);
for i=1:nFrames-1
    for j = 1:maxPeaks
        if (freqTrajTag(j, i) == 2)
            if (freqs(j, i)/pi*ymax  >= 4)
                num_born(i) = num_born(i) + 1;
            end
        end
    end
end
num_born_threshold = 3
range_born = 10
born_line_tag = zeros(1, nFrames);
for i=1:nFrames-1
    x_noline = 0;
    if num_born(i) > num_born_threshold
        for j = 1:range_born
            if i - j > 0
                if (born_line_tag(i-j) == 1)
                    x_noline = 1;
                end
            end
        end
        if x_noline == 0
            born_line_tag(i) = 1;
            xline(R/fs*i);
            hold on
        end
    end
end

for i=1:nFrames-1
    A_born_line = [];
    for j = 1:maxPeaks
        if (freqTrajTag(j, i) == 2)
            if (freqs(j, i)  >= 4)
            A_born_line = [A_born_line, freqs(j, i)];
            end
        end
    end
    A_born_line = [A_born_line, zeros(1, (max(num_born)) - length(A_born_line))];
    A_born = [A_born; A_born_line];
end
pre_born_line = 0;
cut_line = [];
for n=1:nFrames
    if n == 1
        if (born_line_tag(n) == 1)
            pre_born_line = n;
        end
    else
        if (born_line_tag(n) == 1)
            for point = 1:round((n - pre_born_line)*1/50)
                cut_line = [cut_line, round((pre_born_line + n)/2) - round((n - pre_born_line)*1/6) - 1 + point];
            end
            cut_line = [cut_line, round((pre_born_line + n)/2)];
            for point = 1:round((n - pre_born_line)*1/50)
                cut_line = [cut_line, round((pre_born_line + n)/2) + point];
            end
            pre_born_line = n;
        end
    end
end
            

tmp = 0;
y_record_breakword = [];
for m=1:nFrames-1
    tt = (m-1)*R+1: (m+1)*R;
    for i=1:length(cut_line)
        if m == cut_line(i)
            if i == 1
                y_new = [y_sum(1:(m-1)*R)];
                y_record_breakword = [zeros(length(tt)*10,1); y_sum(1:(m-1)*R)];
                tmp = (m+1)*R+1;
            else
                if (tmp < (m-1)*R)
                    if (length(y_new) > R)
                        y_new(end-R:end) = y_new(end-R:end) + y_sum(tmp:tmp + R);
                    end
                    y_new = [y_new; y_sum(tmp+R+1:(m-1)*R)];
                    y_record_breakword = [y_record_breakword; zeros(length(tt)*10,1); y_sum(tmp:(m-1)*R)];
                end
                tmp = (m+1)*R+1;
            end
        end
    end
            
end
for i=1:nFrames-1
    % A
    Ay = [];
    
    for j = 1:maxPeaks
        if (freqTrajTag(j, i) == 0)
            Ay = [Ay, freqs(j, i)];
        end
    end
    Ay = Ay/pi*ymax;
    Ax = ones(1, length(Ay));
    Ax = Ax*R/fs*i;
    scatter(Ax, Ay, [],  [0 1 0], 'marker', 'd'); %Green
    hold on
    Ay = [];
    for j = 1:maxPeaks
        if (freqTrajTag(j, i) == 2)
            Ay = [Ay, freqs(j, i)];
        end
    end
    Ay = Ay/pi*ymax;
    Ax = ones(1, length(Ay));
    Ax = Ax*R/fs*i;
    
    scatter(Ax, Ay, [],  [1 0 0], 'marker', 'd'); %Red
    hold on
end


figure(3)
S = mySpecgram(x_noise,w,M*3/4,N,param); 
set(gcf,'position',[400,120,800,600]); title('noise spectrogram');

figure(4)
plot((0:length(x)-1)/fs, x); hold on;
plot((0:length(x_noise)-1)/fs, x_noise);

figure(5)
subplot(211); plot((0:length(x_noise)-1)/fs, x_noise);
xlabel('time (s)'); ylabel('the noise part');
subplot(212); plot((0:length(noiseResynth)-1)/fs, noiseResynth);
xlabel('time (s)'); ylabel('synthesized noise part');


setFontSizeForAll(14);


