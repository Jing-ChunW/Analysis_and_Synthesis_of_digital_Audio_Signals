% To replace <MyFindpeaks2020.m> and be called by AsasHwSinMod_2022.m
% May 6, 2022
% Yi-Wen Liu
function [amps, freqs, phs] = mySinAnalysis(X, maxNumPeaks, fNyq, freqadjust)
minSep = 80; % Hz

M = size(X(:),1);
df = fNyq/M;
D = floor(minSep/df); % Minimum separation of peaks by bins

amps = -199*ones(maxNumPeaks,1);
freqs = zeros(maxNumPeaks,1);
phs = zeros(maxNumPeaks,1);
MagdB = 20*log10(abs(X));
MagdB = MagdB(:);
ph_unwrap = unwrap(angle(X));

ind=find( (MagdB(1:M)>[MagdB(1)+100;MagdB((1:M-1))]) ...
    & (MagdB(1:M)>=[MagdB(2:M); MagdB(M)+100]) );   %find the location of peaks

peaks = zeros(maxNumPeaks, 3);
peaks_list = [MagdB(ind), ind]; % Two columns
peaks_sorted = sortrows(peaks_list, 1,'descend'); % sorted according to magnitude
J = length(ind);

if ~isempty(ind)
    % Removing peaks that are too close. Larger ones get the right to
    % remove smaller elements in the list.
    for j = 1:J-1
        peakLoc = peaks_sorted(j,2);
        if peakLoc > 0 % a surviving peak, not removed yet
            for k = j+1:J
                smallerPeakLoc = peaks_sorted(k,2);
                if (smallerPeakLoc > 0) && (abs(peakLoc - smallerPeakLoc) <= D)
                    peaks_sorted(k,2) = 0; % mark the row to be removed
                end
            end
        end
    end
    survivingRows = find(peaks_sorted(:,2) > 1);
    peaks_sorted = peaks_sorted(survivingRows,:);
    J = length(survivingRows);
    for kk = 1:min(maxNumPeaks,J)
        ii = peaks_sorted(kk,2);
        L0 = MagdB(ii);
        A = (MagdB(ii-1)+MagdB(ii+1)-2*L0)/2;
        B = (MagdB(ii+1)-MagdB(ii-1))/2;
        q = -B/(2*A);
        if (ii-1+ freqadjust) > 0
            freqs(kk) = ((ii-1 + freqadjust)+q)*pi/M;
        else
            freqs(kk) = ((ii-1)+q)*pi/M;
        end
        amps(kk) = L0-(B^2)/(4*A);
        if q >= 0
            phs(kk) = (1-q)*ph_unwrap(ii)+q*ph_unwrap(ii+1);
        else
            phs(kk) = (1+q)*ph_unwrap(ii)+(-q)*ph_unwrap(ii-1);
        end
        
    end
    peaks(:,1) = amps;
    peaks(:,2) = freqs;
    peaks(:,3) = phs;
end
%% Return the list in the order of ascending frequency
%peaks = sortrows(peaks,2);	
%amps = peaks(:,1);
%freqs = peaks(:,2);
%phs = peaks(:,3);
return