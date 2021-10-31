
% clean up any existing windows
close all
clear all

% This is the sample we're going to examine
% fname = 'heed_f.wav'; sampleStart = 1350;
fname = 'heed_f.wav'; sampleStart = 2600;
% fname = 'heed_m.wav'; sampleStart = 2019;

safeFname = replace(fname, '_', '\_');

lpcOrder = 24;

%% Read in a sound sample
[data, sampleRate] = audioread([Config.SAMPLES_FOLDER, '/', fname]);

%% Establish working parameters

% define some well known values
ms1 = sampleRate / 1000;                 % maximum speech Fx at 1000Hz
ms2 = sampleRate / 500;                  % maximum speech Fx at 500Hz
ms20 = sampleRate / 50;                  % minimum speech Fx at 50Hz
colours = 'mcrgy';

fResult = figure;
% set the figure size, making sure the figure titlebar is at least displayed onscreen
tmp = fResult.Position(2);
fResult.Position(2) = (fResult.Position(2) + fResult.Position(4)) - 600;
fResult.Position(3) = 1200;
fResult.Position(4) = 350;
hold on;

legendText = cell(1, 0);
%% Loop for a range of lpcOrder values
for sliceLen = [0.02 0.04 0.06 0.08 0.1 0.12]
    
    % take slice
    samplePeriod = 1/sampleRate;
    sampleCount = floor(sliceLen / samplePeriod);
    sampleEnd = sampleStart + sampleCount - 1;

    %% grap a small slice for autoregressive modelling
    rawSample = data(sampleStart:sampleEnd);

    arSample = rawSample .* hamming(sampleCount);
    
    %% take fft of signal time->signal domain
    fftResult = fft(arSample);

    fftPlot2 = abs(fftResult / sampleCount);    % Take magnitude to yield the double sampling plot (fft yields a symetric plot)
    fftPlot1 = fftPlot2(1:(sampleCount/2)+1);   % Take the first half of the plot to yield a single sampling plot
    fftPlot1 = 20*log10(fftPlot1);              % Scale to dB

    % sampleCount/2  entries normalised 0..1
    f = sampleRate*(0:(sampleCount/2))/sampleCount;

    %% Apply lpc
     [lpcResult, lpcError] = lpc(arSample, lpcOrder);
%     [lpcResult, lpcError] = arcov(arSample, lpcOrder);

    % Derive / plot frequency response (first half)
    [lpcH, lpcW] = freqz(1, lpcResult, sampleCount, sampleRate);

    %% Fundamental freq determination (Autocorrelation)

    % Calculate auto-correlation (max lag => 1 cycle @ 50 Hz)  
    xcorrResult = xcorr(arSample, ms20, 'coeff');   

    % delay axis -1 cycle => +1 cycle @ 50 Hz
    delays = (-ms20:ms20) / sampleRate;          % times of delays

    % Look at region corresponding to positive delays
    % i.e. a window of 1 cycle ahead at 50Hz
    positiveXcorr = xcorrResult(ms20+1:2*ms20+1);

    % Find the max correlation below 500 Hz - xix is index of max value
    [maxXCorr, xix] = max(positiveXcorr(ms2:end));

    % (ms2 + xix - 1) => number of samples between peaks, divide into sampleRate to yield frequency
    xcorrFx = sampleRate / (ms2 + xix -1);

    %% Derive cepstrum (DFT of log spectrum)
    % ift of log of magnitude of ft
    cepstrumResult = ifft(log(abs(fftResult) + eps));

    % derive axis between 1ms (=1000Hz) and 20ms (=50Hz)
    quefrequency = (ms1:ms20) / sampleRate;

    %% Fundamental freq determination (Cepstrum)
    % find the highest peak in the cepstrum between 1khz and 50Hz; cix is index of max value
    [maxCepstrum, cix] = max(abs(cepstrumResult(ms1:ms20)));

    % (ms1 + cix - 1) is offset from start, divide into sampleRate to yield frequency
    cepstrumFx = sampleRate / (ms1 + cix - 1);

    % select the xcorrFx as the one we'll save
    fx = xcorrFx;

    %% Formant Frequencies
    % Find the roots of the polynomial returned by lpc
    arRoots = roots(lpcResult);

    % only keep positive roots; this is ok as they are conjugate pairs.
    arRoots = arRoots(imag(arRoots) > 0);

    % extract the angles for each root
    [frqs,indices] = sort(atan2(imag(arRoots), real(arRoots)) * sampleRate/(2 * pi));

    bandwidth = -1/2 * (sampleRate / (2 * pi)) * log(abs(arRoots(indices)));
    ii = 1;
    formants = zeros(1, 3);
    for jj = 1:length(frqs)
         if (frqs(jj) > 90 && bandwidth(jj) < 150)
             formants(ii) = frqs(jj);
             fprintf(' %d => %0.2f (%0.2f)\n', ii, frqs(jj), bandwidth(jj));
             ii = ii + 1;
             if (ii > numel(formants))
                 break;
             end
         end
    end

    %% write sample to results
    audiowrite([Config.RESULTS_FOLDER, '/', fname], arSample, sampleRate);

%     %% Save data obtained for synthesis
%     fout = sprintf("%s/%s-%d.mat", Config.RESULTS_FOLDER, fname(1:end-4), lpcOrder);
%     save(fout, 'lpcResult', 'fx');
    
    %% Plot LPC curve
        
    plot(lpcW, 20 * log10(abs(lpcH)));  % Scale to dB   
    legendText = [legendText; sprintf('Len %0.2g', sliceLen)];

    %% Log useful stuff to console
    fprintf('File: %s\n', fname);
    fprintf('Sample start: %d, end: %d, count: %d\n', sampleStart, sampleEnd, sampleCount);
    fprintf('LPC: %d, error %0.3g\nLPC terms:', lpcOrder, lpcError);
    fprintf(' %0.3g', lpcResult);
    fprintf('\n');
    fprintf('Cepstrum: max = %g; Fx = %g Hz\n', maxCepstrum, cepstrumFx);
    fprintf('Correlation: max = %g; Fx = %g Hz\n', maxXCorr, xcorrFx);
end
legend(legendText);
xlabel('Frequency (Hz)')
ylabel('Filter Response (dB)')
