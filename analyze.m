%% Analyse a speech sample
% Dawn Raison 6609229
% - Generate LPC coefficients
% - Determine First 3 formant frequencies
% - Determine fundamental frequency
% - Create plots for documentation

% clean up any existing windows
close all
clear all

% This is the sample we're going to examine
fname = 'heed_f.wav'; sampleStart = 2600;
% fname = 'heed_m.wav'; sampleStart = 2019;

% allow filename to be used in title elements, etc.
safeFname = replace(fname, '_', '\_');

%% Read in a sound sample
[data, sampleRate] = audioread([Config.SAMPLES_FOLDER, '/', fname]);

% define some well known values
ms1 = sampleRate / 1000;                 % maximum speech Fx at 1000Hz
ms2 = sampleRate / 500;                  % maximum speech Fx at 500Hz
ms20 = sampleRate / 50;                  % minimum speech Fx at 50Hz

%% Establish working parameters
lpcOrder = 24;

colours = 'mcrgy';

%% grap a small slice for autoregressive modelling

% 100 mSec sample
samplePeriod = 1/sampleRate;
sampleCount = floor(0.1 / samplePeriod);
sampleEnd = sampleStart + sampleCount - 1;

rawSample = data(sampleStart:sampleEnd);

% apply window to the extracted sample
arSample = rawSample .* hamming(sampleCount);

%% take fft of signal time->signal domain
fftResult = fft(arSample);

% Take magnitude to yield the double sampling plot (symetric plot)
fftPlot2 = abs(fftResult / sampleCount);
% Take the first half of the plot to yield a single sampling plot
fftPlot1 = fftPlot2(1:(sampleCount/2)+1);
% Scale to dB
fftPlot1 = 20*log10(fftPlot1);              

% sampleCount/2  entries normalised 0..1
f = sampleRate*(0:(sampleCount/2))/sampleCount;

%% Apply lpc
[lpcResult, lpcError] = lpc(arSample, lpcOrder);
% alternative is to use arcov
% [lpcResult, lpcError] = arcov(arSample, lpcOrder);

% Derive / plot frequency response (first half)
[lpcH, lpcW] = freqz(1, lpcResult, sampleCount, sampleRate);

%% Derive residual by applying an inverse filter
% use whole plot to match full FFT.
[hh, ww] = freqz(1, lpcResult, sampleCount, 'whole');
freqResidual = fftResult./hh;
timeResidual = real(ifft(freqResidual));

%% Fundamental freq determination (Autocorrelation)

% Calculate auto-correlation (max lag => 1 cycle @ 50 Hz)  
xcorrResult = xcorr(arSample, ms20, 'coeff');   

% delay axis -1 cycle => +1 cycle @ 50 Hz
delays = (-ms20:ms20) / sampleRate;          % times of delays

% Look at region corresponding to positive delays
% i.e. a window of 1 cycle ahead at 50Hz
positiveXcorr = xcorrResult(ms20+1:2*ms20+1);

% Find the max correlation below 500 Hz
% xix is index of max value
[maxXCorr, xix] = max(positiveXcorr(ms2:end));

% (ms2 + xix - 1) => number of samples between peaks
% divide into sampleRate to yield frequency
xcorrFx = sampleRate / (ms2 + xix -1);

%% Derive cepstrum (DFT of log spectrum)
% ift of log of magnitude of ft
cepstrumResult = ifft(log(abs(fftResult) + eps));

% derive axis between 1ms (=1000Hz) and 20ms (=50Hz)
quefrequency = (ms1:ms20) / sampleRate;

%% Fundamental freq determination (Cepstrum)
% find the highest peak in the cepstrum between 1khz and 50Hz
% cix is index of max value
[maxCepstrum, cix] = max(abs(cepstrumResult(ms1:ms20)));

% (ms1 + cix - 1) is offset from start
% divide into sampleRate to yield frequency
cepstrumFx = sampleRate / (ms1 + cix - 1);

% select the xcorrFx as the one we'll save
fx = xcorrFx;

%% Formant Frequencies
% Find the roots of the polynomial returned by lpc
% See https://uk.mathworks.com/help/signal/ug/formant-estimation-with-lpc-coefficients.html
arRoots = roots(lpcResult);

% only keep positive roots; this is ok as they are conjugate pairs.
arRoots = arRoots(imag(arRoots) > 0);

% extract the angles for each root
[frqs,indices] = sort(atan2(imag(arRoots), real(arRoots)) ...
    * sampleRate/(2 * pi));

bandwidth = -1/2 * (sampleRate / (2 * pi)) ...
    * log(abs(arRoots(indices)));

ii = 1;
formants = zeros(1, 3);
for jj = 1:length(frqs)
     if (frqs(jj) > 90 && bandwidth(jj) < 150)
         formants(ii) = frqs(jj);
         fprintf(' %d => %0.2f (%0.2f)\n', ...
             ii, frqs(jj), bandwidth(jj));
         ii = ii + 1;
         if (ii > numel(formants))
             break;
         end
     end
end

%% write sample to results
audiowrite([Config.RESULTS_FOLDER, '/', fname], ...
    arSample, sampleRate);

%% Save data obtained for synthesis
fout = sprintf("%s/%s-%d.mat", ...
    Config.RESULTS_FOLDER, fname(1:end-4), lpcOrder);
save(fout, 'lpcResult', 'fx');

%% Log useful stuff to console
fprintf('File: %s\n', fname);
fprintf('Sample start: %d, end: %d, count: %d\n', ...
    sampleStart, sampleEnd, sampleCount);
fprintf('LPC error %0.3f\nLPC terms:', lpcError);
fprintf(' %0.3f', lpcResult);
fprintf('\n');
fprintf('Cepstrum: max = %g; Fx = %g Hz\n', maxCepstrum, cepstrumFx);
fprintf('Correlation: max = %g; Fx = %g Hz\n', maxXCorr, xcorrFx);

%% Display
CANVAS_ROWS=4;
CANVAS_COLS=3;

% create figure
fResult = figure();

% set the figure size, making sure the figure titlebar
% is at least displayed onscreen
tmp = fResult.Position(2);
fResult.Position(2) = (fResult.Position(2) + fResult.Position(4)) - 600;
fResult.Position(3) = 1024;
fResult.Position(4) = 600;

% create a tiled layout for the figure
tlResult = tiledlayout(fResult, CANVAS_ROWS, CANVAS_COLS, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
tlResult.Title.String = ...
    sprintf('Analysis for %d to %d of sample %s; order %d', ...
    sampleStart, sampleStart + sampleCount - 1, safeFname, lpcOrder);

% the amplitude plot of the whole sample
nexttile(tlResult, 1, [1 3]);
plot(data);
legend('Amplitude');
title( 'Amplitude (whole sample)');
xline(sampleStart, 'LineWidth', 2, 'LineStyle', '-.', ...
    'DisplayName', 'Sample Start', 'Color', colours(1));
xline(sampleEnd, 'LineWidth', 2, 'LineStyle', '-.', ...
    'DisplayName', 'Sample End', 'Color', colours(2));
xlabel('Sample');
ylabel('Amplitude');
axis tight;

% Plot raw sample
nexttile(tlResult, 4);
plot(rawSample);
title( 'Amplitude (slice)');
xlabel('Sample');
ylabel('Amplitude');
legend('Original Sample');
axis tight;

% Plot windowed sample
nexttile(tlResult, 7);
plot(arSample);
title( 'Amplitude (windowed slice)');
xlabel('Sample');
ylabel('Amplitude');
legend('Windowed sample');
axis tight;

% plot residual
nexttile(tlResult, 10);
plot(timeResidual);
title( 'Residual');
legend('Residual');
axis tight;

% plot spectogram
nexttile(tlResult, 6);
segmentlen = 100;
noverlap = 90;
fftSize = 256;
spectrogram(arSample, segmentlen, noverlap, fftSize, sampleRate, 'yaxis');
title('Signal Spectrogram');
axis tight;

% plot amplitude spectrum
nexttile(tlResult, 5);
plot(f, fftPlot1);
title('Amplitude Spectrum f(t)');
xlabel('Frequency (Hz)')
ylabel('Ampliture (dB)')
legend('Speech spectrum');
axis tight;

% lpc filter response
nexttile(tlResult, 8);
plot(lpcW, 20 * log10(abs(lpcH)));  % Scale to dB
title('Freq. response of filter');
xlabel('Frequency (Hz)')
ylabel('Filter Response (dB)')
legend('lpc spectrum');
axis tight;

for ii = 1:numel(formants)
    x = formants(ii);
    xline(x, 'LineWidth', 1, 'LineStyle', '-.', 'DisplayName', ...
        sprintf('Formant %d (%.1f Hz)', ii, x), 'Color', colours(ii));
end

% cepstrum
nexttile(tlResult, 12);
plot(quefrequency, abs(cepstrumResult(ms1:ms20)));
title( 'Cepstrum');
legend('Cepstrum');
xlabel('Quefrency (s)');
ylabel('Amplitude');
axis tight;
xline(quefrequency(cix), 'LineWidth', 1, 'LineStyle', '-.', ...
    'DisplayName', 'Max', 'Color', 'r');

% autocorrelation
nexttile(tlResult, 9);
plot(delays, xcorrResult);
title('Autocorrelation');
legend('Autocorrelation');
xlabel('Delay (s)');
ylabel('Correlation coeff.');
axis tight;
xline(delays(ms20 + 1), 'LineWidth', 1, 'LineStyle', '-.', ...
    'DisplayName', 'Centre', 'Color', 'g');
xline(delays(ms20 + ms2 + xix), 'LineWidth', 1, 'LineStyle', '-.', ...
    'DisplayName', 'Max', 'Color', 'r');
