%% Synthesise a speech sample
% Dawn Raison 6609229
% - Create an impulse train
% - filter it against LPC to produce new sound
% - Display plots for documentation

% clean up any existing windows
close all
clear all

%% Load analysis data for this sample / order
fname = 'heed_f.wav';
lpcOrder = 24;
sampleRate = 24000;

% load LPC params calculated in analyse.m
fin = sprintf("%s/%s-%d.mat", Config.RESULTS_FOLDER, fname(1:end-4), lpcOrder);
load(fin, 'lpcResult', 'fx');

% visual check the correct LPC values were loaded.
fprintf('LPC terms:');
fprintf(' %0.3f', lpcResult);
fprintf('\n');

%% Make an impulse train
% time vector 0 - 1 second
tv=0:1/sampleRate:1;
slotCount = numel(tv);

% empty inpluse train
impulseTrain=zeros([1 slotCount]);

% Set a periodic 1 bit to represent the fundamental frequency
impulseTrain(1:sampleRate/fx:end)=1;

%% Filter the impulse train using the LPC
shapedTrain = filter(1, lpcResult, impulseTrain)';

%% Store generated sound

% normalise 1:-1
normalisedShape = shapedTrain / (max(abs(shapedTrain)));
% write to wav file.
audiowrite(sprintf("%s/%s-%d.wav", Config.RESULTS_FOLDER, ...
    fname(1:end-4), lpcOrder), normalisedShape, sampleRate);

%% Display plots

% display the first 2401 slots (100ms) to ensure the plots are legible.
dispSlots = 600;

fResult = figure;
tlResult = tiledlayout(fResult, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
% set the figure size, making sure the figure titlebar
% is at least displayed onscreen
tmp = fResult.Position(2);
fResult.Position(2) = (fResult.Position(2) + fResult.Position(4)) - 600;
fResult.Position(3) = 1200;
fResult.Position(4) = 350;

% Display the start of the raw impluse train
nexttile;
tvDisp = tv(1:dispSlots);
impulseTrainDisp = impulseTrain(1:dispSlots);
plot(tvDisp,impulseTrainDisp);
xlabel('Sample');
ylabel('Amplitude');
axis off;

% Display the start of the created sound
nexttile;
normalisedShapeDisp = normalisedShape(1:dispSlots);
plot(tvDisp,normalisedShapeDisp);
xlabel('Sample');
ylabel('Amplitude');
axis off;

