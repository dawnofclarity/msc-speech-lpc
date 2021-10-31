%% Config
% Singleton Class to hold global configuration
%
classdef Config < handle
    properties (Constant)
        % Edit the following line to the folder containing
        % the sound samples
        SAMPLES_FOLDER = 'Y:\sapr\samples';

        % Folder that holds the results...
        RESULTS_FOLDER = 'Y:\sapr\results';
    end
end
 