function [ft_data] = xdf2ft(EEG_stream, lay)
% Create a fieldtrip suitable data structure from the  
% EEG stream of the LSL (xdf file).
% input: EEG stream from the LSL
% output: fieldtrip compatible structure for ft_preprocessing.

    dataRaw = EEG_stream.time_series; % EEG signals
    timeEEG = EEG_stream.time_stamps; % Corresponding time stamps
    
    % Extract channel labels from the structs
    ft_data.label = cellfun(@(x) x.label, EEG_stream.info.desc.channels.channel, 'UniformOutput', false);
    lay.label(ismember(lay.label, {'AFz', 'TP9'})) = [];
    ft_data.label(1:64) = lay.label(1:64);

    % Continue with the rest of the FieldTrip structure preparation
    ft_data.fsample = str2double(EEG_stream.info.nominal_srate); % Sampling rate
    ft_data.time = {timeEEG - timeEEG(1)}; % Time in seconds (start from zero)
    ft_data.trial = {dataRaw}; % EEG data as a cell array
    
    % Calculate the number of samples
    numSamples = size(dataRaw, 2); % Assuming data is channels x samples
    
    % Add the sampleinfo field to ft_data
    ft_data.sampleinfo = [1, numSamples]; % Start at sample 1 and end at the last sample
    
    % FieldTrip expects time and trial to be cells (even for continuous data)
    ft_data.time = {timeEEG - timeEEG(1)}; % Time in seconds
    ft_data.trial = {dataRaw};          % EEG data as a cell array


    if any(strcmp(ft_data.label, 'Markers'))
        idx = find(strcmp(ft_data.label, 'Markers'));
        ft_data.trial{1}(idx,:) = [];
        ft_data.label(idx) = [];
    end

end