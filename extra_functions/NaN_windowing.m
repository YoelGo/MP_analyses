function [data_win] = NaN_windowing(cfg, data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split data into windows and expluce any window that contains NaN for
% fft or LAVI purposes.
%
% variables:
% cfg.winlen = length of each divided window (sec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====== split each trial into windows ======
fsample  = data.fsample;      % sampling rate
winSamp  = round(cfg.winlen*fsample);

trl_all = [];

begSamp = all_table.beg_sample(row);
endSamp = all_table.end_sample(row);
starts = begSamp:winSamp:endSamp-winSamp;
for st = starts
    en = st+winSamp-1;
    trl_all = [trl_all; st en 0];
end

if isempty(trl_all)
    warning('No 5s windows for subj %s, bps=%d, loc=%d', subjID, bpsid, locv);
    % log empty row
    newRow = {subjID, bpsid, locv, 0, 0, 0};
    windowLog = [windowLog; newRow];
end

total_windows = size(trl_all,1);

cfg = [];
cfg.trl = trl_all;
data_win = ft_redefinetrial(cfg, data);

% ====== exclude windows containing NaNs ======
keep = true(1,numel(data_win.trial));
for tr = 1:numel(data_win.trial)
    if any(isnan(data_win.trial{tr}(:)))
        keep(tr) = false;
    end
end
removed_windows = sum(~keep);
remaining_windows = total_windows - removed_windows;

data_win.trial = data_win.trial(keep);
data_win.time  = data_win.time(keep);

% log into table
newRow = {subjID, bpsid, locv, total_windows, removed_windows, remaining_windows};
windowLog = [windowLog; newRow];

if isempty(data_win.trial)
    warning('All windows contained NaNs for subj %s, bps=%d, loc=%d', subjID, bpsid, locv);
end
end