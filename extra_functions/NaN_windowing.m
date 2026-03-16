function [data_win, windowSummary] = NaN_windowing(cfg, data)
% Split FieldTrip data into windows and remove NaN windows
%
% cfg.winlen  - window length (sec)
% cfg.overlap - percent overlap (0–99). If 0 or missing → no overlap.

fs   = data.fsample;
wlen = round(cfg.winlen * fs);

% -------- overlap handling --------
if ~isfield(cfg,'overlap') || cfg.overlap == 0
    stepSamp = wlen;   % exact non-overlapping windows
else
    if cfg.overlap < 0 || cfg.overlap >= 100
        error('cfg.overlap must be between 0 and <100');
    end
    stepSamp = round(wlen * (1 - cfg.overlap/100));
    stepSamp = max(stepSamp,1); % safety
end
% ----------------------------------

data_win = [];
data_win.label   = data.label;
data_win.fsample = fs;
data_win.trial   = {};
data_win.time    = {};
data_win.sampleinfo = [];
data_win.trialinfo  = [];

nTrials = numel(data.trial);
totalW = zeros(nTrials,1);
keptW  = zeros(nTrials,1);

winCount = 0;

for t = 1:nTrials

    trialData = data.trial{t};
    nSamples  = size(trialData,2);

    starts = 1:stepSamp:(nSamples - wlen + 1);
    totalW(t) = numel(starts);

    for w = 1:numel(starts)

        beg = starts(w);
        fin = beg + wlen - 1;

        segment = trialData(:,beg:fin);

        % reject NaN windows
        if any(isnan(segment(:)))
            continue
        end

        keptW(t) = keptW(t) + 1;
        winCount = winCount + 1;
        data_win.trial{winCount} = segment;
        data_win.time{winCount}  = (0:wlen-1)/fs;
        data_win.sampleinfo(winCount,:) = [beg fin];
        data_win.trialinfo(winCount,:)  = [t w];
    end 
end

% ----- summary table -----
removedW = totalW - keptW;
percentRemoved = 100 * removedW ./ max(totalW,1);

windowSummary = table( ...
    (1:nTrials)', totalW, keptW, removedW, percentRemoved, ...
    'VariableNames', {'trial','total_windows','kept_windows','removed_windows','percent_removed'});

if isempty(data_win.trial)
    warning('No valid windows remained after NaN rejection.');
end
display(windowSummary)
end
