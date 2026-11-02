%% loading in LSL file and preprocessing the data
clc; clear all; close all;
% set relevant paths %
gitPath = 'C:\Users\yoelgo\Documents\GitHub\'; % path of the local GitHub
cd([gitPath 'MP_analyses\MirrorPods_ExpAnalysis\']); % path of this specific project in Github
addpath([gitPath 'MP_analyses\extra_functions\']); % add the function in the "extra_functions" folder in the github project
env = setupEnviroment11(); % a custom function which sets all the paths and prepares the variable enviroment (env) for this experiment
clear gitPath

%%
n = 2;
dat = load(fullfile(env.data.cleanFiles(n).folder,...
    env.data.cleanFiles(n).name), 'dat_after_ICA').dat_after_ICA;
ID = extractBefore(env.data.cleanFiles(n).name, '_EEG');
env.ID = ID;
%% ===================== settings =====================
trlRow     = 19;                        % row of trlTable1 to analyse
winLen     = 10;                         % window length (s)
overlapPct = 90;                        % window overlap (%)
channels   = {'FC1','FC2'};
Hz         = 40;                        % frequency for the time course
xfreq      = [35 45];                   % FFT x-axis range
col        = [0.85 0.33 0.10];
chStr      = strjoin(channels, ', ');

cfg = [];
cfg.trial = trlRow;
trlDat = ft_selectdata(cfg, dat)

% divide into windows, drop windows with NaNs
fs       = dat.fsample;
wlen     = round(winLen * fs);                              % window length in samples
stepSamp = max(round(wlen * (1 - overlapPct/100)), 1);      % step between window starts
X        = trlDat.trial{1};
starts   = 1:stepSamp:(size(X,2) - wlen + 1);
tWin     = (starts - 1 + wlen/2) / fs;                      % centre time of each window (s from trial start)

winDat = [];
winDat.label   = trlDat.label;
winDat.fsample = fs;
winDat.trial   = {};
winDat.time    = {};
keep = false(size(starts));
for w = 1:numel(starts)
    seg = X(:, starts(w):starts(w)+wlen-1);
    if any(isnan(seg(:)))
        continue
    end
    keep(w) = true;
    winDat.trial{end+1} = seg;
    winDat.time{end+1}  = (0:wlen-1) / fs;
end
fprintf('%d of %d windows kept\n', sum(keep), numel(keep));
if ~any(keep)
    error('All windows of this trial contain NaNs.');
end

% FFT per window
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'pow';
cfg.taper      = 'hanning';
cfg.foi        = 2:1/winLen:45;
cfg.pad        = 'nextpow2';
cfg.keeptrials = 'yes';
cfg.channel    = channels;
freqW = ft_freqanalysis(cfg, winDat);

P = mean(freqW.powspctrm, 2);                               % avg over channels
P = reshape(P, size(P,1), []);                              % windows x freq
f = freqW.freq(:)';
[~, fi] = min(abs(f - Hz));                                 % index of the 40 Hz bin

p40       = NaN(size(starts));                              % NaN = window removed
p40(keep) = P(:, fi);

% ===================== figure =====================
figure('Color','w','Position',[100 100 1400 550]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('ID = %s, Participant 2  |  trail condition  |  channels: %s', ID, chStr), ...
      'FontSize', 16, 'FontWeight','bold', 'Interpreter','none');

% --- left: average FFT over windows
axL = nexttile; hold(axL,'on');
plotMeanSEM(axL, f, P, col);
xline(axL, Hz, '--k', 'LineWidth', 1);
xlim(axL, xfreq);
xlabel(axL, 'Frequency (Hz)'); ylabel(axL, 'Power (\muV^2)');
title(axL, sprintf('FFT (average of %d windows)', sum(keep)), 'FontSize', 14);

% --- right: 40 Hz value of each window over time
axR = nexttile; hold(axR,'on');
plot(axR, tWin, p40, 'Color', col, 'LineWidth', 1.2);
xlabel(axR, 'Time (s)'); ylabel(axR, sprintf('Power at %g Hz (\\muV^2)', Hz));
title(axR, sprintf('%g Hz power per %g s window (%g%% overlap)', Hz, winLen, overlapPct), 'FontSize', 14);

set([axL axR], 'FontSize', 12, 'LineWidth', 1, 'TickDir','out', 'Box','off');
linkaxes([axL axR], 'y');
% same exponent notation on both y axes
e = floor(log10(max([axL.YLim(2) axR.YLim(2)])));
axL.YAxis.Exponent = e;
axR.YAxis.Exponent = e;

fname = sprintf('nonlinked_%s_EEG2_%gsecWindow_FFT_%gHz_trail.jpg', ID, winLen, Hz);
saveas(gcf, fullfile("C:\Users\yoelgo\Documents\GitHub\MP_analyses\MirrorPods_ExpAnalysis\", fname));
%% ===================== local functions (must stay at END of the script) =====================
function h = plotMeanSEM(ax, f, X, col)
    m = mean(X, 1, 'omitnan');
    s = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X(:,1))));
    fill(ax, [f fliplr(f)], [m+s fliplr(m-s)], col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h = plot(ax, f, m, 'Color', col, 'LineWidth', 2);
end