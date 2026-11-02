%% loading in LSL file and preprocessing the data
clc; clear all; close all;
rootDir = 'C:\Users\yoelgo\Documents\GitHub\MP_analyses\MirorPods_FT_Pilot\'; 
addpath(genpath(rootDir));
cd(rootDir)
env = setupEnviroment11();

conds = {'4bps drumming over 40Hz', '(2)4bps O2 -40db',...    % 1 2
         '(3)4bps O6 -30db', '(4)4bps O6 -40db',...    % 3 4
         '(5)4bps O19 -30db', '(6)4bps O19 -40db',...  % 5 6
         '(7)4bps O20 -30db(30Hz)', '(8)4bps O20 -40db(30Hz)',...   % 7 8
         'Barbara replication, 39Hz 41Hz tagging'...                  % 9
         '2Hz beep over 40Hz', 'only 40Hz',...  % 10 11
         '(12)Beep O6 -20db', '(13)noBeep O6 -20db',...  % 12 13
         '(14)Beep O19 -25db', '(15)noBeep O19 -25db',...   % 14 15
         '(16)Beep O20 -25db(30Hz)', '(17)noBeep O20 -25db(30Hz)'}      % 16 17
%% ===================== settings =====================
n          = 1;                         % file number in env.data.cleanFiles
winLen     = 15;                         % window length (s)
overlapPct = 50;                        % window overlap (%)
condSel    = [10 11];                   % conditions to plot
channels   = {'FC1','FC2'};
Hz         = 40;                        % frequency for the time course
xfreq      = [35 45];                   % FFT x-axis range (max 45, see cfg.foi below)
cols       = [0.85 0.33 0.10;           % orange
              0.00 0.45 0.74;           % blue
              0.47 0.67 0.19;           % green
              0.49 0.18 0.56];          % purple
chStr      = strjoin(channels, ', ');

%% load data
try
    dat = load(fullfile(env.data.cleanFiles(n).folder,...
        env.data.cleanFiles(n).name), 'dat_after_ICA2').dat_after_ICA2;
catch
    dat = load(fullfile(env.data.cleanFiles(n).folder,...
        env.data.cleanFiles(n).name), 'dat_after_ICA').dat_after_ICA;
end

%% windowing
cfg = [];
cfg.winlen  = winLen;
cfg.overlap = overlapPct;
[windowDat, windowSummary] = NaN_windowing(cfg, dat);

%% FFT per window
cfg = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foi         = 2:0.5:45;
cfg.pad         = 'nextpow2';
cfg.keeptrials  = 'yes';
freq_win = ft_freqanalysis(cfg, windowDat);   % window-level

%% average per condition
C = unique(freq_win.trialinfo(:,1));
for c = 1:numel(C)
    freq_cond{c} = ft_freqdescriptives( ...
        struct('avgoverrpt','yes'), ...
        ft_selectdata(struct('trials',find(freq_win.trialinfo(:,1)==C(c))), freq_win));
end

%% ===================== FFT + 40 Hz per window over time =====================
wlen     = round(winLen * dat.fsample);                         % window length in samples
stepSamp = max(round(wlen * (1 - overlapPct/100)), 1);          % same step as NaN_windowing
[~, fi]  = min(abs(freq_win.freq - Hz));                        % index of the 40 Hz bin

figure('Color','w','Position',[100 100 1400 550]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('S%02d, channels: %s', n, chStr), 'FontSize', 16, 'FontWeight','bold');
axL = nexttile; hold(axL,'on');
axR = nexttile; hold(axR,'on');
h = gobjects(1,numel(condSel));

for k = 1:numel(condSel)
    c = condSel(k);

    % --- FFT of every kept window of this condition
    cfg = [];
    cfg.trials      = find(freq_win.trialinfo(:,1) == c);
    cfg.channel     = channels;
    cfg.avgoverchan = 'yes';
    x = ft_selectdata(cfg, freq_win);

    f = x.freq(:)';
    P = reshape(x.powspctrm, size(x.powspctrm,1), []);          % windows x freq

    % --- left: average over windows
    h(k) = plotMeanSEM(axL, f, P, cols(k,:));

    % --- right: 40 Hz value of each window at its centre time
    starts = 1:stepSamp:(size(dat.trial{c},2) - wlen + 1);      % all windows, as in NaN_windowing
    tWin   = dat.time{c}(starts + floor(wlen/2));               % centre time of each window
    p40    = NaN(size(starts));                                 % NaN = window rejected
    wIdx   = freq_win.trialinfo(cfg.trials, 2);                 % window number of each kept window
    p40(wIdx) = P(:, fi);

plot(axR, tWin, p40, 'Color', cols(k,:), 'LineWidth', 1.2);
end

% --- cosmetics: FFT
xline(axL, Hz, '-k', 'LineWidth', 1, 'HandleVisibility','off');
xlim(axL, xfreq);
xlabel(axL, 'Frequency (Hz)'); ylabel(axL, 'Power (\muV^2)');
title(axL, 'FFT (average over windows)', 'FontSize', 14);
legend(axL, h, conds(condSel), 'Location','northwest', 'Box','off', 'Interpreter','none');

% --- cosmetics: time course
xlabel(axR, 'Time (s)'); ylabel(axR, sprintf('Power at %g Hz (\\muV^2)', Hz));
title(axR, sprintf('%g Hz power per %g s window', Hz, winLen), 'FontSize', 14);

set([axL axR], 'FontSize', 12, 'LineWidth', 1, 'TickDir','out', 'Box','off');
linkaxes([axL axR], 'y');

fname = sprintf('15secWindow__S%02d_FFT_%gHz_time_cond%s.jpg', n, Hz, strjoin(string(condSel), '_'));
saveas(gcf, fullfile("C:\Users\yoelgo\Desktop\", fname));

%% ===================== local functions (must stay at END of the script) =====================
function h = plotMeanSEM(ax, f, X, col)
    m = mean(X, 1, 'omitnan');
    s = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X(:,1))));
    fill(ax, [f fliplr(f)], [m+s fliplr(m-s)], col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h = plot(ax, f, m, 'Color', col, 'LineWidth', 2);
end