%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Pilot Study, Frequency Tagging, UCL, 04/26 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading in LSL file and preprocessing the data
clc; clear all; close all;
dir = 'C:\Users\yoelgo\Documents\GitHub\MP_analyses\MirorPods_FT_Pilot\'; 
addpath(genpath(dir));
cd(dir)
env = setupEnviroment11();

conds = {'(1)4bps O2 -30db', '(2)4bps O2 -40db',...    % 1 2
         '(3)4bps O6 -30db', '(4)4bps O6 -40db',...    % 3 4
         '(5)4bps O19 -30db', '(6)4bps O19 -40db',...  % 5 6
         '(7)4bps O20 -30db(30Hz)', '(8)4bps O20 -40db(30Hz)',...   % 7 8
         '(9)Barbara 39Hz 41Hz'...                  % 9
         '(10)Beep O2 -20db', '(11)noBeep O2 -20db',...  % 10 11
         '(12)Beep O6 -20db', '(13)noBeep O6 -20db',...  % 12 13
         '(14)Beep O19 -25db', '(15)noBeep O19 -25db',...   % 14 15
         '(16)Beep O20 -25db(30Hz)', '(17)noBeep O20 -25db(30Hz)'}      % 16 17
%%
n = 2;
try
dat = load(fullfile(env.data.cleanFiles(n).folder,...
    env.data.cleanFiles(n).name), 'dat_after_ICA2').dat_after_ICA2;
catch
    dat = load(fullfile(env.data.cleanFiles(n).folder,...
    env.data.cleanFiles(n).name), 'dat_after_ICA').dat_after_ICA;
end

%%
cfg = [];
cfg.winlen  = 5;
cfg.overlap = 50;
[windowDat, windowSummary] = NaN_windowing(cfg, dat);
%%
cfg = [];
cfg.trials = 1:size(windowDat.trial,2)
windowDat2 =ft_selectdata(cfg, windowDat)

%%
cfg = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foi         = 2:0.5:45;
cfg.pad         = 'nextpow2';
cfg.keeptrials  = 'yes';
freq_win = ft_freqanalysis(cfg, windowDat);   % <-- keep this as window-level

%%
C = unique(freq_win.trialinfo(:,1));
for c = 1:numel(C)
    freq_cond{c} = ft_freqdescriptives( ...
        struct('avgoverrpt','yes'), ...
        ft_selectdata(struct('trials',find(freq_win.trialinfo(:,1)==C(c))), freq_win));
end
%%
cfg = [];
%cfg.xlim = [1.3 1.7];
cfg.xlim = [35 45];
%cfg.zlim = [0 0.25];
cfg.layout = env.EEG.lay;
cfg.showlabels = 'yes'
%cfg.parameter = ''; % the default 'avg' is not present in the data
figure; ft_multiplotER(cfg,freq_cond{13}); colorbar

%% plot fft
idx = [10 12 14];
channels = {'FC1', 'FC2'};
Hz = 40;
%xfreq = [Hz-12 Hz+5]
xfreq = [35 45]
figure;
hold on;

xline(Hz);
%xline(30)

if numel(channels) > 1
    cfg = [];
    cfg.channel = channels;
    cfg.avgoverchan = 'yes';
else
    cfg = [];
    cfg.channel = channels;
end

h = gobjects(1,length(idx));
for i = 1:length(idx)
    x1 = ft_selectdata(cfg, freq_cond{idx(i)});
    h(i) = plot(x1.freq, x1.powspctrm, 'LineWidth', 1.2);
end



%xline(30, 'LineStyle','-');
if any(idx == 9)
    xline(39);
    xline(41);
end
xlim(xfreq);

if numel(channels) > 1
    ch_str = strjoin(channels, ', ');
    title(sprintf('fft, condition: %gHz\nElectrodes: (%s)', Hz, ch_str), 'FontSize', 18);
else
    title(sprintf('fft, condition: %gHz\nElectrode: %s', Hz, channels{1}), 'FontSize', 18);
end

if isequal(idx, [1 2 3])
    ttl = 'YoelConditions1';
elseif isequal(idx, [4 5 6])
    ttl = 'YoelConditions2';
elseif idx(end) == 9
    ttl = 'BarbaraRep';
elseif isequal(idx, 10:2:15)      % = [10 12 14]
    ttl = 'DaniellaConditions_Beep';
elseif isequal(idx, [10 11])      % = [10 12 14]
    ttl = 'DaniellaConditions_Beep_vs_nobeep';
elseif isequal(idx, 11:2:15)      % = [11 13 15]
    ttl = 'DaniellaConditions_NoBeep';
elseif isequal(idx, [16 17])
    ttl = 'DaniellaConditions_30Hz';
elseif isequal(idx, [7 8])
    ttl = 'YoelConditions_30Hz';
else
    %error('No defined condition to save for idx = %s', mat2str(idx));
end
lgd = legend(h, conds(idx), 'Location', 'northeastoutside');

% find line objects inside legend and thicken them
lgd_lines = findobj(lgd, 'Type', 'Line');
set(lgd_lines, 'LineWidth', 2.5);

%saveas(gcf, fullfile(env.paths.pltOut, ['S02_' ttl '.jpg']));
xlabel('Frequency (Hz)')
ylabel('Power')
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAVI analysis
LAVIdf = [];
LAVI = [];
env.cfgLAVI.foi = 0.5:0.25:12;
ftLAVI = freq_cond{1};
ftLAVI.freq = env.cfgLAVI.foi;
% Calculate LAVI for each electrode
for t = 1:length(dat.trial)
    for e=1:length(dat.label)
        disp(['Trial: ' num2str(t) ', Electrode: ' num2str(e)]);
        LAVI(e,:) = Prepare_LAVI(env.cfgLAVI, dat.trial{t}(e,:));
    end
    ftLAVI.powspctrm = LAVI;
    LAVIdf{t} = ftLAVI; 
end


%% plot LAVI
idx = [10:2:17]
channels = {'Cz'};
Hz = 40;
%xfreq = [Hz-20 Hz+5]
xfreq = [1 10]
figure;
hold on;

if numel(channels) > 1
    cfg = [];
    cfg.channel = channels;
    cfg.avgoverchan = 'yes';
else
    cfg = [];
    cfg.channel = channels;
end

h = gobjects(1,length(idx));
for i = 1:length(idx)
    x1 = ft_selectdata(cfg, LAVIdf{idx(i)});
    h(i) = plot(x1.freq, x1.powspctrm);
end

xline(Hz);
if any(idx == 9)
    xline(39);
    xline(41);
end

xlim(xfreq);

if numel(channels) > 1
    ch_str = strjoin(channels, ', ');
    title(sprintf('fft, condition: %gHz\nElectrodes: (%s)', Hz, ch_str), 'FontSize', 18);
else
    title(sprintf('fft, condition: %gHz\nElectrode: %s', Hz, channels{1}), 'FontSize', 18);
end

legend(h, conds(idx), 'Location', 'northeastoutside');

xlabel('Frequency (Hz)')
ylabel('EEG rhythmicity (LAVI)')
axis square;
%saveas(gcf, fullfile(env.paths.pltOut, 'LAVI_DaniellaConditions.jpg'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITPC analysis
winSize = 1/4;   % seconds, size of window for each trial.
winSizeSamp = round(winSize*dat.fsample);

ITC_df = [];
ITC_subsamp_df = [];


for i = 1:numel(dat.trial)
    cfg = [];
    cfg.trials = i;
    curDat = ft_selectdata(cfg,dat);

    sampsStart = curDat.sampleinfo(1):winSizeSamp:curDat.sampleinfo(end);
    sampsEnd   = curDat.sampleinfo(1)+winSizeSamp+1:winSizeSamp:curDat.sampleinfo(end);
    if length(sampsEnd) > length(sampsStart)
        sampsEnd   = sampsEnd(1:length(sampsStart));
    elseif length(sampsStart) > length(sampsEnd)
        sampsStart = sampsStart(1:length(sampsEnd));
    end

    curTrig_Table = table(sampsStart', sampsEnd', zeros(size(sampsEnd,2),1), ...
        'VariableNames',{'beg_sample', 'end_sample', 'offset'});
    cfg = [];
    cfg.trl = curTrig_Table{:,:};
    trlDat = ft_redefinetrial(cfg,curDat);

    % Drop NaN windows
    keep = true(1,numel(trlDat.trial));
    for tr = 1:numel(trlDat.trial)
        if any(isnan(trlDat.trial{tr}(:)))
            keep(tr) = false;
        end
    end
    trlDat.trial = trlDat.trial(keep);
    trlDat.time  = trlDat.time(keep);
    trlDat.sampleinfo  = trlDat.sampleinfo(keep);

    % compute itc for the drumming block
    cfg = [];
    cfg.foi = 15:0.5:45;
    cfg.pad = 10;
    itc = ITPC_function(cfg, trlDat);
    itc.cond = conds{i};
    ITC_df{i} = itc;

end

%%
idx = [9];
channels = {'Cz'};
Hz = 40;
xfreq = [15 Hz+5]
figure;
hold on;

if numel(channels) > 1
    cfg = [];
    cfg.channel = channels;
    cfg.avgoverchan = 'yes';
else
    cfg = [];
    cfg.channel = channels;
end

h = gobjects(1,length(idx));
for i = 1:length(idx)
    x1 = ft_selectdata(cfg, ITC_df{idx(i)});
    h(i) = plot(x1.freq, x1.itpc, 'LineWidth', 1.2);
end

%xline(Hz);

if any(idx == 9)
    xline(39);
    xline(41);
else
    xline(40);
    xline(30, 'LineStyle','--');
end

xlim(xfreq);

if numel(channels) > 1
    ch_str = strjoin(channels, ', ');
    title(sprintf('ITPC window size: %gms\nElectrodes: (%s)', round(winSize*1000),ch_str), 'FontSize', 18);
else
    title(sprintf('ITPC, window size: %gms\nElectrode: %s', round(winSize*1000), channels{1}), 'FontSize', 18);
end

lgd = legend(h, conds(idx), 'Location', 'northeastoutside');

% find line objects inside legend and thicken them
lgd_lines = findobj(lgd, 'Type', 'Line');
set(lgd_lines, 'LineWidth', 2.5);

xlabel('Frequency (Hz)')
ylabel('ITPC')
axis square;
if idx(end) == 4
    ttl = 'YoelConditions1';
elseif idx(end) == 8
    ttl = 'YoelConditions2';
elseif idx(end) == 9
    ttl = 'BarbaraRep';
elseif idx(1) == 10
    ttl = 'DaiellaConditions_Beep';
elseif idx(1) == 11
    ttl = 'DaniellaConditions_NoBeep'
else
    %error 'no defined condition to save'
    ttl = '30HzConds'
end

saveas(gcf, fullfile(env.paths.pltOut, ['S01NEW_ITPC41' ttl '.jpg']));

%%
cfg = [];
cfg.xlim = [1.3 1.7];
%cfg.xlim = [3.8 4.2];
cfg.zlim = [0 0.25];
cfg.layout = env.EEG.lay;
cfg.parameter = 'itpc'; % the default 'avg' is not present in the data
figure; ft_topoplotER(cfg,itpcGA{1}); colorbar

%%
figure;
ft_plot_layout(env.EEG.lay);

%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

conds = {'(1)4bps O2 -30db', '(2)4bps O2 -40db',...    % 1 2
         '(3)4bps O6 -30db', '(4)4bps O6 -40db',...    % 3 4
         '(5)4bps O19 -30db', '(6)4bps O19 -40db',...  % 5 6
         '(7)4bps O20 -30db(30Hz)', '(8)4bps O20 -40db(30Hz)',...   % 7 8
         '(9)Barbara 39Hz 41Hz'...                  % 9
         '(10)Beep O2 -20db', '(11)noBeep O2 -20db',...  % 10 11
         '(12)Beep O6 -20db', '(13)noBeep O6 -20db',...  % 12 13
         '(14)Beep O19 -25db', '(15)noBeep O19 -25db',...   % 14 15
         '(16)Beep O20 -25db(30Hz)', '(17)noBeep O20 -25db(30Hz)'}      % 16 17




conds = {'4bps drumming over 40Hz', '(2)4bps O2 -40db',...    % 1 2
         '(3)4bps O6 -30db', '(4)4bps O6 -40db',...    % 3 4
         '(5)4bps O19 -30db', '(6)4bps O19 -40db',...  % 5 6
         '(7)4bps O20 -30db(30Hz)', '(8)4bps O20 -40db(30Hz)',...   % 7 8
         '(9)Barbara replication, 39Hz 41Hz tagging'...                  % 9
         '2Hz Beep over 40Hz', 'no beep over 40Hz',...  % 10 11
         '(12)Beep O6 -20db', '(13)noBeep O6 -20db',...  % 12 13
         '(14)Beep O19 -25db', '(15)noBeep O19 -25db',...   % 14 15
         '(16)Beep O20 -25db(30Hz)', '(17)noBeep O20 -25db(30Hz)'}      % 16 17
%%
clear title xlabel ylabel legend   % in case any of these got shadowed by variables

%% ===================== settings =====================
condPair  = [10 11];                       % 10 = Beep O2, 11 = noBeep O2
condNames = {'Beep (40 Hz)', 'No beep (40 Hz)'};
channels  = {'FC1','FC2'};
Hz        = 40;
xfreq     = [25 45];
cols      = [0.85 0.33 0.10;               % Beep  - orange
             0.00 0.45 0.74];              % noBeep - blue
chStr     = strjoin(channels, ', ');

% ===================== FFT figure =====================
figure('Color','w','Position',[100 100 650 550]);
ax = axes; hold(ax,'on');
h = gobjects(1,2);

for k = 1:2
    cfg = [];
    cfg.trials      = find(freq_win.trialinfo(:,1) == condPair(k));
    cfg.channel     = channels;
    cfg.avgoverchan = 'yes';
    x = ft_selectdata(cfg, freq_win);

    f = x.freq(:)';
    P = reshape(x.powspctrm, size(x.powspctrm,1), []);   % windows x freq

    h(k)    = plotMeanSEM(ax, f, P, cols(k,:));
    nWin(k) = size(P,1);
end

xline(ax, Hz, '-k', 'LineWidth', 1, 'HandleVisibility','off');
xlim(ax, xfreq);
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Power (\muV^2)');
title(ax, sprintf('FFT, channels: %s', chStr), 'FontSize', 16);
legend(ax, h, compose('%s (n=%d)', string(condNames(:)), nWin(:)), ...
       'Location','northwest', 'Box','off');
set(ax, 'FontSize', 12, 'LineWidth', 1, 'TickDir','out', 'Box','off');
axis(ax, 'square');

% exportgraphics(gcf, fullfile(env.paths.pltOut, 'S02_Beep_vs_noBeep_FFT.png'), 'Resolution', 300);
%% ===================== TFR figure (continuous data) =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear title xlabel ylabel legend   % in case any of these got shadowed by variables

%% ===================== settings =====================
condPair  = [1 3 5];                       % 10 = Beep O2, 11 = noBeep O2
condNames = {'Beep (40 Hz)', 'No beep (40 Hz)'};
channels  = {'FC1','FC2'};
Hz        = 40;
xfreq     = [30 45];                       % FFT x-axis range (max 45, see cfg.foi of freq_win)
cols      = [0.85 0.33 0.10;               % Beep   - orange
             0.00 0.45 0.74];              % noBeep - blue
chStr     = strjoin(channels, ', ');

% ===================== FFT figure =====================
figure('Color','w','Position',[100 100 650 550]);
ax = axes; hold(ax,'on');
h = gobjects(1,2);
nWin = zeros(1,2);

for k = 1:2
    cfg = [];
    cfg.trials      = find(freq_win.trialinfo(:,1) == condPair(k));
    cfg.channel     = channels;
    cfg.avgoverchan = 'yes';
    x = ft_selectdata(cfg, freq_win);

    f = x.freq(:)';
    P = reshape(x.powspctrm, size(x.powspctrm,1), []);   % windows x freq

    h(k)    = plotMeanSEM(ax, f, P, cols(k,:));
    nWin(k) = size(P,1);
end

xline(ax, Hz, '-k', 'LineWidth', 1, 'HandleVisibility','off');
xlim(ax, xfreq);
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Power (\muV^2)');
title(ax, sprintf('FFT,  channels: %s', chStr), 'FontSize', 16);
legend(ax, h, compose('%s', string(condNames(:))), ...
       'Location','northwest', 'Box','off');
set(ax, 'FontSize', 12, 'LineWidth', 1, 'TickDir','out', 'Box','off');
axis(ax, 'square');


%% ===================== TFR figure (continuous data) =====================
winLen = 3;                          % analysis window length in s (4 s -> ~0.25 Hz resolution)
hw     = winLen/2;

cfgF = [];
cfgF.method     = 'mtmconvol';
cfgF.output     = 'pow';
cfgF.taper      = 'hanning';
cfgF.foi        = 35:0.25:45;
cfgF.t_ftimwin  = ones(size(cfgF.foi)) * winLen;
cfgF.keeptrials = 'no';                     % finer steps to match the resolution

cfgF.keeptrials = 'no';

tfr = cell(1,2);
for k = 1:2
    cfg = [];
    cfg.trials  = condPair(k);       % trial n in dat = condition n
    cfg.channel = channels;
    d = ft_selectdata(cfg, dat);

    % replace NaNs with 0 so the FFT-based convolution works
    t       = d.time{1};
    nanMask = any(isnan(d.trial{1}), 1);
    d.trial{1}(:, nanMask) = 0;

    cfgF.toi = (t(1)+hw) : 0.1 : (t(end)-hw);
    tfr{k}   = ft_freqanalysis(cfgF, d);

    % put NaN back wherever the analysis window overlapped NaN samples
    cNan = [0 cumsum(nanMask)];
    for ti = 1:numel(tfr{k}.time)
        i1 = find(t >= tfr{k}.time(ti) - hw, 1, 'first');
        i2 = find(t <= tfr{k}.time(ti) + hw, 1, 'last');
        if cNan(i2+1) - cNan(i1) > 0
            tfr{k}.powspctrm(:,:,ti) = NaN;
        end
    end
end

Pdb  = cellfun(@(x) 10*log10(squeeze(mean(x.powspctrm, 1))), tfr, 'UniformOutput', false);
clim = [min([Pdb{:}],[],'all') max([Pdb{:}],[],'all')];

figure('Color','w','Position',[100 100 1400 420]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Time-frequency, channels:  %s', chStr), 'FontSize', 16, 'FontWeight','bold');

for k = 1:2
    ax = nexttile;
    imagesc(tfr{k}.time, tfr{k}.freq, Pdb{k}, 'AlphaData', ~isnan(Pdb{k}));
    axis xy; caxis(clim);
    set(ax, 'Color', [0.8 0.8 0.8]);           % NaN segments show as gray
    cb = colorbar; cb.Label.String = 'Power (dB)';
    yline(Hz, '--w', 'LineWidth', 1);
    title(condNames{k}, 'FontSize', 14);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    set(ax, 'FontSize', 12, 'TickDir','out');
end

%% ===================== local functions (must stay at END of the script) =====================
function h = plotMeanSEM(ax, f, X, col)
    m = mean(X, 1, 'omitnan');
    s = std(X, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(X(:,1))));
    fill(ax, [f fliplr(f)], [m+s fliplr(m-s)], col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h = plot(ax, f, m, 'Color', col, 'LineWidth', 2);
end

%%
%% ===================== FFT figures =====================
conds = {'4bps drumming over 40Hz', '(2)4bps O2 -40db',...    % 1 2
         '(3)4bps O6 -30db', '(4)4bps O6 -40db',...    % 3 4
         '(5)4bps O19 -30db', '(6)4bps O19 -40db',...  % 5 6
         '(7)4bps O20 -30db(30Hz)', '(8)4bps O20 -40db(30Hz)',...   % 7 8
         'Barbara replication, 39Hz 41Hz tagging'...                  % 9
         '2Hz beep over 40Hz', 'only 40Hz',...  % 10 11
         '(12)Beep O6 -20db', '(13)noBeep O6 -20db',...  % 12 13
         '(14)Beep O19 -25db', '(15)noBeep O19 -25db',...   % 14 15
         '(16)Beep O20 -25db(30Hz)', '(17)noBeep O20 -25db(30Hz)'}      % 16 17
channels = {'FC1', 'FC2'};
chStr     = strjoin(channels, ', ');
xfreq = [35 45]
figSets  = { 1, [10 11], 9 };           % conditions per figure
figMarks = { 40, 40, [39 41] };         % dashed line(s) per figure (Hz)
cols     = [0.85 0.33 0.10;             % orange
            0.00 0.45 0.74];            % blue

for fIdx = 1:numel(figSets)
    cIdx = figSets{fIdx};

    figure('Color','w','Position',[100 100 650 550]);
    ax = axes; hold(ax,'on');
    h = gobjects(1,numel(cIdx));

    for k = 1:numel(cIdx)
        cfg = [];
        cfg.trials      = find(freq_win.trialinfo(:,1) == cIdx(k));
        cfg.channel     = channels;
        cfg.avgoverchan = 'yes';
        x = ft_selectdata(cfg, freq_win);

        f = x.freq(:)';
        P = reshape(x.powspctrm, size(x.powspctrm,1), []);   % windows x freq

        h(k) = plotMeanSEM(ax, f, P, cols(k,:));
    end

    for m = figMarks{fIdx}
        xline(ax, m, '-k', 'LineWidth', 1, 'HandleVisibility','off');
    end

    xlim(ax, xfreq);
    xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Power (\muV^2)');
    title(ax, sprintf('S01 FFT, channels:  %s', chStr), 'FontSize', 16);
    legend(ax, h, conds(cIdx), 'Location','northwest', 'Box','off', 'Interpreter','none');
    set(ax, 'FontSize', 12, 'LineWidth', 1, 'TickDir','out', 'Box','off');
    axis(ax, 'square');
    fname = sprintf('S01_FFT_cond%s.jpg', strjoin(string(cIdx), '_'));
    saveas(gcf, fullfile('C:\Users\yoelgo\Documents\GitHub\MP_analyses\MirorPods_FT_Pilot\', fname));
end

%%
channels = {'FC1', 'FC2'};
chStr     = strjoin(channels, ', ');
figSets  = { 1, [10 11], 9 };           % conditions per figure
figMarks = { 40, 40, [39 41] };         % dashed line(s) per figure (Hz)

winLen = 10;                            % analysis window length in s
cfgF = [];
cfgF.method     = 'mtmconvol';
cfgF.output     = 'pow';
cfgF.taper      = 'dpss';
cfgF.foi        = 35:0.2:45;
cfgF.tapsmofrq  = ones(size(cfgF.foi)) * 0.3;
cfgF.t_ftimwin  = ones(size(cfgF.foi)) * winLen;
cfgF.keeptrials = 'no';

for fIdx = 1:numel(figSets)
    cIdx = figSets{fIdx};

    % compute TFR for each condition in this figure
    tfr = cell(1,numel(cIdx));
    Pdb = cell(1,numel(cIdx));
    for k = 1:numel(cIdx)
        tfr{k} = computeTFR(dat, cIdx(k), channels, cfgF, winLen);
        Pdb{k} = squeeze(mean(tfr{k}.powspctrm, 1));
    end
    allP = cellfun(@(p) p(:), Pdb, 'UniformOutput', false);
    clim = prctile(vertcat(allP{:}), [2 98]);                   % shared color scale

    figure('Color','w','Position',[100 100 700*numel(cIdx) 450]);
    tl = tiledlayout(1,numel(cIdx),'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('S01 TFR,  %s', chStr), 'FontSize', 16, 'FontWeight','bold');

    for k = 1:numel(cIdx)
        ax = nexttile;
        imagesc(tfr{k}.time, tfr{k}.freq, Pdb{k}, 'AlphaData', ~isnan(Pdb{k}));
        axis xy; caxis(clim);
        set(ax, 'Color', [0.8 0.8 0.8]);       % NaN segments show as gray
        cb = colorbar;
        cb.Label.String = 'Power (\muV^2)';
        for m = figMarks{fIdx}
            yline(m, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
        end
        title(conds{cIdx(k)}, 'FontSize', 14, 'Interpreter','none');
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        set(ax, 'FontSize', 12, 'TickDir','out');
    end

    fname = sprintf('S02_TFR_cond%s.jpg', strjoin(string(cIdx), '_'));
    saveas(gcf, fullfile(env.paths.pltOut, fname));
end

function tfr = computeTFR(dat, cond, channels, cfgF, winLen)
    hw = winLen/2;

    cfg = [];
    cfg.trials  = cond;                 % trial n in dat = condition n
    cfg.channel = channels;
    d = ft_selectdata(cfg, dat);

    % replace NaNs with 0 so the FFT-based convolution works
    t       = d.time{1};
    nanMask = any(isnan(d.trial{1}), 1);
    d.trial{1}(:, nanMask) = 0;

    cfgF.toi = (t(1)+hw) : 0.1 : (t(end)-hw);
    tfr      = ft_freqanalysis(cfgF, d);

    % put NaN back wherever the analysis window overlapped NaN samples
    cNan = [0 cumsum(nanMask)];
    for ti = 1:numel(tfr.time)
        i1 = find(t >= tfr.time(ti) - hw, 1, 'first');
        i2 = find(t <= tfr.time(ti) + hw, 1, 'last');
        if cNan(i2+1) - cNan(i1) > 0
            tfr.powspctrm(:,:,ti) = NaN;
        end
    end
end
%
figs = findobj('Type', 'figure');
for i = 1:numel(figs)
    saveas(figs(i), fullfile('C:\Users\yoelgo\Documents\GitHub\MP_analyses\MirorPods_FT_Pilot\', sprintf('S02_TFR_figure%d.jpg', figs(i).Number)));
end