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
n = 1;
try
dat = load(fullfile(env.data.cleanFiles(n).folder,...
    env.data.cleanFiles(n).name), 'dat_after_ICA2').dat_after_ICA2;
catch
    dat = load(fullfile(env.data.cleanFiles(n).folder,...
    env.data.cleanFiles(n).name), 'dat_after_ICA').dat_after_ICA;
end


%%
cfg = [];
cfg.winlen  = 10;
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
idx = [1 2 10 11]
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


if all(idx == [1,2,3])
    ttl = 'YoelConditions1';
elseif all(idx == [4,5,6])
    ttl = 'YoelConditions2';
elseif idx(end) == 9
    ttl = 'BarbaraRep';
elseif all(idx == 10:2:15)
    ttl = 'DaiellaConditions_Beep';
elseif all(idx == 11:2:15)
    ttl = 'DaniellaConditions_NoBeep'

elseif all(idx == [16 17])
    ttl = 'DaniellaConditions_30Hz';
elseif all(idx == [7 8])
    ttl = 'YoelConditions_30Hz';
else
    error no defined condition to save
end

lgd = legend(h, conds(idx), 'Location', 'northeastoutside');

% find line objects inside legend and thicken them
lgd_lines = findobj(lgd, 'Type', 'Line');
set(lgd_lines, 'LineWidth', 2.5);

saveas(gcf, fullfile(env.paths.pltOut, ['S02_' ttl '.jpg']));
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