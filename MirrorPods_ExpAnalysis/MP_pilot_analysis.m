%% loading in LSL file and preprocessing the data
clc; clear all; close all;
cd('C:\Users\yoelgo\Documents\GitHub\MP_analyses\');
addpath('C:\Users\yoelgo\Documents\GitHub\MP_analyses\extra_functions\');
env = setupEnviroment11();

%%
n = 1;
dat = load(fullfile(env.data.cleanFiles(n).folder,...
    env.data.cleanFiles(n).name), 'dat_after_ICA').dat_after_ICA;
%%
cfg = [];
cfg.winlen  = 5;
cfg.overlap = 1;
[windowDat, windowSummary] = NaN_windowing(cfg, dat);
%%
cfg = [];
cfg.trials = 5:size(windowDat.trial,2)
windowDat2 =ft_selectdata(cfg, windowDat)

%%
cfg = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foi         = 1:0.25:40;
cfg.pad         = '5';
cfg.keeptrials  = 'yes';
freq_win = ft_freqanalysis(cfg, windowDat2);   % <-- keep this as window-level

%%
conds = unique(freq_win.trialinfo(:,1));
for c = 1:numel(conds)
    freq_cond{c} = ft_freqdescriptives( ...
        struct('avgoverrpt','yes'), ...
        ft_selectdata(struct('trials',find(freq_win.trialinfo(:,1)==conds(c))), freq_win));
end

%%
Hz = 20;

if Hz == 30
    idx = [1 2];
else
    idx = [3 4];
end
cfg = [];
cfg.parameter = 'powspctrm';
cfg.xlim      = [1 35];
%cfg.ylim      = [0 0.3];
cfg.channel = {'EEG', '-FT10', '-F7', '-T8'};
cfg.layout    = env.EEG.lay
cfg.showlabels = 'yes';
cfg.showlegend = 'yes';
figure; ft_multiplotER(cfg, freq_cond{idx(1)}, freq_cond{idx(2)});

%% plot fft
channels = {'Cz'};
Hz = 20;
xfreq = [Hz-10 Hz+10]
%xfreq = [1 Hz+10]
figure;

if Hz == 30
    idx = [1 2];
else
    idx = [3 4];
end
if size(channels) > 1
    cfg = [];
    cfg.channel = channels;
    cfg.avgoverchan = 'yes';
else
    cfg = [];
    cfg.channel = channels;
end
x1 = ft_selectdata(cfg,freq_cond{idx(1)});
x2 = ft_selectdata(cfg,freq_cond{idx(2)});

y1 = plot(x1.freq,x1.powspctrm);

hold on

% no-beep
y2 = plot(x2.freq,x2.powspctrm);


% vertical marker (no handle needed in legend)
xline(Hz);
xline(112/60);
xlim(xfreq);
if size(channels,2) > 1
    ch_str = strjoin(channels, ', ');
    title(sprintf('fft, condition: %d\nElectrodes: (%s)', Hz, ch_str), ...
        'FontSize', 18);
else
title(sprintf('fft, condition: %gHz\nElectrode: %s', Hz, strjoin(channels, ', ')), ...
    'FontSize', 18);
end
legend([y1 y2], {'beep','no-beep'})
%set(gca, 'XScale', 'log', 'YScale', 'log');
%set(gca, 'YScale', 'log');

xlabel('Frequency (Hz)')
ylabel('Power')
axis square;

saveas(gcf,['ID' int2str(n) '_' int2str(Hz) 'Hz.png'])

%% LAVI analysis
LAVIdf = [];
LAVI = [];
env.cfgLAVI.foi = 0.5:0.5:40;
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

%