%% loading in LSL file and preprocessing the data
clc; clear all; close all;
% set relevant paths %
gitPath = 'C:\Users\yoelgo\Documents\GitHub\'; % path of the local GitHub

cd([gitPath 'MP_analyses\MirrorPods_ExpAnalysis\']); % path of this specific project in Github
addpath([gitPath 'MP_analyses\extra_functions\']); % add the function in the "extra_functions" folder in the github project
env = setupEnviroment11(); % a custom function which sets all the paths and prepares the variable enviroment (env) for this experiment
clear gitPath
%% Load single participant
% ***see the participant number at env.data.rawFiles***
n = 55;

% find the participant ID
ID = extractBefore(env.data.rawFiles(n).name, '_data');
env.ID = ID;
env.paths.curData = [env.data.rawFiles(n).folder '\'];

% load data
dat1 = load_xdf([env.paths.curData  env.data.rawFiles(n).name]);

%% divide into the different LSL streams
% EEG of player 1 and 2
EEG1    = dat1{cellfun(@(x) strcmp(x.info.name, 'Player1 EEG'), dat1)};
EEG2    = dat1{cellfun(@(x) strcmp(x.info.name, 'Player2 EEG'), dat1)};

% convert LSL stream to fieldtrip data
ftEEG1 = LSL2ft(EEG1);
ftEEG2 = LSL2ft(EEG2);

% change the name of the EOG and ECG channels
idx = strcmp(ftEEG1.label, 'TP10');  % old name
ftEEG1.label{idx} = 'ECG';      % new name
ftEEG2.label{idx} = 'ECG';      % new name
idx = strcmp(ftEEG1.label, 'TP9');  % old name
ftEEG1.label{idx} = 'EoG';      % new name
ftEEG2.label{idx} = 'EoG';      % new name
clear idx

% take ECG data for inspection
cfg = [];
cfg.channel = {'ECG'};
ECG1 = ft_selectdata(cfg,ftEEG1);
ECG2 = ft_selectdata(cfg,ftEEG2);

% markers
markers = dat1(cellfun(@(x) strcmp(x.info.name, 'MaxMarkers'), dat1));
markers = markers{1};

%% first look at the data
cfg = [];
cfg.ylim  = [-130 130];
cfg.channel = {'all' '-ECG'}; % look at data without ECG
cfg.preproc.demean = 'yes';
cfg.blocksize = 50;
ft_databrowser(cfg, ftEEG1)
%% divide into blocks
idxBeg = find(contains(markers.time_series, 'start'));
idxEnd = find(contains(markers.time_series, 'stop'));

sampBeg  = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps(idxBeg));
sampEnd = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps(idxEnd));

trlTable1 = table(sampBeg', sampEnd',...
    zeros(numel(sampBeg),1),...
    'VariableNames',{'beg_sample', 'end_sample', 'offset'});

sampBeg  = arrayfun(@(t) nearest(EEG2.time_stamps, t), markers.time_stamps(idxBeg));
sampEnd = arrayfun(@(t) nearest(EEG2.time_stamps, t), markers.time_stamps(idxEnd));

trlTable2 = table(sampBeg', sampEnd',...
    zeros(numel(sampBeg),1),...
    'VariableNames',{'beg_sample', 'end_sample', 'offset'});

% find the different conditions from the markers
conditions = extractBefore(markers.time_series(idxBeg),'start');
% find the participant (P1/P2/dyad)
P = repmat("0",size(conditions));
P(contains(conditions,"P1"))="P1";
P(contains(conditions,"P2"))="P2";
P(contains(conditions,"Dyad"))="Dyad";
% find the condition type
type = strtrim(erase(erase(erase(string(conditions),"P1"),"P2"),"Dyad"));
type = strtrim(erase(erase(erase(string(conditions),"P1"),"P2"),"Dyad"));
type(type=="") = "0";

trlTable1.P = P';
trlTable1.type = type';

trlTable2.P = P';
trlTable2.type = type';

%% redefine trials
% redifine trials for EEG1
cfg = [];
cfg.trl = trlTable1{:,1:3};
datTrl1 = ft_redefinetrial(cfg, ftEEG1);

% redifine trials for EEG2
cfg = [];
cfg.trl = trlTable2{:,1:3};
datTrl2 = ft_redefinetrial(cfg, ftEEG2);

clear EEG1 EEG2

datTrl1.trialinfo = trlTable1(:,4:5);   % keep P and type as table
datTrl2.trialinfo = trlTable2(:,4:5);   % keep P and type as table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% basic preproc
% first keep the EOG to add later
cfg = [];
cfg.channel = 'EoG';
cfg.demean = 'yes';
EOG1 = ft_preprocessing(cfg,datTrl1);
EOG2 = ft_preprocessing(cfg,datTrl2);

clc;
cfg            = [];
cfg.channel    = {'EEG', '-ECG','-EoG', '-GyroX', '-GyroY', '-GyroZ'} ; % remove EOG and ECG for refferencing
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.refchannel = {'all'}; % A1 and A2 are the ear-clips
% Notch filters to remove 50Hz and harmonics
cfg.lpfreq      = 70;
cfg.lpfilter    = 'yes';
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.5;
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase

pEEG1 = ft_preprocessing(cfg, datTrl1); %p(preprocessed)EEG
pEEG2 = ft_preprocessing(cfg, datTrl2); %p(preprocessed)EEG

% add EoG back
pEEG1 = ft_appenddata([],pEEG1,EOG1);
pEEG2 = ft_appenddata([],pEEG2,EOG2);
%% 
% participant 1 or 2
I = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if I == 1
    curDat = pEEG1;
else
    curDat = pEEG2;
end

%% Manual: View Data
cfg = [];
cfg.ylim  = [-50 50];
cfg.continuous = 'yes';
cfg.blocksize = 30;

man_artifact = ft_databrowser(cfg,curDat)

save([env.paths.artifacts  ID '_man_artifact' num2str(I)], "man_artifact");
%load([env.paths.artifacts  ID '_man_artifact'], "man_artifact");
%% Manual: Remove Artifacts
cfg = []; 
cfg.artfctdef.reject           = 'nan';
cfg.artfctdef.visual.artifact = man_artifact.artfctdef.visual.artifact;
mEEG = ft_rejectartifact(cfg,curDat);

cfg = [];
cfg.channel = {'all', '-EoG','-GyroX', '-GyroY', '-GyroZ'};
mEEG = ft_preprocessing(cfg,mEEG);

%% check for further noisy channels
cfg        = [];
cfg.metric = 'var';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
mEEG       = ft_rejectvisual(cfg,mEEG);

mEEG.label(1)         % dicsover the channel from the number
%% inspct the bad channel
cfg = [];
cfg.ylim  = [-50 50];
cfg.continuous = 'yes';
cfg.blocksize = 30;
ft_databrowser(cfg,mEEG)

%% remove noisy channels to be interpolated later.
badCh = {'T7'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('badCh','var') && ~isempty(badCh{1})
    cfg = [];
    cfg.channel = badCh;
    badChDat = ft_selectdata(cfg,mEEG);

    cfg = [];
    cfg.channel = {'all', '-EoG', '-T7', '-T8'};  %%%% YOEL FIX
    mEEG = ft_selectdata(cfg, mEEG);
else
    warning('No bad channels to interpolate?');
end
%% Run ICA
cfg = [];
cfg.method  = 'runica';
cfg.channel = {'all' '-EoG'}; % EEG no EoG, no badch
cfg.numcomponent = 20;
%cfg.runica.maxsteps = 100;     % only uncomment if the ICA is too long and doesn't converge
comp = ft_componentanalysis(cfg, mEEG);
%% view ICA components
% view time seriers and topopraphy of ICs
cfg = [];
cfg.viewmode = 'component';
cfg.allowoverlap = 'yes';
cfg.continuous = 'yes';
cfg.blocksize = 30;
cfg.layout = env.EEG.lay;
ft_databrowser(cfg,comp);

%% reject components
cfg = [];
cfg.component = [1 2]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat_after_ICA = ft_rejectcomponent(cfg, comp);

%% reintroduce the bad channels and interpolate them
% reintroduce the bad channels
dat_after_ICA = ft_appenddata([], dat_after_ICA, badChDat);

% interpolate the bad channels using a custom function. see notation iside.
cfg = [];
cfg.badch = badCh;
cfg.section = {'all', 'all'};
cfg.elec = env.EEG.elec;
cfg.layout = env.EEG.lay;
cfg.blocksize = 30; 
dat_after_ICA = fixChannels(cfg, dat_after_ICA);
%%
cfg        = [];
cfg.metric = 'zvalue';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
dat_after_ICA = ft_rejectvisual(cfg,dat_after_ICA);

%% view the data again
cfg = [];
cfg.ylim  = [-30 30];
cfg.continuous = 'yes';
cfg.blocksize = 30;
man_artifact = ft_databrowser(cfg,dat_after_ICA)

%%
cfg = []; 
cfg.artfctdef.reject           = 'nan';
cfg.artfctdef.visual.artifact = man_artifact.artfctdef.visual.artifact;
dat_after_ICA = ft_rejectartifact(cfg,dat_after_ICA);

%% save data after preproc
save([env.paths.cleanData ID '_EEG' num2str(I) '_clean.mat'], "dat_after_ICA");
save([env.paths.cleanData ID '_trlTable1'], "trlTable1");
save([env.paths.cleanData ID '_trlTable2'], "trlTable2");

disp(['Saved All Data!']);