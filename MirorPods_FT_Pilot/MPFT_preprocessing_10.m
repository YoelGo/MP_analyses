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

%% Load single participant
n = 1;
ID = extractBefore(env.data.rawFiles(n).name, '_');
env.ID = ID;
env.paths.curData = [env.data.rawFiles(n).folder '\'];
%
dat1 = load_xdf([env.paths.curData  env.data.rawFiles(n).name]);

% divide into streams
markers    = dat1{cellfun(@(x) strcmp(x.info.name, 'My_Markers'), dat1)};
EEG    = dat1{cellfun(@(x) strcmp(x.info.name, 'actiCHamp-21020490'), dat1)};
audio  = dat1{cellfun(@(x) strcmp(x.info.name, 'MyAudioStream'), dat1)};

%%
ts  = markers.time_series(:);
tst = markers.time_stamps(:);

keep = ~contains(ts, 'exp');
ts  = ts(keep);
tst = tst(keep);

if iscell(tst)
    tst = cellfun(@double, tst);
else
    tst = double(tst);
end

isStart = endsWith(ts, '_start');

block    = strings(0,1);
beg_samp = [];
end_samp = [];
offset   = [];
Omer     = [];
DB       = [];
metro    = strings(0,1);

for i = find(isStart)'
    s = ts{i};
    base = erase(s, '_start');
    endMarker = base + "_end";

    j = find(strcmp(ts, endMarker), 1);
    if isempty(j)
        continue
    end

    tBeg = tst(i);
    tEnd = tst(j);

    b = nearest(EEG.time_stamps, tBeg);
    e = nearest(EEG.time_stamps, tEnd);
    if iscell(b), b = b{1}; end
    if iscell(e), e = e{1}; end

    block(end+1,1)    = string(base);
    beg_samp(end+1,1) = double(b);
    end_samp(end+1,1) = double(e);
    offset(end+1,1)   = 0;

    tok = regexp(s,'Omer(\d+)','tokens','once');
    if isempty(tok)
        Omer(end+1,1) = -1;
    else
        Omer(end+1,1) = str2double(tok{1});
    end

    tok = regexp(s,'(-\d+)(?=d[bB]|DB|dbNoise)','tokens','once');
    if isempty(tok)
        DB(end+1,1) = -1;
    else
        DB(end+1,1) = str2double(tok{1});
    end

    tok = regexp(s,'(D[1-4]|4bps|noMetro)','match','once');
    if isempty(tok)
        metro(end+1,1) = "-1";
    else
        metro(end+1,1) = string(tok);
    end
end

trl = [beg_samp(:) end_samp(:) offset(:)];

trl_table = table(block, beg_samp, end_samp, offset, Omer, DB, metro, ...
    'VariableNames', {'block','beg_sample','end_sample','offset','Omer','DB','metro'});
trl_table = movevars(trl_table, 1, 'After', width(trl_table));

clear e DB b base beg_samp j keep metro Omer offset ts tst tok tEnd tBeg s isStart i end_samp endMarker block

% convert to ft data and redefine trials
ftEEG = LSL2ft(EEG);

cfg = [];
cfg.trl = trl;
datTrl = ft_redefinetrial(cfg, ftEEG);

datTrl.hdr.chantype(1:64) = {'EEG'};
datTrl.hdr.chantype(65:68) = {'AUX'};

trl
%%
elecs = {'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1',...
         'Pz','P3','P7','O1','Oz','O2','P4','P8','TP10','CP6','CP2','Cz',...
         'C4','T8','FT10','FC6','FC2','F4','F8','Fp2''AF1','AF3','AFz',...
         'F1','F5','FT7','FC3','C1','C5','TP7','CP3','P1',...
         'P5','PO7','PO3','POz','PO4','PO8','P6','P2','CPz','CP4','TP8','C6',...
         'C2','FC4','FT8','F6','AF8','AF4','F2','Iz'};
elecs = [elecs, elecs2]
%datTrl.label(1:32) = elecs

datTrl.label(1:64) = elecs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% basic preproc
clc;
cfg            = [];
cfg.channel    = 'EEG';
cfg.detrend    = 'yes';
cfg.demean     = 'yes';
cfg.continuous = 'no';
cfg.reref      = 'yes';
cfg.refchannel = 'all';
% Notch filters to remove 50Hz and harmonics
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase

pEEG = ft_preprocessing(cfg, datTrl); %p(preprocessed)EEG

% Process eye data and append.
cfg = [];
cfg.channel = {'AUX_1' 'AUX_2' 'AUX_3' 'AUX_4'}
EoG = ft_selectdata(cfg, datTrl)
pEoG = processEoG(EoG)

pEEG = ft_appenddata([], pEEG, pEoG)

clear EoG pEoG

%% Manual: View Data
cfg = [];
cfg.ylim  = [-50 50];
cfg.blocksize = 30;
man_artifact = ft_databrowser(cfg,pEEG)
%save([env.paths.artifacts  ID '_man_artifact'], "man_artifact");
load([env.paths.artifacts  ID '_man_artifact'], "man_artifact");
%% Manual: Remove Artifacts
cfg = []; 
cfg.artfctdef.reject           = 'nan';
cfg.artfctdef.visual.artifact = man_artifact.artfctdef.visual.artifact;
mEEG = ft_rejectartifact(cfg,pEEG);

% remove EoG
cfg = [];
cfg.channel = {'all' '-eogV', '-eogH'};
mEEG = ft_selectdata(cfg, mEEG)


%% Semi automatic artifact rejection: do you use this method???
cfg        = [];
cfg.metric = 'zvalue';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
mEEG       = ft_rejectvisual(cfg,mEEG);
%% Run ICA
cfg = [];
cfg.method  = 'runica';
cfg.channel = 'EEG';
cfg.numcomponent = 20;
%cfg.trials = 9:17
%cfg.runica.maxsteps = 100;
comp = ft_componentanalysis(cfg, mEEG);
%% view ICA components
% view time seriers and topopraphy of ICs
cfg = [];
cfg.viewmode = 'component';
cfg.allowoverlap = 'yes';
%cfg.continuous = 'yes';
cfg.blocksize = 30;
%cfg.channel = 1:10;
cfg.layout = env.EEG.lay;
ft_databrowser(cfg,comp);

%% reject components
cfg = [];
cfg.component = [3];
dat_after_ICA = ft_rejectcomponent(cfg, comp);

%%
cfg        = [];
cfg.metric = 'zvalue';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
dat_after_ICA = ft_rejectvisual(cfg,dat_after_ICA);

%%
channel_number = [21 61];
for i=1:length(channel_number)
disp([dat_after_ICA.label{channel_number(i)}]);
end
%%
cfg = [];
cfg.badch = {'TP9'};
cfg.section = {'all'};
cfg.elec = env.EEG.elec;
cfg.layout = env.EEG.lay;
cfg.blocksize = 30; 
dat_after_ICA2 = fixChannels2(cfg, dat_after_ICA);
%% view the data again
cfg = [];
cfg.ylim  = [-30 30];
cfg.blocksize = 30;
man_artifact = ft_databrowser(cfg,dat_after_ICA2)


%% save data after preproc
save([env.paths.cleanData ID '_EEG2_clean.mat'], "dat_after_ICA2");
disp(['Saved All Data!']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
A = EEG;
A.time_series = audio.time_series;
A.time_stamps = audio.time_stamps;
A.info.desc.channels.channel = A.info.desc.channels.channel(1);
A.info.effective_srate = 8000
%
ftAudio = LSL2ft(A);
%%
for k = 1:numel(ftAudio.trial)
    ftAudio.trial{k} = double(ftAudio.trial{k});
end

cfg = [];
cfg.resamplefs = 500;   % target sampling rate
ftAudio= ft_resampledata(cfg, ftAudio);

%%
cfg =[];

pAudio = ft_preprocessing(cfg, ftAudio)
ft_databrowser([],pAudio)
%%
[y, Fs] = audioread("Z:\Experiments\Yoel\Pilots\MirrorPods_FrequencyTagging_Pilot\Audio_Files\Experiment_Stimulus\Omer02_-20DB_noMetro_D1.wav");

