%% loading in LSL file and preprocessing the data
clc; clear all; close all;
cd('C:\Users\yoelgo\Documents\GitHub\MP_analyses\');
addpath('C:\Users\yoelgo\Documents\GitHub\MP_analyses\extra_functions\');
env = setupEnviroment11();

%% Load single participant
n = 7;
ID = extractBefore(env.data.rawFiles(n).name, '_data');
env.ID = ID;
env.paths.curData = [env.data.rawFiles(n).folder '\'];
dat1 = load_xdf([env.paths.curData  env.data.rawFiles(n).name]);

clear csv;
%% divide into streams
% EEG of player 1 and 2
EEG1    = dat1{cellfun(@(x) strcmp(x.info.name, 'Player1 EEG'), dat1)};
EEG2    = dat1{cellfun(@(x) strcmp(x.info.name, 'Player2 EEG'), dat1)};

% markers
markers = dat1(cellfun(@(x) strcmp(x.info.name, 'MaxMarkers'), dat1));
markers = markers{1};

% Find txt files containing "Noise"
txt = dir(fullfile(env.data.rawFiles(n).folder, '**', '*Noise*.txt'));

txtFiles = cell(numel(txt),1);
%
clc;
fs = str2double(EEG1.info.nominal_srate);
gap = max(1,round(0.7*fs));

for f = 1:numel(txt)

    % ---------- read ----------
    T = readtable(fullfile(txt(f).folder, txt(f).name));
    T.(T.Properties.VariableNames{1}) = str2double(string(T.(T.Properties.VariableNames{1})));
  
    % ---------- basic columns ----------
    T = table( ...
        T.Var1/1000*fs, ...
        T.Var14, ...
        T.Var10, ...
        'VariableNames',{'samples500','bpm','beeps'});

    % ---------- rhythm ----------
    b = T.beeps>0;
    T.rhythm = movmin(movmax(b,[gap 0]),[0 gap]);
    T.rhythm(1) = T.rhythm(2);

    % ---------- keep only state changes ----------
    change = [false; diff(T.rhythm)~=0 | diff(T.bpm)~=0];
    change(1)=true; change(end)=true;

    Tabv = T(change,{'samples500','bpm','rhythm'});
    Tabv.Properties.VariableNames{1} = 'beg_sample';
    % ---------- build intervals ----------

    Tabv.end_sample = [Tabv.beg_sample(2:end)-10; 0];
    Tabv(end,:) = [];
    Tabv.duration = (Tabv.end_sample-Tabv.beg_sample)/fs;
    Tabv.section  = repmat(f,height(Tabv),1);
    Tabv = movevars(Tabv, 'end_sample', 'After', 'beg_sample');
    Tabv.beg_sample(Tabv.beg_sample==0) = 1;
    % ---------- keep only real blocks ----------
    txtFiles{f} = Tabv(Tabv.duration>1.5,:);

end
%%
idxBeg  = find(strcmp(markers.time_series, 'Background Noise start'));
idxEnd = find(strcmp(markers.time_series, 'Background Noise stop'));
% convert marker time to EEG samples 
sampBeg  = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps(idxBeg));
sampEnd = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps(idxEnd));
trlTable = table(sampBeg', sampEnd',...
    zeros(numel(sampBeg),1), (1:numel(sampBeg))',...
    'VariableNames',{'beg_sample', 'end_sample', 'offset', 'section'});
trlTable.duration = (trlTable.end_sample - trlTable.beg_sample)/fs;

trlTableNew = table([], [],...
    [], [],[],[],...
    'VariableNames',txtFiles{1}.Properties.VariableNames);

for i=1:height(trlTable)
    x = txtFiles{i};
    x.beg_sample = x.beg_sample+trlTable{i,"beg_sample"};
    x.end_sample = x.end_sample+trlTable{i,"beg_sample"};
    trlTableNew = vertcat(trlTableNew, x);
end
%% convert to ft data and redefine trials
% convert LSL stream to fieldtrip data
ftEEG = LSL2ft(EEG1);
% redifine trials
cfg = [];
cfg.trl = trlTable{:,:};
datTrl = ft_redefinetrial(cfg, ftEEG);
%clear EEG1

%% basic preproc
clc;
cfg            = [];
cfg.channel    = 'EEG';
cfg.detrend    = 'yes';
cfg.continuous = 'no';
cfg.hpfilter    = 'yes';
cfg.hpfreq      = 0.5;
cfg.reref      = 'yes';
cfg.refchannel = {'all'}; % A1 and A2 are the ear-clips
% Notch filters to remove 50Hz and harmonics
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase

%cfg.dftfilter = 'yes';
%cfg.dftfreq = [50 100]; % line noise removal

pEEG = ft_preprocessing(cfg, datTrl); %p(preprocessed)EEG

%% high amplitude artifact detection
cfg = [];
cfg.continuous                   = 'yes';
cfg.artfctdef.zvalue.channel     = 'all';
cfg.artfctdef.zvalue.cutoff      = 50;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.zscore      = 'yes';
cfg.artfctdef.zvalue.interactive = 'yes';
[cfg, z_artifact] = ft_artifact_zvalue(cfg, pEEG);
% reject atrifact 
cfg = []; 
cfg.artfctdef.reject            = 'nan';
cfg.artfctdef.visual.artifact   = z_artifact;
zEEG = ft_rejectartifact(cfg,pEEG);

%% Manual: View Data
cfg = [];
cfg.ylim  = [-30 30];
cfg.blocksize = 30;
man_artifact = ft_databrowser(cfg,zEEG)
save([env.paths.artifacts  ID '_man_artifact'], "man_artifact");
%load([env.paths.artifacts  ID '_man_artifact'], "man_artifact");
%% Manual: Remove Artifacts
cfg = []; 
cfg.artfctdef.reject           = 'nan';
cfg.artfctdef.visual.artifact = man_artifact.artfctdef.visual.artifact;
mEEG = ft_rejectartifact(cfg,zEEG);

%% Run ICA
cfg = [];
cfg.method  = 'runica';
cfg.channel = {'EEG'};
cfg.numcomponent = 20;
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
cfg.component = [1,2];
dat_after_ICA = ft_rejectcomponent(cfg, comp);

%% view the data again
cfg = [];
cfg.ylim  = [-30 30];
cfg.blocksize = 720;
man_artifact = ft_databrowser(cfg,dat_after_ICA)


%% save data after preproc
save([env.paths.cleanData ID '_clean.mat'], "dat_after_ICA");
disp(['Saved All Data!']);













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DUMP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% divide into trials based on the markers
begIdx  = find(strcmp(markers.time_series, 'Background Noise start'));
stopIdx = find(strcmp(markers.time_series, 'Background Noise stop'));

% convert TIME → SAMPLE
begSamp  = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps(begIdx));
stopSamp = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps(stopIdx));
allSamp  = arrayfun(@(t) nearest(EEG1.time_stamps, t), markers.time_stamps);
section = zeros(numel(markers.time_series),1);


for k = 1:numel(begIdx)
    section(begIdx(k):stopIdx(k)) = k;
end

trltable = table( ...
    allSamp', ...
    zeros(numel(markers.time_stamps),1), ...
    zeros(numel(markers.time_stamps),1), ...
    section,...
    markers.time_series(:),...
    'VariableNames', {'beg_sample', 'end_sample', 'offset', 'section','type'} );

% make every 2nd ection the end sample of the previous one.
trltable.end_sample(1:end-1) = trltable.beg_sample(2:end);

% find the bpm values from the vector
bpm    = bpmInfo.time_series(13,:);
bpm    = bpm([true, diff(bpm) ~= 0]);
% remove any number that is not in the bpm values (112/129)
bpm    = bpm(find(bpm == 112 | bpm == 129));
% add 'no bpm' conditions
bpm = reshape([bpm; zeros(size(bpm))], 1, []);

trltable.bpm = [0 bpm(1:end-1)]';
%
% remove unwanted lines
trltable = trltable(1:end-1,:);
trltable = trltable(trltable.beg_sample ~= trltable.end_sample,:);
% fix overlapping problem caused by one trial ending when other begins.
trltable.beg_sample(2:end) = trltable.beg_sample(2:end)+1;
%trltable.type = [];

% add block duration (sec)
trltable.duration = (trltable.end_sample - trltable.beg_sample)/EEG1.info.effective_srate;