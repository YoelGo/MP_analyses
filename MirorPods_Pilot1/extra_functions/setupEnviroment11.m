function [env] = setupEnviroment11()

%%%%% these paths need to be adjusted %%%%%
dataAyelet          = 'Z:\'
yoelpath            = [dataAyelet 'Experiments\Yoel\'];
toolsPath           = [yoelpath 'tools\'];
analysis_dir        = [yoelpath 'MP_pilot_analysis\'];
fieldtrip_path      = [toolsPath 'fieldtrip-20250106\'];
LAVI_path           = [toolsPath 'LAVI\'];
%layout_path        = [toolsPath 'landaulab_layout_62_P9-10.mat'];
layout_path         = [toolsPath 'mbt32ch_layout.mat'];

expdir              = [dataAyelet 'Experiments\MP_JointAction\'];
datadir             = [expdir 'data\'];

env = [];

%%
%%%%% initiate fieldtrip %%%%%
env.paths.fieldtrip_path = fieldtrip_path;
addpath(fieldtrip_path);
addpath([fieldtrip_path 'external\xdf\']); % add specific LSL (.xdf) file functions
ft_defaults;

%%
%%%%% innitiate env (envelope) variable %%%%%
%%% set paths
env.paths.rawData   = datadir;
env.paths.cleanData = [analysis_dir 'preproc\clean\'];
env.paths.artifacts = [analysis_dir 'preproc\artifacts\'];
env.paths.pltOut    = [analysis_dir 'plots\'];
env.paths.LAVIpath        = LAVI_path;

%%% set data variables
%env.data.cleanID         = cellfun(@(x) regexprep(x.ID, '_.*', ''), env.data.dfLAVI, 'UniformOutput', false);
env.data.rawFiles        = dir(fullfile(env.paths.rawData, '**', '*.fif'));
env.data.rawFileNames    = {env.data.rawFiles.name};    
env.data.cleanFiles      = dir(fullfile(env.paths.cleanData, '**', '*.mat'));
env.data.cleanFileNames  = {env.data.cleanFiles.name};    

%%% set EEG variables
env.EEG.elec        = ft_read_sens([fieldtrip_path 'template\electrode\standard_1020.elc']);
env.EEG.fsample     = 512;

% load EEG layout
cfg          = [];
cfg.layout   = fullfile(layout_path);
env.EEG.lay  = ft_prepare_layout(cfg);

% set frequency of interest (foi)
foi = 10.^(-0.1:0.05:1.2);
must = [1.5, 4, 9.5];
for i=1:2
    for mi = 1:3
        ind = nearest(foi,must(mi));
        foi(ind) = must(mi);
    end
end
env.EEG.foi = foi;

%%
%%% set LAVI variables
env.cfgLAVI.foi   = env.EEG.foi;
env.cfgLAVI.lag   = 1.5;
env.cfgLAVI.width = 5;
env.cfgLAVI.fs    = env.EEG.fsample;

%%
%%% set plotting variables
env.plt.darkGrey     = [84, 90, 99]/255;
env.plt.lightGrey    = [200, 203, 208]/255;
env.plt.lightBlue    = [152, 207, 222]/255;
env.plt.lilac        = [185, 150, 197]/255;
env.plt.peach        = [243, 212, 184]/255;

addpath(env.paths.LAVIpath)
end
