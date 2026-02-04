%% loading in LSL file and preprocessing the data
clc; clear all; close all;
cd('C:\Users\yoelgo\Documents\GitHub\MP_analyses\');
addpath('C:\Users\yoelgo\Documents\GitHub\MP_analyses\extra_functions\');
env = setupEnviroment11();

%%
dat = load([env.data.cleanFiles.folder '\' env.data.cleanFiles.name]);
dat = dat.dat_after_ICA;

%%
cfg = [];
cfg.winlen = 5;
windowDat = NaN_windowing(cfg, dat);
%%
% --- original power spectrum (keep trials) ---
cfg = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foi         = env.foi;
cfg.pad         = 'nextpow2';
cfg.keeptrials  = 'no';
subj_freq = ft_freqanalysis(cfg, data);