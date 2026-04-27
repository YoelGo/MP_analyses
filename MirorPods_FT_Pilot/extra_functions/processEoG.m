function [pEoG] = processEoG(EoGdat)
% process EoG data. find the polarities and reference to each other.
% Input: EoG channels.

%%% sort out eye channels
EoGdat.label{1} = 'eogH1';
EoGdat.label{2} = 'eogH2';
EoGdat.label{3} = 'eogV1';
EoGdat.label{4} = 'eogV2';

%extracting EOG signals from horizontal sensors
cfg = [];
cfg.channel = {'eogH1' 'eogH2'};
cfg.refref = 'yes';
cfg.demean = 'yes';
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase
cfg.refchannel = {'eogH1'};
eogH = ft_preprocessing(cfg,EoGdat);
%keep only one channel and rename to eogH
cfg = [];
cfg.channel = {'eogH2'};
eogH = ft_selectdata(cfg,eogH);
eogH.label = {'eogH'};

%extracting EOG signals from vertical sensors
cfg = [];
cfg.channel = {'eogV1' 'eogV2'};
cfg.refref = 'yes';
cfg.demean = 'yes';
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [49 51; 99 101];
cfg.bsfilttype  = 'but';
cfg.bsfiltord   = 3; % 3rd order
cfg.bsfiltdir   = 'twopass'; % zero-phase
cfg.refchannel = {'eogV1'};
eogV = ft_preprocessing(cfg,EoGdat);
%keep only one channel and rename to eogH
cfg = [];
cfg.channel = {'eogV2'};
eogV = ft_selectdata(cfg,eogV);
eogV.label = {'eogV'};

pEoG = ft_appenddata([], eogH, eogV);
end