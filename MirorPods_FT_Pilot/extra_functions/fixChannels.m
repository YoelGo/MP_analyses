function [EEGchFixed] = fixChannels(cfg, EEGdat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove a whole channel or section of a channel 
% variables:
% cfg.badchannel = insert the bad channel(s) name(s) here
% cfg.section    = 'all' to fix all the channel or time section
% cfg.elec       = electrodes layout for interpolation (1020 usually)
% cfg.layout     = EEG electrode cap layout.
%
% example:
% cfg.badchannel = {'F6, 'Cz', 'F8'}
% cfg.section    = {'all', [12, 14], [12, 'last']}
%
% section: 
% all or specific segments from the manual view data {{n1, m1}, {n2, m2}}. 
% note: see the segment number in the data viewer before with the given 
% block size. you can also use 'last' to fix until the end.
%
% badch:
% one or more channels (e.g., {'P8'} or {'P8', 'Cz'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run this cell only if there are channels to remove
% A code I built which can fix a channel (based on adjacent channels) in
% either a specific trial or in all trials.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% cfg field checks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(cfg, 'badch')
    badchannel = cfg.badch;
    if any(~ismember(cfg.badch, EEGdat.label))
    error('One or more of the badchannel names are wrong!!!!!!');
    end
else
    warning('no channels being fixed!')
end

if isfield(cfg, 'section')
    section = cfg.section;
    if length(section) ~= length(cfg.badch) && ~strcmp(section,'all')
    error('Channels and sections are not in the same length.');
    end
end


if ~isfield(cfg, 'layout')
    error('no layout detected');
end

if ~isfield(cfg, 'elec')
    error('no electrode mapping detected');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

badchannel   = cfg.badch;
section      = cfg.section;
elec         = cfg.elec;
layout       = cfg.layout;
blockSz      = cfg.blocksize;

% include the first section selected by user
for i=1:size(section,2)
    if ~strcmp(section{i},'all')
        if section{i}{1} ~= 1
            section{i}{1} = section{i}{1} -1
        else
            continue
        end
    end
end

cfg = [];
cfg.layout      = layout;
cfg.method      = 'triangulation';
neighbours      = ft_prepare_neighbours(cfg, EEGdat);

for i=1:length(badchannel)
    segment = section{i};
    cfg = [];
    cfg.badchannel  = {badchannel{i}};
    cfg.neighbours  = neighbours;
    cfg.method      =  'spline';
    cfg.elec        = elec;
    if strcmp(segment, 'all')
        segment = {0, 'last'};
    end
    if strcmp(segment{2}, 'last')
        segment{2} = ceil(EEGdat.time{1}(end)/blockSz);
    end
    segment{1} = segment{1} * blockSz;
    segment{2} = segment{2} * blockSz;
    cfg2 = [];
    if length(badchannel)>1
        curCh = badchannel;
        curCh(i) = [];
    else
        curCh = badchannel;
    end
    cfg2.channel = ft_channelselection( ...
                    [{'all'}, strcat('-',curCh)], ...
                    EEGdat.label);
    cfg2.latency = [segment{1}, segment{2}];
    tmp = ft_selectdata(cfg2,EEGdat);
    tmp = ft_channelrepair(cfg, tmp);

    [~, t1] = min(abs(EEGdat.time{1} - segment{1}));
    [~, t2] = min(abs(EEGdat.time{1} - segment{2}));
    chan_idx = find(strcmp(badchannel{i}, EEGdat.label));
    EEGdat.trial{1}(chan_idx,t1:t2) =  tmp.trial{1}(chan_idx,:);
end

EEGchFixed = EEGdat;
end