function [EEGchFixed] = fixChannels(cfg, EEGdat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove a whole channel or section of a channel 
% variables:
% cfg.ch      = insert the bad channel(s) name(s) here
% cfg.section = 'all' to fix all the channel or time section
% cfg.elec    = electrodes layout for interpolation
%
% you can enter multiple channels to fix at once, for example:
% cfg.badch   = {'F6, 'Cz'}
% cfg.section = {'all', [12, 14]}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run this cell only if there are channels to remove
% A code I built which can fix a channel (based on adjacent channels) in
% either a specific trial or in all trials. 
badch   = cfg.badch;      % insert the bad channel(s) name(s) here
section = cfg.section;
elec    = cfg.elec;
if ~isempty(badch)
    cfg = [];
    cfg.layout      = env.lay ;
    cfg.method      = 'triangulation';
    neighbours      = ft_prepare_neighbours(cfg, mEEG);
    
    cfg = [];
    cfg.badchannel  = badch;
    cfg.neighbours  = neighbours;
    cfg.method      = 'spline';
    %cfg.elec       = ft_read_sens('C:\Users\fire-\OneDrive - rus10\#LIFE\MA\EEG_Experiment\fieldtrip-20230522\template\electrode\standard_1020.elc')
    cfg.elec        = elec;
    if trl ~= 'all'
        cfg2 = []
        cfg2.trials = trl
        tmp = ft_selectdata(cfg2,mEEG)
        tmp = ft_channelrepair(cfg, tmp);
        mEEG.trial{trl} =  tmp.trial{1}
    else
        mEEG = ft_channelrepair(cfg, mEEG);
    end
end
end