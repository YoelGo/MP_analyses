function [itc] = ITPC_function(Lcfg, trlDat)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% --- defaults ---
if ~isfield(Lcfg, 'foi')
    Lcfg.foi = 1:1:40;   % example default
end

if ~isfield(Lcfg, 'pad')
    Lcfg.pad = 10;
end

if ~isfield(Lcfg, 'channel')
    Lcfg.channel = 'EEG';
end

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'fourier';
cfg.taper     = 'hanning';     % or 'dpss' if you prefer multitaper smoothing
cfg.pad       = Lcfg.pad;
cfg.foi       = Lcfg.foi;      % frequencies of interest
cfg.channel   = Lcfg.channel;
freq          = ft_freqanalysis(cfg, trlDat);

% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition

itc           = [];
itc.label     = freq.label;
itc.freq      = freq.freq;
%itc.time      = freq.time;
itc.dimord    = 'chan_freq_time';

F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);           % number of trials

% compute inter-trial phase coherence (itpc)
itc.itpc      = F./abs(F);         % divide by amplitude
itc.itpc      = sum(itc.itpc,1);   % sum angles
itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension

% compute inter-trial linear coherence (itlc)
itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension

end