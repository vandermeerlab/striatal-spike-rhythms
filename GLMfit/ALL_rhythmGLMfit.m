% ALL_rhythmGLMfit.m
%
% Batch script to run rhythmGLMfit.m on multiple sessions
restoredefaultpath;
addpath(genpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms'));

%% 
cfg = [];
cfg.writeOutput = 1;
cfg.plotOutput = 0;
cfg.output_dir = 'C:\temp\GLMfit'; % store files here
cfg.output_prefix = 'R0_'; % prefix filenames with this (identify runs)
cfg.Target = 'Striatum';
cfg.nMinSpikes = 100;

cfg.nPleats = 2; % number of pleats (cross-validation runs) per cell
cfg.kFold = 2; % folds per pleat

%%
please = [];
please.rats = {'R117', 'R119', 'R132'}; % testing on isidro
%please.rats = {'R117', 'R119', 'R131', 'R132'}; % vStr-only rats
% please.rats = {'R149', 'R152', 'R156', 'R159', 'R169', 'R170', 'R184', 'R192', 'R194'}; % vStr-HC
[cfg.fd, cfg.fd_extra] = getDataPath(please);

%%
for iS = 1:length(cfg.fd) % for each session...
    
    cfg.iS = iS;
    pushdir(cfg.fd{iS});
    
    rhythmGLMfit2(cfg); % do the business
    
    popdir;
    
end % of sessions