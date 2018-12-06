% ALL_rhythmGLMfit.m
%
% Batch script to run rhythmGLMfit.m on multiple sessions

%% 
cfg = [];
cfg.writeOutput = 1;
cfg.plotOutput = 0;
cfg.output_dir = 'C:\temp'; % store files here
cfg.output_prefix = 'R0_'; % prefix filenames with this (identify runs)
cfg_master.nMinSpikes = 200;

cfg.nPleats = 2; % number of pleats (cross-validation runs) per cell
cfg.kFold = 2; % folds per pleat

%%
please = [];
please.rats = {'R117', 'R119', 'R132'};
fc = getDataPath(please);

%%
for iS = 1:length(fc) % for each session...
    
    pushdir(fc{iS});
    
    rhythmGLMfit(cfg); % do the business
    
    popdir;
    
end % of sessions