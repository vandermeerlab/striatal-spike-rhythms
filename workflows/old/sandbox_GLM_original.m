%% setup

% paths: striatal-spike-rhythms

% load data
cd('C:\data\adrlab\R117-2007-06-20');
LoadExpKeys;

% spikes
S = LoadSpikes([]);

% LFP
cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
csc = LoadCSC(cfg);

% position
pos = LoadPos([]);

% reward deliveries
evt = LoadEvents([]);
reward_t = evt.t{2}; % should generalize this with known labels etc..

%% params
cfg_master = [];
cfg_master.dt = 0.001;

cfg_tcx = [];
cfg_tcx.bins = 50:10:600;

%% establish common timebase and other variables common across cells in this session
TVECe = ExpKeys.TimeOnTrack:cfg_master.dt:ExpKeys.TimeOffTrack; % edges
TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers

spd = getLinSpd([],pos);

C = cvpartition(length(TVECc),'KFold',10);

%% loop over all cells 
for iC = 1%length(S.t):-1:1
%iC = 2;

%%% dependent variable: binned spike train %%%
spk_binned = histc(S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);
spk_binned = spk_binned > 0; % binarize

%%% PREDICTOR: x- and y-location %%%

% compute tuning curve
cfg_tc = cfg_tcx;
cfg_tc.interp = 'linear';
[predicted_x, lambda_x, x_binned] = MakeTC_1D(cfg_tc, pos.tvec, pos.data(1,:), TVECc, spk_binned);

cfg_tc.bins = 0:10:480;
[predicted_y, lambda_y, y_binned] = MakeTC_1D(cfg_tc, pos.tvec, pos.data(2,:), TVECc, spk_binned);

%%% TODO: 2-D field predictor, requires 2-D tuning curve function %%%

%%% PREDICTOR: speed %%%
cfg_spd = cfg_tc;
cfg_spd.bins = 0:5:120;
[pred_spd, lambda_spd, spd_binned] = MakeTC_1D(cfg_spd, pos.tvec, spd.data, TVECc, spk_binned);

%%% PREDICTOR: LFP phase (e.g. gamma) %%%
cfg_phi = [];
cfg_phi.fpass = [6 10];
cfg_phi.dt = median(diff(csc.tvec));
cfg_phi.ord = 100;
cfg_phi.bins = -pi:pi/18:pi;
cfg_phi.interp = 'nearest';

phase = ComputePhase(cfg_phi,csc);
[theta_phi, lambda_phi, phi_binned] = MakeTC_1D(cfg_phi, csc.tvec, phase, TVECc, spk_binned);

cfg_phi.fpass = [1.5 4.5];
phase = ComputePhase(cfg_phi,csc);
[delta_phi, lambda_phi, phi_binned] = MakeTC_1D(cfg_phi, csc.tvec, phase, TVECc, spk_binned);

cfg_phi.fpass = [14 25];
phase = ComputePhase(cfg_phi,csc);
[beta_phi, lambda_phi, phi_binned] = MakeTC_1D(cfg_phi, csc.tvec, phase, TVECc, spk_binned);

cfg_phi.fpass = [45 65];
phase = ComputePhase(cfg_phi,csc);
[lowGamma_phi, lambda_phi, phi_binned] = MakeTC_1D(cfg_phi, csc.tvec, phase, TVECc, spk_binned);

cfg_phi.fpass = [70 90];
phase = ComputePhase(cfg_phi,csc);
[highGamma_phi, lambda_phi, phi_binned] = MakeTC_1D(cfg_phi, csc.tvec, phase, TVECc, spk_binned);

% TODO: threshold phase predictors to times when rhythm is actually
% present?

%%% PREDICTOR: conditional intensity function %%%
cif = conv(spk_binned,[0 0 0 1 1],'same'); % 2 ms refractory period
cif = (cif > 0);

%%% TODO: use cif estimated from acorr -- how to convert that to a
%%% predictor? simply convolve? (need to make asymmetric I suppose)
cfg_acf = []; cfg_acf.maxlag = 100; cfg_acf.binsize = cfg_master.dt;
cfg_acf.sided = 'onezero';
[acf,acf_tvec] = ComputeACF(cfg_acf,spk_binned);

cif_full = conv(spk_binned,acf,'same');

%%% PREDICTOR: time relative to reward delivery %%%
t_to_reward = ComputeTTR(reward_t,TVECc);
cfg_ttr = []; cfg_ttr.interp = 'nearest'; cfg_ttr.bins = [-5:0.1:5];
[ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, spk_binned);


%% implement models
%m1 = fitglm(cat(1,predicted_x,predicted_y,ttr,cif')', spk_binned', 'Distribution', 'poisson')
m2 = fitglm(cat(1,predicted_x,predicted_y,ttr,cif',cif_full',pred_spd)', spk_binned', 'Distribution', 'poisson')
m3 = fitglm(cat(1,predicted_x,predicted_y,ttr,cif',cif_full',pred_spd,delta_phi,theta_phi,beta_phi,lowGamma_phi,highGamma_phi)', spk_binned', 'Distribution', 'poisson')

OUT_tstat(iC,:) = m3.Coefficients.tStat;

end

%% some functions
function [predicted_x, lambda_x, x_binned] = MakeTC_1D(cfg, tvec, data, TVECc, spk_binned)
% construct tuning curve (mean value of spk_binned for binned values of
% data
%
% INPUTS:
%
% tvec, data: describe tuning variable (e.g. position, phase)
% TVECc, spk_binned: describe spike train
%
% cfg.interp: option to use for interpolation, 'linear', 'nearest' (for
%   things like phase)
% cfg.bins: bin edges (vector) for constructing tuning curves
%
% OUTPUTS:
%
% x_binned: binned version of data input
% lambda_x: tuning curve (size cfg.bins), gives spike probability as
%   function of binned tuning variable
% predicted_x: predicted value of tuning variable for each TVECc input
%   sample (simply looked up in lambda_x)

x_binned = interp1(tvec, data, TVECc, cfg.interp); % get variable of interest on common timebase 

[x_occ_hist, x_binned] = histc(x_binned, cfg.bins); % bin variable of interest
x_spk = x_binned(spk_binned); % variable of interest at spike times
x_spk_hist = histc(x_spk, 0.5:1:length(cfg.bins) + 0.5); x_spk_hist = x_spk_hist(1:end - 1);
lambda_x = x_spk_hist ./ x_occ_hist;

% based on tuning curve, get lambdas for each time bin
predicted_x = nan(size(x_binned));
keep = x_binned ~= 0;
predicted_x(keep) = lambda_x(x_binned(keep));

end

function phase = ComputePhase(cfg, csc)
% compute csc.data phase in frequency band cfg.fpass

fNQ = (1 ./ cfg.dt) / 2; % Nyquist frequency
Wn = cfg.fpass ./ fNQ;
b  = fir1(cfg.ord, Wn);

data = filtfilt(b, 1, csc.data);
phase = angle(hilbert(data));

end

function t_to_reward = ComputeTTR(reward_t,TVECc)
% find time to closest reward at each time point in TVECc

t_to_reward = TVECc-reward_t(1); % running total
for iR = 2:length(reward_t)
    
    this_dt = TVECc-reward_t(iR);
    swap_idx = find(abs(this_dt) < abs(t_to_reward)); % samples that are closer than running total
    t_to_reward(swap_idx) = this_dt(swap_idx);
end

end

%%
function [acf,tvec] = ComputeACF(cfg,spk_binned)

if isfield(cfg,'maxlag')
    [acf,tvec] = xcorr(spk_binned,spk_binned,cfg.maxlag);
else
    [acf,tvec] = xcorr(spk_binned,spk_binned);
end

tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;

if isfield(cfg,'sided')
    cut_idx = ceil(length(acf)/2);
    switch cfg.sided
        case 'one'
            acf = acf(cut_idx:end);
            tvec = tvec(cut_idx:end);
        case 'onezero'
            acf(1:cut_idx-1) = 0;
    end
end

end