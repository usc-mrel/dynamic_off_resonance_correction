%
% Parameter strcuture for Off-resonance Correction
%
% (c) Yongwan Lim (yongwanl@usc.edu) University of Southern California, 2018.

param = struct;

%% File parameters
param.rawDataPath       = fullfile(pwd, 'lac04012018_17_23_28.mat');
param.miscDataDir       = fullfile(pwd);
param.imgDataDir        = fullfile(pwd);
param.outputDataDir     = fullfile(pwd, 'output');

%% Sequence parameters
param.Nk = 630; % the number of samples per spiral interleave
param.Coils = [1 2 3 4 5 7 8];
param.Nc = length(param.Coils); % the number of coil elements used
param.Ntviews = 3000; % the number of TRs (spiral interleaves)
param.TR = 6.004*1e-3; % repetition time
param.te = 8e-4; % echo time
param.Ts = 4e-6; % sampling time interval which is a reciprocal of bandwidth  
param.N = [83 83]; % image size to be reconstructed 
param.tad = param.Nk*param.Ts; % readout duration
param.t = param.te + [0:param.Nk-1]*param.tad/param.Nk; % time map 

% the number of spirals per frame required for full-sampling (Nyquist)
param.NarmsFull = 13; 

% time frame of field map is updated every ntstep-TR
% whose number determines the frame rate of the final reconstructed image
param.Ntstep = 4;
param.fps = 1/(param.Ntstep*param.TR); % frame rate (for video generation)

% reconstruct T sec time series of image (for example, here T=5)
% # of time frames (synchronize with coil images)
param.Nt =  round(4.8*param.fps); 
% param.Nt = floor(param.Ntviews/param.NarmsFull); 

% check if Nt is out of bound or not
nt_max = floor(max((param.Ntviews-param.NarmsFull)/param.Ntstep + 1));
nt_max = nt_max - 1;

if param.Nt > nt_max
    param.Nt = nt_max;
    fprintf('Nt is out of bound, Nt is set to %d \n',nt_max);
end

%% coil estimation parameter
param.sizehanningfilter = 21; % filter size to smooth coil sensitivity map

%% Field map estimation parameter
param.center_freq_adjust = 0; 
param.mfi_freq_range = -625:125:625;
param.crop_thresh = 1/param.te/2;

%% Reconstruction paramters
% L in the forward models (corresponding to # of MCNUFFTs performed)
param.L = 6;
% the number of histogram bins in constructing forward model in image recon
param.n_hist_bins = 40;

% regularization paramter lambda for temporal finite difference
param.lambda = 0.2;
param.W = TV_Temp();

param.noutite = 4; % the number of outer iterations 
param.nite = 8; % the number of inner iterations 
param.display=1;
