% -------------------------------------------------------------------------------
% Implementation of dynamic off-resonance correction for spiral dynamic MRI 
%
% Reference: Y Lim, SG Lingala, S Narayanan, KS Nayak, Dynamic off-resonance
% correction for spiral real-time MRI of speech, Magn Reson Med, 2018;81:234-246. 
% https://doi.org/10.1002/mrm.27373
%
% (c) Yongwan Lim (yongwanl@usc.edu) University of Southern California, 2018.
% -------------------------------------------------------------------------------

clear all; close all; clc

addpath(genpath('./'))

% Load config file and raw data
run ('./reconParam.m');
data = load(param.rawDataPath);
data.kdata = data.kdata(:,param.Coils,:); % do not use coil element that has low SNR
data.kdata = permute(data.kdata,[1 3 2]);


%% Coil Sensivity Map Estimation 
% Estimate from a temporal averaged and spatial low resolution image
fprintf('Coil Map Estimation \n');

csm = estimate_csm(data, param);

% Plot coil sensitivity map
% csm_plot = to_3d_img(csm, 90, param.Nc, 1);
% figure,  
% fig=subplot(2,1,1); imagesc(abs(csm_plot)); colormap(fig,'gray'); axis image off; title('Magnitude of CSM');
% fig=subplot(2,1,2); imagesc(angle(csm_plot)); colormap(fig,'jet');axis image off; title('Phase of CSM');


%% Dynamic Field Map Estimation
fprintf('Field Map Estimation \n');

[dfm, mask, ] = estimate_dfm(data, csm, param);
% dfm = dfm.*mask; % masking 
% dfm (dfm <0 ) = 0; % Practically, correction may work better with this thresholding

% Plot estimated dynamic field map
% figure, colormap jet;
% for tt=1:param.Nt
%     imagesc(imrotate(dfm(:,:,tt),90)); 
%     axis image off; colorbar; caxis([-1/param.te/2 1/param.te/2]);    
%     title(sprintf('Estimated Field Map - Time frame %i',tt));
%     drawnow;
% end


%% Multi Frequency Interpolation Reconstruction
% (just for a simple comparison purpose only)
fprintf('Multi Frequency Interpolation Reconstruction \n');

[img_w_cor_mfi, ] = recon_mfi(data, csm, dfm, param);

% Plot reconstructed image using MFI method 
% figure, 
% for tt=1:param.Nt
%     temp = [imrotate(img_wo_cor_nufft(:,:,tt), 90) imrotate(img_w_cor_mfi(:,:,tt), 90)] ;
%     fig = subplot(2,1,1); 
%     imagesc(abs(temp)); 
%     axis image off; 
%     colormap(fig,'gray'); 
%     title(sprintf('Reconstructed Images - Time frame %i',tt));
% 
%     fig = subplot(2,1,2); 
%     imagesc( imrotate(dfm(:,:,tt),90) ); 
%     axis image off; 
%     colormap (fig,'jet'); 
%     caxis([param.mfi_freq_range(1) param.mfi_freq_range(end)]); 
%     colorbar;
%     title(sprintf('Estimated Field Map - Time frame %i',tt));
%     
%     drawnow;
% end


%% Preparation for CS Reconstruction
fprintf('Preparation for CS Reconstruction \n');

w = repmat(data.w(:,1), [1 size(data.kloc,2)]); % dcf

% rearrange k-space data into a time-series
for tt=1:param.Nt
    idx = (tt-1)*param.Ntstep+1:(tt-1)*param.Ntstep+param.Ntstep;
    kdatau2(:,:,:,tt) = data.kdata(:,idx,:);
    ku2(:,:,tt) = data.kloc(:,idx);
    wu2(:,:,tt) = w(:,idx);
end

coil_idx = 1:length(param.Coils);
time_idx = param.Nt;
frame_idx = 1:time_idx;
input_csm = csm(:,:,coil_idx);
input_dfm = dfm(:,:,frame_idx);
input_ku = ku2(:,:,frame_idx);
input_wu = wu2(:,:,frame_idx);
input_kdatau = kdatau2(:,:,coil_idx,frame_idx);


%% CS Reconstruction w/o field inhomogeneity correction
fprintf('CS Reconstruction without field inhomogeneity correction \n');

param.y = input_kdatau.*permute(repmat(sqrt(input_wu),[1 1 1 param.Nc]),[1 2 4 3]);
param.E = MCNUFFT(input_ku, input_wu, input_csm);

% zero-filled recon
cs_img_wo_cor = param.E'*param.y;

% rescale regularization parameter lambda w.r.t the zero-filled recon
param.lambda = param.lambda*max(abs(cs_img_wo_cor(:))); 

fprintf('\n Start CG iterations ... \n')
T0 = tic;

figure(100); colormap gray;
for n=1:param.noutite
    [cs_img_wo_cor, cost, l2cost, l1cost, finalcost] = CSL1NlCg(cs_img_wo_cor, param); 
    imshow(abs(imrotate(cs_img_wo_cor(:,:,end),90)),[]); title(sprintf('Iteration - %i',n)); drawnow;
end
cs_img_wo_cor = cs_img_wo_cor/max(abs(cs_img_wo_cor(:)));

param.comp_time_cs_wo_cor = toc(T0);
fprintf('\n Recon is done ...(%4.1fs) \n', param.comp_time_cs_wo_cor)


%% CS Reconstruction w field inhomogeneity correction
fprintf('CS Reconstruction with field inhomogeneity correction \n');

param.E = MCNUFFTSVD_FC(input_ku, input_wu, input_csm, input_dfm, param);
cs_img_w_cor=param.E'*param.y;

fprintf('\n Start CG iterations... \n')
T0 = tic;

figure(101); colormap gray;
for n=1:param.noutite
    [cs_img_w_cor, cost, l2cost, l1cost, finalcost] = CSL1NlCg(cs_img_w_cor, param); 
    imshow(abs(imrotate(cs_img_w_cor(:,:,end),90)),[]); title(sprintf('Iteration - %i',n)); drawnow;
end
cs_img_w_cor = cs_img_w_cor/max(abs(cs_img_w_cor(:)));

param.comp_time_cs_w_cor = toc(T0);
fprintf('\n Recon is done ...(%4.1fs) \n', param.comp_time_cs_w_cor)


%% save result
fprintf('Save result \n');
mkdir(param.outputDataDir)
save(fullfile(param.outputDataDir, 'result.mat'), 'cs_img_wo_cor', 'cs_img_w_cor', 'dfm', 'csm', 'param');

clear img_result
for tt=1:param.Nt
     img_result(:,:,tt) = [imrotate(cs_img_wo_cor(:,:,tt), 90) imrotate(cs_img_w_cor(:,:,tt), 90)] ;
     dfmap_result(:,:,tt) = imrotate(dfm(:,:,tt), 90); 
end

close all;
save_video(fullfile(param.outputDataDir, 'result_video.avi'), img_result, dfmap_result, param.fps);

