function [dfm, mask, c_comp_img] =  estimate_dfm(data, csm, para)
%   [dfm, mask, coil_comp_img] =  estimate_dfm(data, csm, para)
%
%   Estimates dynamic field map from signle echo-time spiral data
%
%   INPUT:
%     - data.w      [samples, full_spirals]   : density compensation func
%     - data.kdata  [samples, spirals, coil]  : spiral k-space 
%     - data.kloc   [samples, spirals]        : k-space coordinate
%     - csm         [x,y,coil]                : Relative coil sensitivity maps
%     - para
%
%   OUTPUT:
%     - dfm         [x,y,frames]    : dynamic field map estimate
%     - mask        [x,y,frames]    : binary mask
%     - c_comp_img  [x,y,frames]    : coil composited dynamic image
%
% (c) Yongwan Lim (yongwanl@usc.edu) University of Southern California, 2018.

% smoothing kernels
kernel_s = gen_2d_kernel(3,2);
kernel_t = [1 2 1]/4;

crop_thresh = para.crop_thresh; % 1/para.te/2;

w = repmat(data.w(:,1), [1 size(data.kloc,2)]); % dcf

% rearrange k-space data into a time-series
for tt=1:para.Nt
    idx = (tt-1)*para.Ntstep+1:(tt-1)*para.Ntstep+para.NarmsFull;
    kdatau(:,:,:,tt) = data.kdata(:,idx,:);
    ku(:,:,tt) = data.kloc(:,idx);
    wu(:,:,tt) = w(:,idx);
end

% MCNUFFT operator to get coil-composite image series
op_MCNUFFT = MCNUFFT(ku,wu,csm);
c_comp_img = op_MCNUFFT'*(kdatau.*permute(repmat(sqrt(wu), 1, 1, 1, para.Nc), [1 2 4 3]));


% Smooth the coil-composited image
c_comp_img_s = zeros(size(c_comp_img));
c_comp_img_s_t = zeros(size(c_comp_img));

normalization_mask = convn(ones(para.N),kernel_s,'same'); % normalization mask
for tt=1:para.Nt
    c_comp_img_s(:,:,tt) = conv2(squeeze(c_comp_img(:,:,tt)),kernel_s,'same')./normalization_mask;
end

normalization_mask = convn(ones(para.Nt, 1),kernel_t','same'); % normalization mask
for ii=1:para.N(1)
    for jj=1:para.N(2)
        c_comp_img_s_t(ii,jj,:) = conv(squeeze(c_comp_img_s(ii,jj,:)),kernel_t,'same')./normalization_mask;
    end
end

% Estimate dynamic field map
dfm = angle(c_comp_img_s_t)/(-2*pi*para.te);

dfm(isnan(dfm)) = 0;
dfm(dfm>=crop_thresh) = crop_thresh;
dfm(dfm<=-crop_thresh) = -crop_thresh;

% Mask for field map 
mask= zeros(size(dfm));
thresh = 0.01;
for tt=1:para.Nt
    temp = abs(c_comp_img_s_t(:,:,tt)).^2;
    temp = temp/max(temp(:));    
    temp(temp<thresh) = 0;
    temp(temp>thresh) = 1;
    se = strel('disk',1);
    temp = imclose(temp,se);
    mask(:,:,tt) = temp;
end

c_comp_img = c_comp_img./max(abs(c_comp_img(:)));
