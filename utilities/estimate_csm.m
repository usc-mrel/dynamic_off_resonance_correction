function csm =  estimate_csm(data, para)
%   [csm] = estimate_csm(data, para)
%
%   Estimates coil sensitivity maps from a set of temporal averaged 
%   and spatial low-pass filtered spiral images
%
%   INPUT:
%     - data.w      [samples, full_spirals]   : density compensation func
%     - data.kdata  [samples, spirals, coil]  : spiral k-space 
%     - data.kloc   [samples, spirals]        : k-space coordinate
%     - para
%
%   OUTPUT:
%     - csm     [x,y,coil]    : Relative coil sensitivity maps
%
%
% (c) Yongwan Lim (yongwanl@usc.edu) University of Southern California, 2018.

% Average k-space data over time
Nframes= length(para.NarmsFull:para.NarmsFull:para.Ntviews); % the number of time frame to be averaged
time_average_kspace = zeros(para.Nk,para.NarmsFull,para.Nc); % time averaged kspace

for i=1:para.NarmsFull
    idx = i:para.NarmsFull:para.Ntviews;
    temp = data.kdata(:,idx(1:Nframes),:);   
    time_average_kspace (:,i,:) = sum(temp,2)./Nframes; 
end

time_average_img = zeros([para.N, para.Nc]); % time averaged image
op_NUFFT = NUFFT(data.kloc(:,1:para.NarmsFull),data.w(:,1:para.NarmsFull), [0,0], para.N);

for cc =1:para.Nc
    time_average_img(:,:,cc) = (op_NUFFT'*( time_average_kspace(:,:,cc).*sqrt(data.w(:,1:para.NarmsFull))));
end

% Spatial low-pass filter
talr_img = zeros([para.N, para.Nc]); % time averaged spatial low resolution image
kernel = gen_2d_kernel(para.sizehanningfilter,2); % spatial filter window (Hanning window)
normalization = conv2(ones(para.N),kernel,'same'); % normalization mask

for cc=1:para.Nc
    talr_img(:,:,cc) = conv2(squeeze(time_average_img(:,:,cc)),kernel,'same')./normalization;
end

%
csm = ismrm_estimate_csm_mckenzie(talr_img,0);

