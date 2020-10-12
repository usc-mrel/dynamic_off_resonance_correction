function [img_mfi, img_mf] = recon_mfi(data, csm, dfm, para)
%   [img_mfi, img_mf] = recon_mfi(data, csm, dfm, para)
%
%   Multi frequency interpolation reconstruction
%
%   INPUT:
%     - data.w      [samples, full_spirals]   : density compensation func
%     - data.kdata  [samples, spirals, coil]  : spiral k-space 
%     - data.kloc   [samples, spirals]        : k-space coordinate
%     - csm         [x,y,coil]                : Relative coil sensitivity maps
%     - dfm         [x,y,frames]              : dynamic field map estimate
%     - para
%
%   OUTPUT:
%     - img_mfi     [x,y,frames]              : multi-frequency reconstructed image
%     - img_mf      [x,y,frames, dem_freq]    : multi-frequency image
%
% (c) Yongwan Lim (yongwanl@usc.edu) University of Southern California, 2018.

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

% Demodulate raw k-space data at a set of frequencies defined by para.mfi.freqrange
ii=1;
for f=para.mfi_freq_range
    % demoudlate the data
    modulated_kdatau(:,:,:,:,ii) = kdatau.*repmat(exp(1i*2*pi*f*para.t'),1,para.NarmsFull,para.Nc,para.Nt);
    ii = ii+1;
end

% MCNUFFT to get the coil-composite image at different frequency
img_mf = zeros([para.N(1), para.N(2), para.Nt, length(para.mfi_freq_range)]);
dcf_sq = permute( repmat( sqrt(wu) , 1, 1, 1, para.Nc), [1 2 4 3]);
for ff=1:length(para.mfi_freq_range)
    T0 = tic;
    img_mf(:,:,:,ff) = op_MCNUFFT'*(modulated_kdatau(:,:,:,:,ff).*dcf_sq);
    T = toc(T0);
    fprintf('Inverse NUFFT at freq = %4.fHz: %4.2fs\n',para.mfi_freq_range(ff),T);
end

%
A = exp(1i*2*pi*para.t.'*para.mfi_freq_range);

invA = pinv(A);
img_mfi = zeros(para.N(1), para.N(2), para.Nt);
for ii=1:para.N(1)
    for jj=1:para.N(2)
        coef = invA * exp(1i*2*pi*para.t.'*squeeze(dfm(ii, jj, :)).');
        img_mfi(ii, jj, :) = diag(coef.' * squeeze(img_mf(ii, jj, :, :)).');             
    end
end
img_mfi = img_mfi./max(abs(img_mfi(:)));

end