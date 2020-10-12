function  res = MCNUFFTSVD_FC(k,w,b1, dfm, para)

% Multicoil NUFFT operator 
% Based on the NUFFT toolbox from Jeff Fessler and the single-coil NUFFT
% operator from Miki Lustig
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng & Ricardo Otazo, NYU, 2012

% Add Field Inhomogeniety Correction (DOI: 10.1002/mrm.27373)
% Yongwan Lim (yongwanl@usc.edu), USC, 2017
%
% dfm : dynamic field map 

Nd = size(b1(:,:,1));
Jd = [6,6];
Kd = floor([Nd*2]);
n_shift = Nd/2;
for tt=1:size(k,3),
	kk=k(:,:,tt);
	om = [real(kk(:)), imag(kk(:))]*2*pi;
	res.st{tt} = nufft_init(om, Nd, Jd, Kd, n_shift,'kaiser');
end
res.adjoint = 0;
res.imSize = size(b1(:,:,1));
res.dataSize = size(k);
res.w = sqrt(w);
res.b1 = b1;
res.para = para;


nt = size(dfm,3);
res.mask = zeros(res.imSize(1), res.imSize(2), para.L, nt);
res.pmask = zeros(para.Nk, size(k,2), para.L, nt);

% Approximate the off-resonance exponential term
for tt=1:nt
%     disp(tt);
    fmap = dfm(:,:,tt);
    E_mat = exp(-1i*2*pi*para.t'*fmap(:)');
    [count, center] = hist(fmap(:),para.n_hist_bins);
    ek = exp(-1i*2*pi*para.t'*center);
    [U, ~, ~] = svd(repmat(sqrt(count), [size(ek,1) 1]).*ek, 'econ');
    B_mat = U(:,1:para.L);
    C_mat = pinv(B_mat)*E_mat;
    C_mat = reshape(C_mat, [para.L res.imSize(1), res.imSize(2)]);
   
    res.mask(:,:,:,tt) = permute(C_mat,[2 3 1]);

    for ff=1:para.L
        res.pmask(:,:,ff,tt) = repmat(B_mat(:,ff), [1,size(k,2)]);
    end
end

res = class(res,'MCNUFFTSVD_FC');

