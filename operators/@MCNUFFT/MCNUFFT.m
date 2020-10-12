function  res = MCNUFFT(k,w,b1)

% Multicoil NUFFT operator
% Based on the NUFFT toolbox from Jeff Fessler and the single-coil NUFFT
% operator from Miki Lustig
% Input
% k: k-space trajectory
% w: density compensation
% b1: coil sensitivity maps
%
% Li Feng & Ricardo Otazo, NYU, 2012
    
Nd = size(b1(:,:,1));
Jd = [6,6];
% Kd = floor([Nd*1.5]);
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
% res.w = (w);
res.b1 = b1;
res = class(res,'MCNUFFT');

