function ress = mtimes(a,bb)

 if a.adjoint,
     % Multicoil non-Cartesian k-space to Cartesian image domain
     % nufft for each coil and time point
     for tt=1:size(bb,4),
         for ch=1:size(bb,3),
             for ff=1:a.para.L
                b = conj(a.pmask(:,:,ff,tt)).*bb(:,:,ch,tt).*a.w(:,:,tt);
                res_temp(:,:,ff) = conj(a.mask(:,:,ff,tt)).*reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize)),a.imSize(1),a.imSize(2));
             end
             res(:,:,ch,tt) = squeeze(sum(res_temp,3));
         end
     end
     % compensate for undersampling factor
%      res=res*size(a.b1,1)*pi/2/size(a.w,2);     
     % coil combination for each time point
     for tt=1:size(bb,4),
%          ress(:,:,tt)=sum(res(:,:,:,tt).*conj(a.b1),3)./sum(abs((a.b1)).^2,3); %#ok<AGROW>
         ress(:,:,tt)=sum(res(:,:,:,tt).*conj(a.b1),3); %#ok<AGROW>
     end
 else
     % Cartesian image to multicoil non-Cartesian k-space 
     for tt=1:size(bb,3),
         for ch=1:size(a.b1,3)
             for ff=1:a.para.L
                res= (a.mask(:,:,ff,tt)).*bb(:,:,tt).*a.b1(:,:,ch); %#ok<AGROW>
                res2(:,:,ff) = (a.pmask(:,:,ff,tt)).*reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);
             end
             ress(:,:,ch,tt) = squeeze(sum(res2,3));
         end
     end
 end

