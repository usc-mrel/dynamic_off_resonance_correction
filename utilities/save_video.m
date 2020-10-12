function save_video(fname, img, fmap, framerate)
%     outputVideo = VideoWriter(fname,'Uncompressed AVI');
    outputVideo = VideoWriter(fname,'Motion JPEG AVI');
    outputVideo.FrameRate = framerate;
    open(outputVideo);

    Nt = size(img, 3);
    max_intensity = 0.6*max(abs(img(:)));
    max_fmap = max(fmap(:));
    hf = figure(100);
%     pause;
    for tt=1:Nt    
        fig = subplot(2,1,1);
        img_temp = img(:,:,tt);
        imagesc(abs(img_temp)); caxis([0 max_intensity]); axis image off; colormap(fig,'gray'); 
        title(sprintf('Image w/o Correction     Image w/ Correction'));
        
        fig = subplot(2,1,2);
        imagesc(fmap(:,:,tt)); axis image off; colormap(fig,'jet'); 
        caxis([-max_fmap max_fmap]); colorbar;
        title(sprintf('Estimated Field Map \nTime frame %i',tt));

        drawnow;
        frame = getframe(hf);
        writeVideo(outputVideo,frame);
    end
    close(outputVideo);

    close (hf)
end