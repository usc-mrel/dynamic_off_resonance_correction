function [kernel] = gen_2d_kernel(wsize, shape) 

    if shape == 1 % rectangular shape
        kernel = ones(wsize, wsize);
        kernel = kernel/sum(kernel(:));

    elseif shape == 2 % circle 
        kernel = ones(wsize,wsize); 
        temp = hanning(wsize);
        w = temp*temp';
        kernel = kernel.*w;
        kernel = kernel./sum(kernel(:));
    else
        disp('shape input should be either of 1 (rectagular) or 2 (circle)');
        exit;
    end
end