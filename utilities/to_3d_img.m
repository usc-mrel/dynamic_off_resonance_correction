function result = to_3d_img (img, rot_angle, num_col, num_row)

result = [];

for bb=1:num_row
    concat_ver = [];
    for aa=1:num_col
        slice_idx = aa+(bb-1)*num_col;
        concat_ver = [concat_ver imrotate(squeeze(img(:,:,slice_idx)),rot_angle)];
    end
    result = [result; concat_ver];
end

end