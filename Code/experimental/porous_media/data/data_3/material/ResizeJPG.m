function [ ] = ResizeJPG(  )
%This function reads an jpg file and returns the hsv img

file1 = 'porosity.jpg';
file2 = 'permeability.jpg';
file3 = 'saturation.jpg';

tmp1_l = rgb2hsv(imread(file1));
tmp2_l = rgb2hsv(imread(file2));
tmp3_l = rgb2hsv(imread(file3));
tmp1_l = imresize(tmp1_l, [100 100]);
tmp2_l = imresize(tmp2_l, [100 100]);
tmp3_l = imresize(tmp3_l, [100 100]);
tmp1_m = imresize(tmp1_l, [50 50]);
tmp2_m = imresize(tmp2_l, [50 50]);
tmp3_m = imresize(tmp3_l, [50 50]);
tmp1_s = imresize(tmp1_l, [20 20]);
tmp2_s = imresize(tmp2_l, [20 20]);
tmp3_s = imresize(tmp3_l, [20 20]);

imwrite(tmp1_l,'porosity_large.jpg','jpg')
imwrite(tmp1_m,'porosity_medium.jpg','jpg')
imwrite(tmp1_s,'porosity_small.jpg','jpg')
imwrite(tmp2_l,'permeability_large.jpg','jpg')
imwrite(tmp2_m,'permeability_medium.jpg','jpg')
imwrite(tmp2_s,'permeability_small.jpg','jpg')
imwrite(tmp3_l,'saturation_large.jpg','jpg')
imwrite(tmp3_m,'saturation_medium.jpg','jpg')
imwrite(tmp3_s,'saturation_small.jpg','jpg')

subplot(3,3,1)
    pcolor(tmp1_l(:,:,3));
    shading interp
subplot(3,3,2)
    pcolor(tmp2_l(:,:,3))
    shading interp
subplot(3,3,3)
    pcolor(tmp3_l(:,:,3))
    shading interp
subplot(3,3,4)
    pcolor(tmp1_m(:,:,3));
    shading interp
subplot(3,3,5)
    pcolor(tmp2_m(:,:,3))
    shading interp
subplot(3,3,6)
    pcolor(tmp3_m(:,:,3));
    shading interp
subplot(3,3,7)
    pcolor(tmp1_s(:,:,3));
    shading interp
subplot(3,3,8)
    pcolor(tmp2_s(:,:,3))
    shading interp
subplot(3,3,9)
    pcolor(tmp3_s(:,:,3))
    shading interp

end

