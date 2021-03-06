%filter
gaussian_filter = [.05, .25, .4, .25, .05];  
gaussian_filter = gaussian_filter'*gaussian_filter;
gaussian_filter2 = 4*gaussian_filter;
N=2;%total pyramid level
lap_pyr = cell(N,1);
Ori_image = uCSRImgSet(:,:,17);%original image
J = NLmeansfilter(Ori_image,5,2,25);%denoise

%construct laplacian pyramid
for n = 1 : N-1
  J_gauss = imfilter(J, gaussian_filter, 'conv', 'same', 'replicate');
  J_gauss_low = J_gauss(1:2:size(J_gauss,1)-1,1:2:size(J_gauss,2)-1); %downsample
  J_gauss_high0=zeros(size(J_gauss,1),size(J_gauss,2));%upsample
  J_gauss_high0(1:2:end,1:2:end)=J_gauss_low;%insert zeros(to even columns and rows)
  J_gauss_high = imfilter(J_gauss_high0, gaussian_filter2, 'conv', 'same', 'replicate');
  lap_pyr{n} = J-J_gauss_high;%laplacian pyramid
  J=J_gauss_low; 
end

%enhancement function
M=140;%upper limit
p=[1.2 0 0 0 0 0];%power for different pyramid level
g=[6 0 0 0 0 0];%gain factor for different pyramid level
pyr_enhanced=enhancement(N,lap_pyr,M,p,g);%enhancement
                
%reconstruction
lap_pyr{N}=J_gauss_low;%set the Gaussian pyramid to be the second level of Laplacian pyramid
out = lap_pyr{N};

for i = N-1 : -1 : 1
    out_up=zeros(size(lap_pyr{i},1),size(lap_pyr{i},2));
    out_up(1:2:end,1:2:end)=out;%upsample,insert zero(to even columns and rows)
    out_filt=imfilter(out_up, gaussian_filter2, 'conv', 'same', 'replicate');
    out=pyr_enhanced{i}+out_filt;%reconstruct
end

out=mat2gray(out);
out=out.^1.3;
%Gaussianfilter = fspecial('gaussian',9,0.5);
%out=imfilter(out, Gaussianfilter, 'conv', 'same', 'replicate');smoothing
figure;imshow([out,mat2gray(Ori_image)]);

%image quality assessment
%max_variation = max(max(mat2gray(Ori_image))) - min(min(mat2gray(Ori_image)));
%psnr_value = psnr(out,mat2gray(Ori_image),max_variation);
%mssim_value = mssim(out,mat2gray(Ori_image),8);
%ambe_value = abs(mean2(out)-mean2(mat2gray(Ori_image)));