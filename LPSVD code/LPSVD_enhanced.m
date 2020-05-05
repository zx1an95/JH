Gaussianfilter=fspecial('gaussian',9,0.5);
gaussian_filter = [.05, .25, .4, .25, .05];  
gaussian_filter = gaussian_filter'*gaussian_filter;
gaussian_filter2=4.*gaussian_filter;
N=2;%total pyramid level
lap_pyr = cell(N,1);
O_image = uCSRImgSet(:,:,2);%original image
J=mat2gray(O_image);

%construct laplacian pyramid
for n = 1 : N-1
  J_gauss = imfilter(J, gaussian_filter, 'conv', 'same', 'replicate');
  J_gauss_low = J_gauss(1:2:size(J_gauss,1)-1,1:2:size(J_gauss,2)-1); 
  J_gauss_high0=zeros(size(J_gauss,1),size(J_gauss,2));%upsample
  J_gauss_high0(1:2:end,1:2:end)=J_gauss_low;%insert zeros
  J_gauss_high = imfilter(J_gauss_high0, gaussian_filter2, 'conv', 'same', 'replicate');
  lap_pyr{n} = J-J_gauss_high;%laplacian pyramid
  J=J_gauss_low; 
end

%SVD
GHE_gauss=histeq(mat2gray(J_gauss_low));
[U1,S1,V1]=svd(mat2gray(J_gauss_low));
[U2,S2,V2]=svd(GHE_gauss);
weighting_factor=max(max(S2))/max(max(S1));
NewS=0.5*(weighting_factor.*S1+1/weighting_factor.*S2);
New_J_gauss_low=U1*NewS*V1';
%enhancement function
M=0.04;
xc=0.015;
p=[1.2 0 0 0 0 0];%for different levels
g=[2 0 0 0 0 0];%for different levels 
pyr_enhanced=enhancement_func(N,lap_pyr,xc,M,p,g);
                
%reconstruction
lap_pyr{N}=New_J_gauss_low;
out = lap_pyr{N};
for i = N-1 : -1 : 1
    out_up=zeros(size(lap_pyr{i},1),size(lap_pyr{i},2));
    out_up(1:2:end,1:2:end)=out;%upsample,insert zero
    out_filt = imfilter(out_up, gaussian_filter2, 'conv', 'same', 'replicate');
    out = mat2gray(pyr_enhanced{i}) + out_filt;
end

out=mat2gray(out);
%different activate function
%out=1.*log((1+out)./(1-out));
%out=exp(out)./sum(sum(exp(out)));
%out=2.^out;
%out=1./(0.2+exp(-out));
out=out.^1.3;
%out=tanh(out);
%out=max(0.18,out);
out=mat2gray(out);
figure;imshow([out,mat2gray(O_image)]);
max_variation2 = max(max(mat2gray(O_image))) - min(min(mat2gray(O_image)));
psnr_value2 = psnr(out,mat2gray(O_image),max_variation2);
mssim_value4 = mssim(out,mat2gray(O_image),8);
ambe_value = abs(mean2(out)-mean2(mat2gray(O_image)));