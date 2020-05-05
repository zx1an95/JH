function mind= mssim(org,sup,win)

% org: Original image
% sup: Super-resolution image
% win: Size of the window

[m,n,p]= size(org);
c1= 0.01;
c2= 0.03;
sz= floor(win/2);

mer= [];
for k= 1:p
    for i= sz+1:m-sz
        for j= sz+1:n-sz
            six= org(i-sz:i+sz,j-sz:j+sz,k);
            six= six(:);
            siy= sup(i-sz:i+sz,j-sz:j+sz,k);
            siy= siy(:);
            meux= mean(six);
            meuy= mean(siy);
            sigx= std(six);
            sigy= std(siy);
            sigxy= sum((six-meux).*(siy-meuy))/(numel(six)-1);
            er= ((2*meux*meuy+c1)*(2*sigxy+c2))/((meux^2+meuy^2+c1)*(sigx^2+sigy^2+c2));
            mer= [mer er];
        end
    end
end

mind= sum(mer)/(numel(mer));
