function [pyr_enhance]=enhancement(N,lap_pyr,M,p,g)
pyr_enhance=cell(N,1);
for x = 1:N-1
    for y = 1:size(lap_pyr{x},1)
        for z=1:size(lap_pyr{x},2)
            if abs(lap_pyr{x}(y,z))<=M 
                pyr_enhance{x}(y,z)=g(x)*lap_pyr{x}(y,z)*...
                (1-(abs(lap_pyr{x}(y,z))/M))^p(x)+lap_pyr{x}(y,z);
            else
                pyr_enhance{x}(y,z)=lap_pyr{x}(y,z);
            end
        end
    end
end 