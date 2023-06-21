%% FINAL PROJECT 
% PARAMETER AND STATE ESTIMATION (CH5115)
% SUBMITTED BY: ED19D402 
% NAME: DEEPANSHU
% This function plots posterior of run length distribution 

function  plot_post(posterior)
post_prob1= 60+log10(posterior);
for i= 1:size(posterior,1)
    for j = 1:size(posterior,2)
        if post_prob1(i,j)==-inf
            post_prob1(i,j)=0;
        end
    end
end
figure,
I = mat2gray(-post_prob1);
J = imrotate(I,90);
imshow(J,[-0.1 0.15]); colorbar; 
title('Posterior of Run Length on Log Scale')
ylabel('Run Length');xlabel('Time (T)')
end