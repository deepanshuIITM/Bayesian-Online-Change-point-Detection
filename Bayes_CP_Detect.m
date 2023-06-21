%% FINAL PROJECT 
% PARAMETER AND STATE ESTIMATION (CH5115)
% SUBMITTED BY: ED19D402 
% NAME: DEEPANSHU
% This function computes posterior matrix, runlength and finds out change
% points by taking inputs as data, initial prior values  

function [R,max_R,run_len] = Bayes_CP_Detect(Alpha,Beta,k,mu,lambda_cp,y)

H = 1/lambda_cp;
y = y(:,1:end); l = length(y); pred_p = [];
R = zeros(l+1,400); R(1,1) = 1; max_R(1) = 1;
for j = 1:length(y)
Alpha = [Alpha   (Alpha(j)+0.5)];
Beta = [Beta    (Beta(j)+(k(j)*(y(j)-mu(j))^2)/(2*(k(j)+1)))];
mu = [mu     ((k(j)*mu(j)+y(j))/(k(j)+1))];
k = [k (k(j)+1)];

pred_p = [pred_p  pred_probability(y(j),k,Alpha,Beta,mu,j)];
 
    R(j+1,2:j+1) = R(j,1:j).*exp(pred_p)*(1-H);
    R(j+1,1) = sum(R(j,1:j).*exp(pred_p)*H);
    R(j+1,:) = R(j+1,:)/sum(R(j+1,:));
    R(j+1,find(R(j+1,:)<10^-15))=0;
    %%
max_R(j+1) = R(j+1,j+1);
run_len = j;
if max_R(j+1)<=10^(-8) 
    break
end
end
end