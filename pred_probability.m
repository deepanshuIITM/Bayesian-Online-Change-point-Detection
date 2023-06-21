%% FINAL PROJECT 
% PARAMETER AND STATE ESTIMATION (CH5115)
% SUBMITTED BY: ED19D402 
% NAME: DEEPANSHU
% This function computes posterior probability of new datapoint 

function pred_prob = pred_probability(x,k,alpha,Beta,mu,t)
kappa = alpha(t)*k(t)/(Beta(t)*(k(t+1)));
dof = 2*alpha(t);
pred_prob = -log(beta(0.5,0.5*dof)) + 0.5*log(kappa/dof)-...
    (0.5*(dof + 1))*log(1 + (kappa*((x-mu(t))^2))/(dof));
%   pred_prob = log(pdf('tLocationScale',x,mu(t),1/sqrt(kappa),dof));
end