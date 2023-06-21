%% FINAL PROJECT 
% PARAMETER AND STATE ESTIMATION (CH5115)
% SUBMITTED BY: ED19D402 
% NAME: DEEPANSHU
%% QUES NO : 1
clc
clear all
close all

%% Import data
load('NMRlogWell')

% initialization 
RL(1) = 0; i=0; len = length(y); post_prob = zeros(len+1,400);
while(sum(RL)<=len-1)
    i=i+1;
    y=y(RL(i)+1:end);
    lambda_cp = 250; Alpha = 20; Beta = 2; mu = 1.15; k = 0.01;
  
    [R,max_R,run_len] = Bayes_CP_Detect(Alpha,Beta,k,mu,lambda_cp,y);
    RL(i+1) = run_len;
    fprintf('\n Change point detected at t=%d',sum(RL));
    data = y(1:RL(i+1));
    t = 1+sum(RL(1:i)):sum(RL(1:i))+length(data);
    post_prob(1+sum(RL(1:i)):sum(RL(1:i))+length(data)+1,:) = R(1:run_len+1,1:400);
    figure(1);
    subplot(2,1,1)
    plot(t,data),hold on
    xline(sum(RL),'r-.','LineWidth',1.5,'Label','CP')
    ylabel('data(y)');xlabel('Time (T)')
    title('Bayesian Online Change Point Detection')
    subplot(2,1,2)
    Mat = 50+log10(R(1:run_len+1,1:400));
    Y = linspace(1,400,400); X = linspace(sum(RL(1:i)),sum(RL),RL(i+1)+1);
    contour(X,Y,Mat'); colorbar; hold on
    title('Posterior of Run Length on Log Scale')
     ylabel('Run Length');xlabel('Time (T)')
end

%% Shades
plot_post(post_prob)
