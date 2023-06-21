%% FINAL PROJECT 
% PARAMETER AND STATE ESTIMATION (CH5115)
% SUBMITTED BY: ED19D402 
% NAME: DEEPANSHU
%% QUES NO : 2
clc
clear all
close all

%% Load data
load('historic')
figure(1);
subplot(2,1,1); parcorr(y);ylabel('\phi_{yy}[l]')

%% LS estimates
sys = ar(y',1); % ar routine uses LS method
subplot(2,1,2)
resid(y',sys);ylabel('\sigma_{ee}[l]');
theta1 = -sys.Report.Parameters.ParVector  ;
var_v = sys.Report.Parameters.FreeParCovariance;  
sigma_e = sys.NoiseVariance; 


%% Recursive LS
load('new')
% Fit Model 1 
thetainit = theta1; Pinit = var_v; stop = 200;
obj_p1 = recursiveLS(1,'InitialParameters',thetainit,'InitialParameterCovariance',Pinit);
thetaest_vec = []; Ptheta = [];K=[]; i=0; RL(1) = 0;
for n = 2:numel(y)
    if n>stop+sum(RL)
        break
    end
    H = y(n-1);    
    [theta,~] = obj_p1(y(n),H);    
    Ptheta(n-1,1:length(thetainit)) = obj_p1.ParameterCovariance;     
    thetaest_vec(n-1,1:length(thetainit)) = theta; 
    K(n-1,1:length(thetainit)) = obj_p1.ParameterCovariance*H';
    if n>=stop
         y= y(n:end);
        while(i<1)
            i=i+1;
            y=y(RL(i)+1:end);
            lambda_cp = 50; Alpha = 10; Beta = 1; mu = 0; k = 1;
            [R,max_R,run_len] = Bayes_CP_Detect(Alpha,Beta,k,mu,lambda_cp,y);
            RL(i+1) = run_len;
            fprintf('\n Change point detected at time %d \n',n+sum(RL));
            data = y(1:RL(i+1));
            t = 1+sum(RL(1:i)):sum(RL(1:i))+length(data);
            figure(2); 
            subplot(3,1,1)
            xline(n+sum(RL),'r-.','LineWidth',1.5,'Label','CP'); hold on
            ylabel('data(y)');xlabel('Time (T)')
            title('Bayesian Online Change Point Detection');hold on
            subplot(3,1,3)
            Mat = -log(R(1:run_len+1,1:100));
            Y1 = linspace(1,100,100); X = n+linspace(sum(RL(1:i)),sum(RL),RL(i+1)+1);
            contour(X,Y1,Mat'); colorbar; hold on
            title('Posterior of Run Length on Log Scale')
            ylabel('Run Length');xlabel('Time (T)')
        end
    end
    load('new')
end
figure(2)
subplot(3,1,1); plot(thetaest_vec); title('Estimated Parameter \theta');
ylabel('\theta');xlabel('Time');grid on; hold on;
 
%% Model Update
y_new = y(stop+run_len:end);
% Fit Model 1 
thetainit1 = -0.4; Pinit1 = var_v;
obj_p2 = recursiveLS(1,'InitialParameters',thetainit1,'InitialParameterCovariance',Pinit1);
thetaest_vec1 = []; Ptheta1 = []; Kn = [];
for n = 2:numel(y_new)
    H1 = y_new(n-1);    
    [theta,~] = obj_p2(y_new(n),H1);    
    thetaest_vec1(n-1,1:length(thetainit1)) = theta;    
    Ptheta1 (n-1,1:length(thetainit1)) = obj_p2.ParameterCovariance; 
    Kn(n-1,1:length(thetainit1)) = obj_p2.ParameterCovariance*H1';
end

%% Plots
figure(2)
subplot(3,1,2)
plot(Ptheta); title('Parameter Variance \sigma_{\theta}'); ylabel('\sigma_\theta');
xlabel('Time');grid on; hold on;

figure(3)
subplot(2,1,1)
plot(thetaest_vec1); title('Updated Parameter \theta_{new}');
ylabel('\theta_{new}');xlabel('Time');grid on; hold on;
subplot(2,1,2)
plot(Ptheta1);title('Updated Parameter Variance \sigma_{\theta_{new}}');
ylabel('\sigma_{\theta_{new}}');xlabel('Time');grid on; hold on;

%% LS estimates
figure(4); 
subplot(2,1,1); parcorr(y_new);ylabel('\phi_{yy_{new}}[l]');
sys1 = ar(y_new',1); % ar routine uses LS method
subplot(2,1,2)
resid(y_new',sys1);ylabel('\sigma_{ee_{new}}[l]');
theta_new = -sys1.Report.Parameters.ParVector  ;
var_v_new = sys1.Report.Parameters.FreeParCovariance;  
sigma_e_new = sys1.NoiseVariance; 

%% Error Calculation
fprintf('Error Variance in base model= %d \n',sigma_e);
fprintf('Error Variance in updated model= %d \n',sigma_e_new);
