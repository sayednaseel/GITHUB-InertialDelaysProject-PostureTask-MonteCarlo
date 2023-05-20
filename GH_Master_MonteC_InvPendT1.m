%Master_MonteC_InvPendT1
% Code to pass InitVal to MonteC_InvPend to get confidence intervals on
% the coefficient and exponent of the inertial delay scaling law for the
% inverted pendulum. Doing only Froude perturbations. When using the 95% confidence interval lower bounds,
%i.e. longest limbs and weakest muscle, the strength limit is 0.0483 Fr
%T1 uses fitlm type 1 regression data
% see ID doc 5-0-2
% see ID paper drafting doc 3-0 Posture task Type 1 regress graphs

clear all;close all;clc;
%% Loading data
graph=1;
%Niter=1e4;% Inertial delays paper
Niter=2;% Test
global g
g=-9.8066;

% scaling data
load('forelimbCIcalcv3','mdl_Llimb');mdl_LlimbF=mdl_Llimb;clear mdl_Llimb;
load('hindlimbCIcalcv3','mdl_Llimb');mdl_LlimbH=mdl_Llimb;clear mdl_Llimb;

% Muscle data for hindlimb
%Alexander, R. M., Jayes, A. S., Maloiy, G. M. O., & Wathuta, E. M. (1981). Allometry of the leg muscles of mammals. Journal of Zoology, 194(4), 539–552. https://doi.org/10.1111/j.1469-7998.1981.tb04600.x
% Muscle mass and fiber length from Table I: Ankle Extensors, All non hoppers
% Moment arm data from Table II: Ankle extensors, All non hoppers.
% Alexander et.al. 1981 reports the 95% confidence intervals on the coefficient "a" and exponent "b".
%The authors give the the 95% confidence intervals as a */ Range value and b +- Range value. 
% Here I find the standard error for a and b as: 
%the whole range of the 95% confidence interval/(2*TcritB) 
n2=33;%
TcritB=tinv(1-.025,n2-2);% Critical value for  95% confidence interval from a 33-2 DOF two tailed T distribution. Assuming DOF to be n-2 instead of n-1, since a and b required to perform a linear fit. 
% Muscle Mass
PowerL.aMmusc=5.1/1000;% Kg
PowerL.bMmusc=0.97;% 
PowerL.aMmuscRange=1.16;
PowerL.bMmuscRange=0.05;
%Note: 13th oct 2017,rearranging to match the strcture of coefCI
PowerL.CIMmusc=[PowerL.aMmusc/PowerL.aMmuscRange PowerL.aMmusc*PowerL.aMmuscRange;PowerL.bMmusc-PowerL.bMmuscRange PowerL.bMmusc+PowerL.bMmuscRange];
PowerL.bMmuscSE=(PowerL.CIMmusc(2,2)-(PowerL.CIMmusc(2,1)))/(2*TcritB);
PowerL.aMmuscSE=(PowerL.CIMmusc(1,2)-(PowerL.CIMmusc(1,1)))/(2*TcritB);
% Fiber length
PowerL.aFLmusc=10.6/1000;% meters
PowerL.bFLmusc=0.14;
PowerL.aFLmuscRange=1.18;
PowerL.bFLmuscRange=0.08;
PowerL.CIFLMusc=[PowerL.aFLmusc/PowerL.aFLmuscRange PowerL.aFLmusc*PowerL.aFLmuscRange;PowerL.bFLmusc-PowerL.bFLmuscRange PowerL.bFLmusc+PowerL.bFLmuscRange];
PowerL.bFLMuscSE=(PowerL.CIFLMusc(2,2)-(PowerL.CIFLMusc(2,1)))/(2*TcritB);
PowerL.aFLMuscSE=(PowerL.CIFLMusc(1,2)-(PowerL.CIFLMusc(1,1)))/(2*TcritB);
% Moment arm 
PowerL.aAMusc=9.4/1000; % meters
PowerL.bAMusc=0.38;
PowerL.aAMuscRange=1.07;
PowerL.bAMuscRange=0.03;
PowerL.CIAMusc=[PowerL.aAMusc/PowerL.aAMuscRange  PowerL.aAMusc*PowerL.aAMuscRange;PowerL.bAMusc-PowerL.bAMuscRange PowerL.bAMusc+PowerL.bAMuscRange ];
PowerL.bAMuscSE=(PowerL.CIAMusc(2,2)-(PowerL.CIAMusc(2,1)))/(2*TcritB);
PowerL.aAMuscSE=(PowerL.CIAMusc(1,2)-(PowerL.CIAMusc(1,1)))/(2*TcritB);


AnkleExt.TABLE(1,:)=[PowerL.aMmusc,PowerL.aMmuscRange,PowerL.CIMmusc(1,:),PowerL.bMmusc,PowerL.bMmuscRange,PowerL.CIMmusc(2,:)];
AnkleExt.TABLE(2,:)=[PowerL.aFLmusc,PowerL.aFLmuscRange,PowerL.CIFLMusc(1,:),PowerL.bFLmusc,PowerL.bFLmuscRange,PowerL.CIFLMusc(2,:)];
AnkleExt.TABLE(3,:)=[PowerL.aAMusc,PowerL.aAMuscRange,PowerL.CIAMusc(1,:),PowerL.bAMusc,PowerL.bAMuscRange,PowerL.CIAMusc(2,:)];

%}
%% Running the simulations 

%{1
%InitVal2=-0.01:-0.01:-0.49;modenam='Froude perturbation';labelx='Froude number';% in degrees. ORIGINAL
InitVal2=-0.01:-0.1:-0.49;modenam='Froude perturbation';labelx='Froude number';% in degrees. Test

Exp=[-3,log10(0.005),-2,-1,0,1,2,3,log10(5000),4];
MonteC.DelayVals=zeros(Niter,length(Exp),length(InitVal2));


 % If the inverted pendulum falls beyong 90 deg, I consider it a bad
 % simulation. In dataset MonteC.DelayVals, I ignore the bad sims by
 % replacing its movement time with NaN within MonteC_InvPend, which is
 % then ignored by the Fitlm command. in dataset MonteC.DelayVals2, I
 % include the bad sim movmeent times. 
 tic
for n=1:length(InitVal2)
    
   [OP]=GH_MonteC_InvPendT1(InitVal2(n),Niter,Exp,mdl_LlimbF,mdl_LlimbH,PowerL);
    
    MonteC.DelayVals(:,:,n)=OP.Tdelay;% in this dataset, I have replaced bad sim movement times with NaNs, which are ignored by fitlm
    MonteC.a(n,:)=OP.amean;
    MonteC.astd(n,:)=OP.astd;
    MonteC.aCI(n,:)=OP.aCI;
    MonteC.b(n,:)=OP.bmean;
    MonteC.bstd(n,:)=OP.bstd;
    MonteC.bCI(n,:)=OP.bCI;  
    MonteC.badsims(n,:)=OP.badsims;
  
    MonteC.DelayVals2(:,:,n)=OP.Tdelay2;% in this dataset, I have included the bad sim movement times. 
    MonteC.a2(n,:)=OP.a2mean;
    MonteC.a2std(n,:)=OP.a2std;
    MonteC.a2CI(n,:)=OP.a2CI;
    MonteC.b2(n,:)=OP.b2mean;
    MonteC.b2std(n,:)=OP.b2std;
    MonteC.b2CI(n,:)=OP.b2CI;      
    
    disp(['Perturbation size :  ' num2str(InitVal2(n)) ' ' labelx]);
end

%}


%%
%{
t=datetime;
coderutime=toc;
notes={'Using Type 1 regression data, 10000 monte carlo sims from 0.01 to 0.49 Fr. '};
save('Master-MonteC-InvPendT1');
%}
%% Graphing
%clear all;close all;clc

%load('Master-MonteC-InvPendT1');% Use this once you finish the run till 0.49 Froude number
%load('Master-MonteC-InvPend');
%graph=1;


if graph==1
    
load('InvPendT1-Froude_pert.mat','PowerLaw','PowerLawG','strengthlimG','InitVal');% loading results based on fitlm Model 1 regression
PowerLawB=PowerLaw;InitValB=InitVal;
indG=find(InitVal>strengthlimG,1,'last');
clear PowerLaw InitVal;    
    
% load('InvPend2-Froude_pert.mat','PowerLaw','PowerLawG','strengthlimG','InitVal');% loading results based on Kilbourne Model II data 
% indG=find(InitVal>strengthlimG,1,'last');


Exponentmeans=[mean(MonteC.bCI(:,1)) mean(MonteC.b) mean(MonteC.bCI(:,2))];
% VARIATION OF COEFF AND EXPONENT WITH PERTURBATION SIZE-------------------
nam=['InvPend Exponent and coefficient variations with 95 % Confidence intervals'];
figure('name',nam)
subplot(2,1,1)
hold on;
plot(-InitValB,PowerLawB(:,1).*1000,'k-','LineWidth',2)
plot(-InitValB(1:indG),PowerLawG(1:indG,1).*1000,'g-','LineWidth',2)
plot(-InitVal2,MonteC.a(:,1).*1000,'r--','LineWidth',1)
plot(-InitVal2,MonteC.aCI(:,1).*1000,'r.-','LineWidth',1)
plot(-InitVal2,MonteC.aCI(:,2).*1000,'r.-','LineWidth',1)
plot(-InitVal2,MonteC.a2(:,1).*1000,'b--','LineWidth',1)
plot(-InitVal2,MonteC.a2CI(:,1).*1000,'b.-','LineWidth',1)
plot(-InitVal2,MonteC.a2CI(:,2).*1000,'b.-','LineWidth',1)
axis([ 0 max(-InitValB) 0 max(MonteC.a2CI(:,2).*1000)]);
xlabel(labelx)
ylabel('Coefficient (ms)')
grid on;
%title(nam);
%legend('Kilbourne Mean','Geometric Scaling','Monte Carlo mean','95% CI')
legend('Type 1 fitlm Scaling','Geometric Scaling','a-Ignore bad sims','95% CI','95% CI','a2-Consider bad sims','95% CI2','95% CI2')

%title('Confidence Intervals Posture Task')

subplot(2,1,2)
hold on;
plot(-InitValB,PowerLawB(:,2),'k-','LineWidth',2)
plot(-InitValB(1:indG),PowerLawG(1:indG,2),'g-','LineWidth',2)
plot(-InitVal2,MonteC.b(:,1),'r--','LineWidth',1)
plot(-InitVal2,MonteC.bCI(:,1),'r.-','LineWidth',1)
plot(-InitVal2,MonteC.bCI(:,2),'r.-','LineWidth',1)
plot(-InitVal2,MonteC.b2(:,1),'b--','LineWidth',1)
plot(-InitVal2,MonteC.b2CI(:,1),'b.-','LineWidth',1)
plot(-InitVal2,MonteC.b2CI(:,2),'b.-','LineWidth',1)
axis([ 0 max(-InitValB) 0 max(MonteC.b2CI(:,2))]);
xlabel(labelx)
ylabel('Exponent')
grid on;

nam='No of bad sims';
figure('name',nam)
plot(-InitVal2,MonteC.badsims,'k.-')
xlabel('Perturbation size (Froude number)')
ylabel('failed sims')
grid on;
title(nam)




nam='Model I  vs Monte Carlo mean';
figure('name',nam)
subplot(2,1,1)
hold on;
plot(-InitValB,PowerLawB(:,1).*1000,'c.-','LineWidth',2)
plot(-InitVal2,MonteC.a(:,1).*1000,'r--','LineWidth',1)
grid on;
title(nam)
legend('Model I-fitlm','Monte Carlo mean')
subplot(2,1,2)
hold on;
plot(-InitValB,PowerLawB(:,2),'c.-','LineWidth',2)
plot(-InitVal2,MonteC.b(:,1),'r--','LineWidth',1)
grid on;


nam=['ID paper fig - InvPend Exponent and coefficient variations with 95 % Confidence intervals'];
figure('name',nam)
subplot(2,1,1)
hold on;
plot(-InitVal2,MonteC.a(:,1).*1000,'r-','LineWidth',1)
plot(-InitVal2,MonteC.aCI(:,1).*1000,'r-','LineWidth',1)
plot(-InitVal2,MonteC.aCI(:,2).*1000,'r-','LineWidth',1)
axis([ 0 max(-InitValB) 0 max(MonteC.aCI(:,2).*1000)]);
ylabel('a (ms)')

subplot(2,1,2)
hold on;

plot(-InitVal2,MonteC.b(:,1),'r-','LineWidth',1)
plot(-InitVal2,MonteC.bCI(:,1),'r-','LineWidth',1)
plot(-InitVal2,MonteC.bCI(:,2),'r-','LineWidth',1)

axis([ 0 max(-InitValB) 0 0.5]);
ylabel('b')
xlabel('Dimensionless velocity')





end

disp(['Posture task mean exponent values ::' num2str(mean(MonteC.bCI(:,1))) ' ' num2str(mean(PowerLawB(:,2))) ' ' num2str(mean(MonteC.bCI(:,2)))   ])

%{
% VARIATION OF COEFF AND EXPONENT WITH PERTURBATION SIZE----Paper---------------
nam=['InvPend Exponent and coefficient variations with 95 % Confidence intervals'];
figure('name',nam)
subplot(2,1,1)
hold on;
% plot(-InitVal,PowerLaw(:,1).*1000,'k-','LineWidth',2)
% plot(-InitVal2,MonteC.a(:,1).*1000,'r--','LineWidth',1)
% plot(-InitVal2,MonteC.aCI(:,1).*1000,'r.-','LineWidth',1)
% plot(-InitVal2,MonteC.aCI(:,2).*1000,'r.-','LineWidth',1)
plot(-InitVal2,MonteC.a2(:,1).*1000,'k-','LineWidth',2)
plot(-InitVal(1:indG),PowerLawG(1:indG,1).*1000,'g-','LineWidth',2)

plot(-InitVal2,MonteC.a2CI(:,1).*1000,'r.-','LineWidth',2)
plot(-InitVal2,MonteC.a2CI(:,2).*1000,'r.-','LineWidth',2)
axis([ 0 max(-InitVal) 0 max(MonteC.a2CI(:,2).*1000)]);
xlabel(labelx)
ylabel('Coefficient (ms)')
grid on;
%title(nam);
%legend('Kilbourne Mean','Geometric Scaling','Monte Carlo mean','95% CI')
%legend('Kilbourne Scaling','Geometric Scaling','a-Ignore bad sims','95% CI','95% CI','a2-Consider bad sims','95% CI2','95% CI2')
legend('Measured scaling','Geometric Scaling','95% CI')

%title('Confidence Intervals Posture Task')

subplot(2,1,2)
hold on;
%plot(-InitVal,PowerLaw(:,2),'k-','LineWidth',2)
% plot(-InitVal2,MonteC.b(:,1),'r--','LineWidth',1)
% plot(-InitVal2,MonteC.bCI(:,1),'r.-','LineWidth',1)
% plot(-InitVal2,MonteC.bCI(:,2),'r.-','LineWidth',1)
plot(-InitVal2,MonteC.b2(:,1),'k','LineWidth',2)
plot(-InitVal(1:indG),PowerLawG(1:indG,2),'g-','LineWidth',2)
plot(-InitVal2,MonteC.b2CI(:,1),'r.-','LineWidth',2)
plot(-InitVal2,MonteC.b2CI(:,2),'r.-','LineWidth',2)
axis([ 0 max(-InitVal) 0 max(MonteC.b2CI(:,2))]);
xlabel(labelx)
ylabel('Exponent')
%}
%% COMPS presentation grapsh
%{
%close all;
% VARIATION OF COEFF AND EXPONENT WITH PERTURBATION SIZE----Paper---------------
nam=['COMPS presentation graph-InvPend'];
figure('name',nam)
subplot(2,1,1)
hold on;
% plot(-InitVal,PowerLaw(:,1).*1000,'k-','LineWidth',2)
plot(-InitVal2,MonteC.a(:,1).*1000,'r-','LineWidth',1)
plot(-InitVal2,MonteC.aCI(:,1).*1000,'r--','LineWidth',1)
plot(-InitVal2,MonteC.aCI(:,2).*1000,'r--','LineWidth',1)
% plot(-InitVal2,MonteC.a2(:,1).*1000,'k-','LineWidth',2)
% plot(-InitVal(1:indG),PowerLawG(1:indG,1).*1000,'g-','LineWidth',2)
% 
% plot(-InitVal2,MonteC.a2CI(:,1).*1000,'r.-','LineWidth',2)
% plot(-InitVal2,MonteC.a2CI(:,2).*1000,'r.-','LineWidth',2)
axis([ 0 max(-InitVal) 0 max(MonteC.a2CI(:,2).*1000)]);
%xlabel(labelx)
ylabel('Coefficient (ms)')
%grid on;
%title(nam);
%legend('Kilbourne Mean','Geometric Scaling','Monte Carlo mean','95% CI')
%legend('Kilbourne Scaling','Geometric Scaling','a-Ignore bad sims','95% CI','95% CI','a2-Consider bad sims','95% CI2','95% CI2')
legend('Mean','+ 95% CI','- 95% CI')

%title('Confidence Intervals Posture Task')

subplot(2,1,2)
hold on;
%plot(-InitVal,PowerLaw(:,2),'k-','LineWidth',2)
plot(-InitVal2,MonteC.b(:,1),'r-','LineWidth',1)
plot(-InitVal2,MonteC.bCI(:,1),'r--','LineWidth',1)
plot(-InitVal2,MonteC.bCI(:,2),'r--','LineWidth',1)
% plot(-InitVal2,MonteC.b2(:,1),'k','LineWidth',2)
% plot(-InitVal(1:indG),PowerLawG(1:indG,2),'g-','LineWidth',2)
% plot(-InitVal2,MonteC.b2CI(:,1),'r.-','LineWidth',2)
% plot(-InitVal2,MonteC.b2CI(:,2),'r.-','LineWidth',2)
%axis([ 0 max(-InitVal) 0 max(MonteC.b2CI(:,2))]);
xlabel('Dimensionless velocity')
ylabel('Exponent')
axis([ 0 max(-InitVal) 0  0.5]);

%grid on;
%}