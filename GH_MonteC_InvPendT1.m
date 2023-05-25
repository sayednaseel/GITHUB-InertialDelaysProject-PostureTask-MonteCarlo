function [OP]=GH_MonteC_InvPendT1(InitVal2,Niter,Exp,mdl_LlimbF,mdl_LlimbH,PowerL)

%Code to run Monte Carlo sims on the inverted pendulum. Scaling initial
%perturbation velocity based on Froude number, not constant angular
%velocity. 
% Inputs:: 
% InitVal2 is the Froude number.
% Niter is the number of iterations.
% mdl_LlimbF &H are models created using fitlm for the Kilbourne 2013 scaling data for forelimb and hindlimb length. 
% PowerL is the scaling information for the ankle extensor muscles from Alexander 1981
% 
% if a simualation fails, i.e. inverted pendulum falls below 90 degrees, OP.fail outputs a zero for that row. 
% see IDdoc 5-0-2, 5-0-3 and update 4-1 for more details. 
    
graph=0;
textshow=1;
global g

%% To run as script, comment out function definition on top and uncomment below. 
%{
clear all;close all;clc;
global g
graph=1;
textshow=1;
g=-9.8066;
Niter=1e3;% number of Monte carlo sims
Exp=[-3,log10(0.005),-2,-1,0,1,2,3,log10(5000),4];ind0=find(Exp==0);
%Inverted pendulum test. Initval is Initial angle. PLa and PLb are results
% from previous simulations (Master_InvPend.m) without considering the
% confidence intervals. 
Testmode='Froude perturbation';
InitVal2=-0.49;% Initval in Froude number
PLa=0.080819*1000;% Actual coeff for InitVal from Master_InvPendT1.m
PLb=0.38568;% Actual expoenent for InitVal 0.49

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
%}
%% Predefined matrices
M=10.^Exp;
nMusc=33;% Alexander used 33 species to calculate his muscle features scaling. 
dofMusc=nMusc-2;
%Initializing matrices
IPvals.LlimbF=zeros(Niter,length(M));% m
IPvals.LlimbH=zeros(Niter,length(M));% m
IPvals.Llimb=zeros(Niter,length(M));% m
IPvals.MOI=zeros(Niter,length(M));% m
IPvals.Vpert=zeros(Niter,length(M));
IPvals.Wpert=zeros(Niter,length(M));
IPvals.Mmusc=zeros(Niter,length(M));% muscle mass, kg
IPvals.FLmusc=zeros(Niter,length(M));% muscle fiber length, m
IPvals.Amusc=zeros(Niter,length(M));% muscle moment arm, m
IPvals.Tmusc=zeros(Niter,length(M));% muscle torque. Nm
Coeff.aMmusc=zeros(Niter,1);
Coeff.bMmusc=zeros(Niter,1);
Coeff.aFLmusc=zeros(Niter,1);
Coeff.bFLmusc=zeros(Niter,1);
Coeff.aAmusc=zeros(Niter,1);
Coeff.bAmusc=zeros(Niter,1);
%% Simulation


tic
for nn=1:Niter
% Inertial properties

IPvals.LlimbF(nn,:)=10.^random(mdl_LlimbF,Exp');
IPvals.LlimbH(nn,:)=10.^random(mdl_LlimbH,Exp');

IPvals.Llimb(nn,:)=(IPvals.LlimbF(nn,:)+IPvals.LlimbH(nn,:))./2;
IPvals.MOI(nn,:)=M.*IPvals.Llimb(nn,:).^2;% Point mass moment of inertia 

%Initial Conditions
    IPvals.Vpert(nn,:)=InitVal2.*sqrt(IPvals.Llimb(nn,:).*-g);
    IPvals.Wpert(nn,:)=IPvals.Vpert(nn,:)./IPvals.Llimb(nn,:);
    yinit=[zeros(length(M),1) IPvals.Wpert(nn,:)'];


%Muscle properties randomly sampled from a T distribution
Coeff.aMmusc(nn,1)=PowerL.aMmusc+trnd(dofMusc).*PowerL.aMmuscSE;
Coeff.bMmusc(nn,1)=PowerL.bMmusc+trnd(dofMusc).*PowerL.bMmuscSE;
IPvals.Mmusc(nn,:)=Coeff.aMmusc(nn,1)*M.^Coeff.bMmusc(nn,1);% mass of muscle (kg). Triceps All non hopper

Coeff.aFLmusc(nn,1)=PowerL.aFLmusc+trnd(dofMusc).*PowerL.aFLMuscSE;
Coeff.bFLmusc(nn,1)=PowerL.bFLmusc+trnd(dofMusc).*PowerL.bFLMuscSE;
IPvals.FLmusc(nn,:)=Coeff.aFLmusc(nn,1)*M.^Coeff.bFLmusc(nn,1);% fiber length muscle (m). Triceps All non hoppers
rho=1060;%density of muscle in kg/m^3. (Mendeys and keys 1960, Alexander 1977)
CSmusc=IPvals.Mmusc(nn,:)./rho./IPvals.FLmusc(nn,:).*(100^2);% Cross sectional area of muscle in cm^2
Fmusc=CSmusc.*20;% assuming 20N/cm^2 
Coeff.aAmusc(nn,1)=PowerL.aAMusc+trnd(dofMusc).*PowerL.aAMuscSE;
Coeff.bAmusc(nn,1)=PowerL.bAMusc+trnd(dofMusc).*PowerL.bAMuscSE;
IPvals.Amusc(nn,:)=Coeff.aAmusc(nn,1)*M.^Coeff.bAmusc(nn,1);% Moment arm value from Triceps, all non hoppers Table II
IPvals.Tmusc(nn,:)=4*Fmusc.*IPvals.Amusc(nn,:);% Torque by the muscle (Nm)    

%--------------------------------------------------------------------------
for i=1:length(M)
       
    [t,y,Tswitch]=GH_InvPendT1(IPvals.Tmusc(nn,i),IPvals.Llimb(nn,i),M(i),yinit(i,:));
    
% If the inverted pendulum has fallen beyond 90 degs, record a bad simulation
% and enter NaN as movement time. NaNs are ignored by fitlm. 
    if min(radtodeg(y(:,1)))>-90
        OP.Tdelay(nn,i)=t(end);
        OP.goodsim(nn,i)=1;
    elseif min(radtodeg(y(:,1)))<-90
        OP.Tdelay(nn,i)=NaN;
        OP.goodsim(nn,i)=0;
    end
    OP.Tswitch(nn,i)=Tswitch;
    OP.Tdelay2(nn,i)=t(end);%In Tdelay2, I include the bad sims
    
    if textshow==1
    disp(['Simulation number:: ' num2str(nn) ';Mass (kg):: ' num2str(M(i)) ';Torque (Nm):: ' num2str(IPvals.Tmusc(nn,i))]);
    end
end
%--------------------------------------------------------------------------
%% finding scaling of inertial delay
% OP.Tdelay, I ignore the bad simulations. 
mdl=fitlm(log10(M),log10(OP.Tdelay(nn,:)),'linear');% fitlm ignores NaNs 
OP.b(nn)=mdl.Coefficients.Estimate(2);
OP.loga(nn)=mdl.Coefficients.Estimate(1);

% in OP.Tdelay, I include the bad simulations
mdl2=fitlm(log10(M),log10(OP.Tdelay2(nn,:)),'linear');%  
OP.b2(nn)=mdl2.Coefficients.Estimate(2);
OP.loga2(nn)=mdl2.Coefficients.Estimate(1);

end
toc

%% Statistics on a and b of the inertial delay scaling 
OP.bmean=mean(OP.b);
OP.bstd=std(OP.b);
%OP.bSE=OP.bstd/sqrt(Niter);% since running monte carlo sims are like
%repeatedly sampling the data, and fitting a line is like taking a mean,
%you dont need to use the formula for estimated standard error "std(x)/sqrt(n)".
%The standard deviation of the sample mean distribution itself is the standard error
Tval=tinv([0.025 0.975],Niter-1);
OP.bCI=OP.bmean+Tval*OP.bstd;
OP.amean=mean(10.^OP.loga);
OP.astd=std(10.^OP.loga);
%OP.aSE=OP.astd/sqrt(Niter);
OP.aCI=OP.amean+Tval*OP.astd;


OP.b2mean=mean(OP.b2);
OP.b2std=std(OP.b2);
OP.b2CI=OP.b2mean+Tval*OP.b2std;
OP.a2mean=mean(10.^OP.loga2);
OP.a2std=std(10.^OP.loga2);
%OP.aSE=OP.astd/sqrt(Niter);
OP.a2CI=OP.a2mean+Tval*OP.a2std;



Badsimcheck1=sum(OP.goodsim,2);
Badsimcheck2=Badsimcheck1<length(M);
OP.badsims=sum(Badsimcheck2);% counts number of sims where muscles could not stop InvPend from falling beyond 90 degs
%% Saving dataset
%{
t=datetime;
Notes='MonteC_InvPendT1 run of 1000 sims after changing DOF to n-1, ID paper codes, Froude number from 0.01 to 0.49';
save MonteCsimT1;
%}

%% Graphing and printing results to screen
if graph==1
nam='Variation in Exponent and Coeff';
figure('name',nam)
subplot(2,1,1)
hold on;
plot(OP.b,'r.-')
plot(PLb*ones(Niter,1),'k-')
ylabel('Exponent')
grid on;
subplot(2,1,2)
hold on;
plot(10.^OP.loga*1000,'g.-')
plot(PLa*ones(Niter,1),'k-')
ylabel('Coeff (ms)')
grid on;
legend('Monte Carlo values','Actual scaling')


nam='Histogram';
figure('name',nam)
subplot(2,1,1)
hold on;
h1=histogram(OP.b);
plot([PLb PLb],[0 max(h1.BinCounts)+20],'k-','LineWidth',4)
grid on;
xlabel('Exponent')
ylabel('n')
subplot(2,1,2)
hold on;
h1=histogram(1000*10.^OP.loga);
plot([PLa PLa],[0 max(h1.BinCounts)+20],'b-','LineWidth',4)

grid on;
xlabel('Coefficient (ms)')
ylabel('n')


disp(['------------------Final Output-------------'])
disp(['Scaling of InvPend inertial delay using mean values::' num2str(PLa) '*M^' num2str(PLb) 'ms'  ]);
disp(['Coefficeint of InvPend inertial delay using Monte Carlo::' num2str(OP.amean*1000) ' +-95%CI [' num2str(OP.aCI*1000) ']'  'ms']);
disp(['Exponent of InvPend inertial delay using Monte Carlo::' num2str(OP.bmean) ' +-95%CI [' num2str(OP.bCI) ']'  ]);  

end