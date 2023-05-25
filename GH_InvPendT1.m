function [t,y,Tswitch]=GH_InvPendT1(B,L,M,yi)
% Function to find the optimum time to switch direction of bang bang torque
% simulate the results and output inertial delay
global g
graph=0;
%% Test case : Comment out function definition above and uncomment below section. Aso comment out the "end" at line 165
%{
clear all;close all;clc;
graph=1;
g=-9.8066;
% Scaling values based on Mass
%M=0.005;XX=1000;Tunit='(ms)';% for 5 gram
%M=1;XX=1000;Tunit='(ms)';% for 1 kg
%M=5000;XX=1;Tunit='(s)';% for 5 gram
M=10000;XX=1;Tunit='(s)';

% Limb inertial properties from fitlm Type 1 regression, all final values in SI units
L_HlimbA=0.163;L_HlimbB=0.357;% Mammalian fore limb
L_Hlimb=L_HlimbA*M.^L_HlimbB;
L_FlimbA=0.161;L_FlimbB=0.384;% Mammalian fore limb
L_Flimb=L_FlimbA*M.^L_FlimbB;
%Llimb=L_Hlimb;% old method
L=(L_Hlimb+L_Flimb)./2;
% Muscle scaling. Alexander et. al.  (1981), Table 1 All non-hoppers, ankle extensors
% Assuming muscle force (100 N)=muscle mass (kg)/ fiber length (m)
Mmusc=(5.1*M.^0.97)/1000;% mass of muscle (kg), using ankle extensors
FLmusc=(10.6*M.^0.14)/1000;% fiber length muscle (m)
rho=1060;%density of muscle in kg/m^3 (Mendez and keys 1960)
CSmusc=Mmusc./rho./FLmusc.*(100^2);% Cross sectional area of muscle in cm^2
Fmusc=CSmusc.*20;% assuming 20N/cm^2 muscle isometric stress. varies from 20-30 N/cm^2
% REfs: Allometry of muscle,tendon and elastic energy storage capacity in
% mammals Shadwick and Pollock 1994
% Close R 1972 Dynamic Properties of mammalian skeletal muscles 
Amusc=(9.4*M.^0.38)/1000;% Moment arm value from Table II ankle extensors-all non hoppers
B=4*Fmusc.*Amusc;% Torque by the 4 leg muscles (Nm)

%mode=0;Iangle=-5;modnam='Angle Perturbation';% set mode to 0 for initial position/0 velocity  
mode=1;Fr=-0.21;modnam='Froude Perturbation';%set mode to 1 for initial velocity/0 position
%mode=1;Fr=-0.49;modnam='Froude Perturbation';%set mode to 1 for initial velocity/0 position

%Fr=-sqrt(2)/3;
%Fr=-0.5199359;

if mode==0 % sim recovery from angle perturbation
    yi(1) = deg2rad(Iangle);
    yi(2) = 0;
elseif mode==1 % sim recovery from angular velocity perturbation
    yi(1) = 0;
    Vpert=Fr.*sqrt(L.*-g);
    Wpert=Vpert./L; 
    yi(2)=Wpert;
end
%}
%% Find optimal switch time
Tswitchi = 0.1; % initial guess at swtich time
options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
Tswitch = fsolve(@(Tswitch)pendObjective(Tswitch,yi,B,L,M), Tswitchi, options);
% simulate using optimal switch time
[t,y] = pendSim(Tswitch,yi,B,L,M);
% creating torque vector
indTswitch=find(abs(t-Tswitch)<1e-5);
Torque_Musc=ones(length(t),1);
Torque_Musc(1:indTswitch)=B;
Torque_Musc(indTswitch:end)=-B;
Torque_Grav=M*-g*L.*sin(y(:,1));

%disp(['Mass:' num2str(M) ' Fall Time:' num2str(t(end))]); 
% plot results
if graph==1
    figure(1); clf; hold on
    plot(t,radtodeg(y));
    plot([Tswitch Tswitch],[-50 50],'k.-')
    grid on
    xlabel('time (s)');
    ylabel('Angle (deg) & Angular velocity (deg/sec)')
    disp('');
    legend('Angle (deg)','Angular velocity (deg/sec)','Tswitch')
    

    nam=['Inverted pendulum-' modnam];
    figure('name',nam)
    subplot(3,1,1)
    plot(t,radtodeg(y(:,1)),'r-','LineWidth',3)
    ylabel('Angle (deg)')
    grid on;
%title('5 gram shrew')
%title('5 ton elephant')

    subplot(3,1,2)
    plot(t,radtodeg(y(:,2)),'b-','LineWidth',3)
    ylabel('Angular Velocity (deg/sec)') 
    grid on;

    subplot(3,1,3)
    hold on;
    plot(t,Torque_Musc,'k-','LineWidth',3)
    plot(t,Torque_Grav,'g-','LineWidth',3)
    legend('Muscle','Gravity')
    grid on;
    ylabel('Torque (Nm)')
    xlabel('Time (s)')


end

% % ID paper graph Posture task kinetics V4.3
% nam=['Inverted pendulum-' modnam '-' num2str(M)];
% figure('name',nam)
% %title('5 gram shrew')
% title('5 ton elephant')
% x0=100;
% y0=100;
% width=240;
% height=360;
% set(gcf,'units','points','position',[x0,y0,width,height])   
%     
%     subplot(3,1,1)
%     plot(t*XX,radtodeg(y(:,1)),'r-','LineWidth',3)
%     ylabel('Angle (deg)')
%     axis([min(t)*XX max(t)*XX min(radtodeg(y(:,1)))+0.1*min(radtodeg(y(:,1))) max(radtodeg(y(:,1)))+0.1*max(radtodeg(y(:,1))) ])
%     
%     subplot(3,1,2)
%     plot(t*XX,radtodeg(y(:,2)),'r-','LineWidth',3)
%     ylabel('Angular Velocity (deg/sec)') 
%     axis([min(t)*XX max(t)*XX min(radtodeg(y(:,2)))+0.1*min(radtodeg(y(:,2))) max(radtodeg(y(:,2)))+0.1*max(radtodeg(y(:,2))) ])
% 
%     
%     subplot(3,1,3)
%     hold on;
%     plot(t*XX,Torque_Musc,'r-','LineWidth',3)
%     ylabel('Torque (Nm)')
%     xlabel(['Time ' Tunit])
%     axis([min(t)*XX max(t)*XX min(Torque_Musc)+0.1*min(Torque_Musc) max(Torque_Musc)+0.1*max(Torque_Musc) ])

if graph==1
% ID paper graph Posture task kinetics V5
nam=['Inverted pendulum-' modnam '-' num2str(M)];
figure('name',nam)

x0=10;
y0=10;
width=6;
height=8;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    
    subplot(3,1,1)
    hold on;
    plot(t*XX,Torque_Musc,'r-','LineWidth',3)
    %ylabel('Torque (Nm)')
    axis([min(t)*XX max(t)*XX min(Torque_Musc)+0.1*min(Torque_Musc) max(Torque_Musc)+0.1*max(Torque_Musc) ])

    
    subplot(3,1,2)
    plot(t*XX,radtodeg(y(:,2)),'r-','LineWidth',3)
    %ylabel('Angular Velocity (deg/sec)') 
    axis([min(t)*XX max(t)*XX min(radtodeg(y(:,2)))+0.1*min(radtodeg(y(:,2))) max(radtodeg(y(:,2)))+0.1*max(radtodeg(y(:,2))) ])  
       
    subplot(3,1,3)
    plot(t*XX,radtodeg(y(:,1)),'r-','LineWidth',3)
   % ylabel('Angle (deg)')
    xlabel(['Time ' Tunit])
    axis([min(t)*XX max(t)*XX min(radtodeg(y(:,1)))+0.1*min(radtodeg(y(:,1))) max(radtodeg(y(:,1)))+0.1*max(radtodeg(y(:,1))) ])
end

%%
% if running as script, comment out end below
end


%% In line function definitions

function f = pendObjective(Tswitch,yi,B,L,M)
% this is the objective function that we are trying to minimize.
[t,y] = pendSim(Tswitch,yi,B,L,M);

% minimize this objective function. Note that if the event crossing
% occurred, the second term in this objective function will be nearly zero.
%f = y(end,1).^2 + y(end,2).^2;
f = (y(end,1)*1000)^2 + (y(end,2))^2;
end

function dy = pendEOM(t,y,B,L,M)
% pendulum equations of motion.
global g
dy(1,1) = y(2);
dy(2,1) = (B+M*9.8066*L*sin(y(1)))/(M*L^2);
end

function [position,isterminal,direction] = pendEvents(t,y)
% the event function. We could halt the sim when either angle or velocity
% is zero, but the optimization works much better when we set the event to
% zero velocity and the objective function to zero angle.
position = y(2); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end

function [t,y] = pendSim(Tswitch,yi,B,L,M)
% Once we have an optimal Tswitch, redo the simulation one last time to
% inspect the results. 

% first simulate from t=0 to t=Tswitch. 
[t,y] = ode45(@(t,y) pendEOM(t,y,B,L,M), [0 Tswitch], yi);

% use the end state and time as the initial state and time for the new sim
yinew = y(end,:);
tspannew = [t(end) 20];
B = -B; % switch the torque direction.

% this sim stops the sim when the angular velocity is zero, and not at a
% specific time.
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'Events',@pendEvents);
[tnew,ynew,te,ye] = ode45(@(t,y) pendEOM(t,y,B,L,M), tspannew, yinew,options);
% one interesting thing to note is that if I set the event to be zero
% angle, the optimization would not converge on zero angular velocity. But
% when I set the even to zero velocity, and searched for the Tswitch that
% yielded zero angle, the optimization converged. 

% final time and states are the two sim results joined together.
t = [t(1:end-1);tnew];
y = [y(1:end-1,:);ynew];

end

