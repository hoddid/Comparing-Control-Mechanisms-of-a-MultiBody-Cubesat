%% Running simulation main function
% This is the main function for runing the simulations
% Fill in the appropriate simulation variables under function setup to
% customize the simulation length, refresh rate of the control system,
% desired angles and gains.
%% function setup
clear mex; %clear the mex memory
%mex OSCAR_MEX.c %create MEX if it isnt already made
clear;clc;
close all;

tic %begin timer
global Fname
Fname = "input.xlsx";

Tfinal_master = 25; %total sim time, seconds
refresh = .001;% How often to recalculate torques (seconds), must be small

%Calculate approximate system runtime (trial and error, computer dependent)
T_real_torun = ((Tfinal_master*refresh))/(refresh)*(1/2);
fprintf("This code will take approximatly %2f Minute(s) to run\n",T_real_torun)
dt = datetime( 'now', 'InputFormat', 'HH:mm:ss' );
dt.Format = 'HH:mm:ss';
dt = char(dt + seconds(T_real_torun*60));
fprintf("This means the Code will be done around %s \n",dt)
fprintf("This will be longer than expected if velocities are high,\nand shorter if velocities are low\n\n")

%Choose Gains
% C1 affects angles, with the first value of C1 affecting the first axis,
% the second the second, etc.
% C2 affects the angular velocities, with the first value affecting the first axis,
% the second the second, etc.
% R affects the control effort. A higher R value will make the system use
% less control. Can be useful for eliminating oscilation.



C1 = [1,1,1]; 
C2 = [150,110,150];
R = 1000000;







%initial torques
Ctorque1 = 0/1000;%.001;
Ctorque2 = 0/1000;%.002;
Ctorque3 = 0/1000;%003;
Tmax = 2/1000;

% Set desired euler Angles

angdes = [0 0 0];

angdes = fliplr(angdes); %flip and make into euler parameters to align
%with how the multibody model is written
EPdes = eul2quat(angdes);

%Desired euler parameters
EP1des = EPdes(1);
EP2des = EPdes(2);
EP3des = EPdes(3);
EP4des = EPdes(4);

%Set Angular Velocity desired rates
U4des = 0;
U5des = 0;
U6des = 0;

%check to ensure refresh rate and tfinal_master line up
%if they dont, slightly adjust the final time to make them fit neatly
%togehter
div_even = mod(Tfinal_master,refresh);
if div_even ~= 0
    disp('Your Tfinal cannot be evenly divided by the refresh rate specified')
    NumSteps = round(Tfinal_master/refresh);
    Tfinal_master = refresh*NumSteps;
    fprintf('The Tfinal has been changed to %f in order to accomadate the step size\n',Tfinal_master)
end

%Run simulink Simulation
outputs = sim('Running_mex_SIMULINK');
%% Plot simulation outputs
%set the figure background colour to be white
set(gcf,'color','w');

T = linspace(0,Tfinal_master,(Tfinal_master/refresh)+1);


%Fetch output variables from simulation outputs
q1 = outputs.Angs(:,1);
q2 = outputs.Angs(:,2);
q3 = outputs.Angs(:,3);

% Plot Euler Angles
figure(1)
hold on
plot(T,q1,'r')
plot(T,q2,'b')
plot(T,q3,'k')
yline(angdes(1),'r-')
yline(angdes(2),'b-')
yline(angdes(3),'k-')
legend('Angle 1','Angle 2','Angle 3','Desired Angle 1','Desired Angle 2','Desired Angle 3')
title('Euler Angles LQR')
xlabel('Time (seconds)')
ylabel('Radians')
hold off






%Get the angular velocities to plot
w1 = outputs.Vel(:,1);
w2 = outputs.Vel(:,2);
w3 = outputs.Vel(:,3);

set(gcf,'color','w');
%Plot Angular Velocities.
figure(2)
hold on
plot(T,(w1),'r') 
plot(T,(w3),'b')
plot(T,(w2),'k')
legend('\omega 1','\omega 2','\omega 3')
title('Angular Velocities LQR')
xlabel('Time (seconds)')
ylabel('Radians per Second')
hold off

%plot Control Moments
figure(3)
hold on
title("Moments")
plot(T,outputs.control_moments(:,1),'r')
plot(T,outputs.control_moments(:,2),'b')
plot(T,outputs.control_moments(:,3),'k')
legend("Moment 1","Moment 2","Moment 3")
ylabel("N-M")
xlabel('Time (seconds)')

hold off

%Calculate total moment output from reaction wheels
CM1 = (outputs.control_moments(:,1));
CM2 = (outputs.control_moments(:,2));
CM3 = (outputs.control_moments(:,3));

CMT1 = trapz(CM1,T);
CMT2 = trapz(CM2,T);
CMT3 = trapz(CM3,T);
fprintf("The Total Control moment output in the 1 Direction was %3f N-M-S\n",CMT1)
fprintf("The Total Control moment output in the 2 Direction was %3f N-M-S\n",CMT2)
fprintf("The Total Control moment output in the 3 Direction was %3f N-M-S\n",CMT3)





tel = toc;
fprintf("The elapsed time was %f minutes\n",tel/60)
%% alarm
%play a sound when completed to notify user
%as code often takes a signifigant amount of time to run
load chirp.mat;
chirp = audioplayer(y, Fs);
play(chirp)
