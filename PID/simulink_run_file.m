%% Running simulation main function
% This is the main function for runing the simulations
% Fill in the appropriate simulation variables under function setup to
% customize the simulation length, refresh rate of the control system,
% desired angles and gains.
%% function setup
clear mex; % Clear the MEX memory
%mex OSCAR_MEX.c%create MEX if it isnt already made
clear;clc;
close all;


tic %begin timer
global Fname
Fname = "input.xlsx";
%Consider using the full file path name here, i.e. C:/users/....
Tfinal_master = 30; %total sim time, seconds
refresh = .005;% How often to recalculate torques (seconds), must be small

%Calculate approximate system runtime (trial and error, computer dependent)
T_real_torun = ((Tfinal_master*refresh)*1000)/(refresh/.001)*(3/8);
fprintf("This code will take approximatly %2f Minute(s) to run\n",T_real_torun)
dt = datetime( 'now', 'InputFormat', 'HH:mm:ss' );
dt.Format = 'HH:mm:ss';
dt = char(dt + seconds(T_real_torun*60));
fprintf("This means the Code will be done around %s \n",dt)
fprintf("This will be longer than expected if velocities are high,\nand shorter if velocities are low\n\n")



%For Angular position errors
Kp = .065; 
Kd = .105;
Kint = 0.01;
Kint_time = .5;

% %for angular velocities
% Kp = .065; 
% Kd = .505;
% Kint = 0.01;
% Kint_time = .5;




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

%initial torques
Ctorque1 = 0/1000;%.001;
Ctorque2 = 0/1000;%.002;
Ctorque3 = 0/1000;%003;
Tmax = 2/1000;

%desired angles
angdes = [0 0 0];

%angdes = [pi/12 pi/6 pi/12];
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
%% Run Sim
outputs = sim('Running_mex_SIMULINK_PID');
%% PLot stuff
set(gcf,'color','w');
%Fetch Time
T = outputs.Time;
%Fetch output variables from simulation outputs
q1 = outputs.angles(:,1);
q2 = outputs.angles(:,2);
q3 = outputs.angles(:,3);

%plot Euler Angles
figure(1)
hold on
plot(T,q1,'r')
plot(T,q2,'b')
plot(T,q3,'k')
yline(angdes(1),'b-')
% yline(angdes(2),'r-')
% yline(angdes(3),'g-')
%legend('q1','q2','q3','ang1 des','ang2des','ang3 des')
legend('Angle 1','Angle 2','Angle 3','Desired Angle')

title('Euler Angles')
xlabel('Time (seconds)')
ylabel('Radians')
hold off

%Get the angular velocities to plot
%some are multiplied by negative 1 to correct for how the model defines
%direction, to ensure that the user can view everything consistantly.
w1 = -1.*outputs.Vels(:,1);
w2 = outputs.Vels(:,2);
w3 = outputs.Vels(:,3);
figure(2)
hold on
plot(T,w1,'r')
plot(T,w2,'b')
plot(T,w3,'k')
legend('\omega 1','\omega 2','\omega 3')
title('Angular Velocities')
xlabel('Time (seconds)')
ylabel('Radians per Second')
hold off

%plot moments
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
CM1 = abs(outputs.control_moments(:,1));
CM2 = abs(outputs.control_moments(:,2));
CM3 = abs(outputs.control_moments(:,3));

CMT1 = trapz(CM1,T);
CMT2 = trapz(CM2,T);
CMT3 = trapz(CM3,T);
fprintf("The Total Control moment output in the 1 Direction was %3f N-M-S\n",CMT1)
fprintf("The Total Control moment output in the 2 Direction was %3f N-M-S\n",CMT2)
fprintf("The Total Control moment output in the 3 Direction was %3f N-M-S\n",CMT3)


% Other potentially useful things to plot
% Ang_err = outputs.Angle_err;
% ang_err1 = Ang_err(:,1);
% ang_err2 = Ang_err(:,2);
% ang_err3 = Ang_err(:,3);

% figure(4)
% hold on
% title("Angle Errors")
% plot(T,ang_err1,'r')
% plot(T,ang_err2,'*')
% plot(T,ang_err3,'k')
% legend("Err1","Err2","Err3")
% ylabel("Radians")
% 
% hold off
% Prop_moment = out.Prop_moment;
% figure(5)
% hold on
% title("Prop gain moment")
% plot(T,Prop_moment(:,1),'*')
% plot(T,Prop_moment(:,2),'r')
% plot(T,Prop_moment(:,3),'k')
% legend('axis 1','axis 2','axis 3')
% %ylim([-2/1000 2/1000])
% hold off

% Deriv_moment = out.Deriv_moment;
% figure(6)
% hold on
% title("Deriv gain moment")
% plot(T,Deriv_moment(:,1),'b')
% plot(T,Deriv_moment(:,2),'r')
% plot(T,Deriv_moment(:,3),'k')
% legend('axis 1','axis 2','axis 3')
% %ylim([-2/1000 2/1000])
% hold off



% Int_moment = out.Int_moment;
% figure(7)
% hold on
% title("Int gain moment")
% plot(T,Int_moment(:,1),'b')
% plot(T,Int_moment(:,2),'r')
% plot(T,Int_moment(:,3),'k')
% legend('axis 1','axis 2','axis 3')
% ylim([-2/1000 2/1000])

%hold off

tel = toc;
fprintf("The elapsed time was %f minutes\n",tel/60)
%% alarm
%play a sound when completed to notify user
%as code often takes a signifigant amount of time to run
load chirp.mat;
chirp = audioplayer(y, Fs);
play(chirp)
