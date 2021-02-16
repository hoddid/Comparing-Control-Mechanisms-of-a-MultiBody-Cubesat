function CMs = running_simulink_function_stepsize(refresh,Tfinal_master)
%% Running simulink file
%clear mex;
%mex OSCAR_MEX.c
%clear;clc;
%close all;
tic
Tfinal_master = Tfinal_master; %total sim time, seconds
%refresh = .1;% How often to recalculate torques (seconds), must be small

%T_real_torun = ((Tfinal_master*refresh)*1000)/(refresh/.001)*(1/8)*2;
%fprintf("This code will take approximatly %2f Minute(s) to run\n",T_real_torun)
%dt = datetime( 'now', 'InputFormat', 'HH:mm:ss' );
%dt.Format = 'HH:mm:ss';
%dt = char(dt + seconds(T_real_torun*60));
%fprintf("This means the Code will be done around %s \n",dt)
%fprintf("This will be longer than expected if velocities are high,\nand shorter if velocities are low\n\n")


%Good for 1.5 vel offset
% C1 = [10,5,10]; % C1 = [4,1,1]; was pretty good, with c2 .01 and R = 3
% C2 = [50,10,50];%.00000001;%.00001; %150
% R = 500;

% % Good for a 1 radian, -1 radian .1 radia offset
% C1 = [70,15,15]; % C1 = [4,1,1]; was pretty good, with c2 .01 and R = 3
% C2 = [.01,7,15];%.00000001;%.00001; %150
% R = 500;

C1 = [70,15,15]; % C1 = [4,1,1]; was pretty good, with c2 .01 and R = 3
C2 = [.01,7,15];%.00000001;%.00001; %150
R = 500;





%initial torques
Ctorque1 = 0/1000;%.001;
Ctorque2 = 0/1000;%.002;
Ctorque3 = 0/1000;%003;
Tmax = 2/1000;

%desired angles
%angdes = [.61 -.65 0.1];
angdes = [0 0 0];

angdes = fliplr(angdes);
%angdes = [pi/12 pi/6 pi/12];
EPdes = eul2quat(angdes);

%Desired euler parameters
EP1des = EPdes(1);
EP2des = EPdes(2);
EP3des = EPdes(3);
EP4des = EPdes(4);


U4des = 0;
U5des = 0;
U6des = 0;
%check to ensure refresh rate and tfinal_master line up
%if they dont, make them
div_even = mod(Tfinal_master,refresh);
if div_even ~= 0
    disp('Your Tfinal cannot be evenly divided by the refresh rate specified')
    NumSteps = round(Tfinal_master/refresh);
    Tfinal_master = refresh*NumSteps;
    fprintf('The Tfinal has been changed to %f in order to accomadate the step size\n',Tfinal_master)
end
options = simset('SrcWorkspace','current');
outputs = sim('Running_mex_SIMULINK',[],options);
%% PLot stuff
EP1 = outputs.fout(:,11);
EP2 = outputs.fout(:,12);
EP3 = outputs.fout(:,13);
EP4 = outputs.fout(:,14);
T = linspace(0,Tfinal_master,(Tfinal_master/refresh)+1);
%T = outputs.Time.';
EPmat = [EP1,EP2,EP3,EP4];
i = 1;
% q = [];
% while i <= length(EPmat)
%     q(i,:) = quat2eul(EPmat(i,:));
%     i = i+1;
% end
q1 = outputs.Angs(1);
q2 = outputs.Angs(2);
q3 = outputs.Angs(3);

figure(1)
hold on
plot(T,q1,'r')
plot(T,q2,'*')
plot(T,q3,'k')
yline(angdes(1),'b-')
legend('Angle 1','Angle 2','Angle 3','Desired Angle')
title('Euler Angles LQR')
xlabel('Time (seconds)')
ylabel('Radians')
hold off

%Save file name
baseFileName = sprintf('STEPSIZE_%.3g_Angles_C11_%.3g_C12_%.3g_C13_%.3g_C21_%.3g_C22_%.3g_C23_%.3g_R_%.3g.fig',refresh,C1(1),C1(2),C1(3),C2(1),C2(2),C2(3),R);
pathName = fullfile('\automated', baseFileName);  
saveas(figure(1),[pwd pathName]);



w1 = outputs.fout(:,21);
w2 = outputs.fout(:,22);
w3 = -1.*outputs.fout(:,23); %because axis is 'upside down', flip so they are all consistent

w1 = -1.*outputs.Vel(:,1);
w2 = outputs.Vel(:,2);
w3 = -1.*outputs.Vel(:,3);


figure(2)
hold on
plot(T,(w1),'r')
plot(T,(w3),'*')
plot(T,(w2),'k')
legend('\omega 1','\omega 2','\omega 3')
title('Angular Velocities LQR')
xlabel('Time (seconds)')
ylabel('Radians per Second')
hold off

%Save file name
baseFileName = sprintf('STEPSIZE_%.3g_Velocities_C11_%.3g_C12_%.3g_C13_%.3g_C21_%.3g_C22_%.3g_C23_%.3g_R_%.3g.fig',refresh,C1(1),C1(2),C1(3),C2(1),C2(2),C2(3),R);
pathName = fullfile('\automated', baseFileName);  
saveas(figure(2),[pwd pathName]);

figure(3)
hold on
title("Moments")
plot(T,smooth(outputs.control_moments(:,1),50),'r')
plot(T,smooth(outputs.control_moments(:,2),50),'*')
plot(T,smooth(outputs.control_moments(:,3),50),'k')
legend("Moment 1","Moment 2","Moment 3")
ylabel("N-M")
xlabel('Time (seconds)')

%Save file name
baseFileName = sprintf('STEPSIZE_%.3g_Moments_C11_%.3g_C12_%.3g_C13_%.3g_C21_%.3g_C22_%.3g_C23_%.3g_R_%.3g.fig',refresh,C1(1),C1(2),C1(3),C2(1),C2(2),C2(3),R);
% Specify some particular, specific folder:
%pathName = fullfile('C:\Users\hoddid\Documents\School_RPI\Thesis\Mex_correct_al\controllers\LQR\automated', baseFileName);  
pathName = fullfile('\automated', baseFileName);  

%figure(3); % Activate the figure again.
%NameAndPath = sprintf('%s',pathName);
%saveas(gcf,NameAndPath);
saveas(figure(3),[pwd pathName]);


hold off

CM1 = (smooth(outputs.control_moments(:,1),100));
CM2 = (smooth(outputs.control_moments(:,2),100));
CM3 = (smooth(outputs.control_moments(:,3),100));

CMT1 = trapz(CM1,T);
CMT2 = trapz(CM2,T);
CMT3 = trapz(CM3,T);
CMs = [CMT1,CMT2,CMT3];


tel = toc;
fprintf("The elapsed time was %f minutes\n",tel/60)
close all;
%% alarm
%play a sound when completed to notify user
%as code often takes 1-3 minutes to run
% load chirp.mat;
% chirp = audioplayer(y, Fs);
% play(chirp)
end