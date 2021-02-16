%% Simulation master code
% This is the run file for a simulation without any control moments
% Note: 
%% setup simulation variables
clear mex;clear;clc;close all; %clear/close memories
% mex OSCAR_MEX.c % Create MEX if not already done
Tfinal_master = 1; %total sim time, seconds
refresh = .02;% How often to recalculate torques (seconds), must be small
Ctorque1 = 0/1000;%.001;
Ctorque2 = 0/1000;%.002;
Ctorque3 = 2/1000;%003;
Tmax = 2/1000;

%% Run simulation
outputs = sim('Running_mex_SIMULINK_No_controller');
%% Plot outputs

%Convert Euler Angles
EP1 = outputs.fout(:,11);
EP2 = outputs.fout(:,12);
EP3 = outputs.fout(:,13);
EP4 = outputs.fout(:,14);
T = linspace(0,Tfinal_master,(Tfinal_master/refresh)+1);
EPmat = [EP1,EP2,EP3,EP4];
i = 1;
q = [];
while i <= length(EPmat)
    q(i,:) = quat2eul(EPmat(i,:));
    i = i+1;
end
q1 = q(:,1);
q2 = q(:,2);
q3 = q(:,3);

%plot Euler Angles
figure(1)
hold on
plot(T,q1,'r')
plot(T,q2,'*')
plot(T,q3,'k')
legend('q1','q2','q3')
title('Euler Angles')
hold off



%Plot Angular Velocities
w1 = outputs.fout(:,21);
w2 = outputs.fout(:,22);
w3 = outputs.fout(:,23);
figure(2)
hold on
plot(T,w1,'r')
plot(T,w2,'*')
plot(T,w3,'k')
legend('w1','w2','w3')
title('Angular Velocities')
hold off




