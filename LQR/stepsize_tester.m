clear;clc;close all;

stepsize = 1;
Tfinal = 20;
CM = [];
CMavgs = [];
timearr = [];
steps = [];
CMmin = 100;
while stepsize >=.0001
    tic;
    CMscurr = running_simulink_function_stepsize(stepsize,Tfinal);
    Timeel = toc
    if mean(CMscurr) <= CMmin
        CMmin = mean(CMscurr);
        CMmin_stepsize = stepsize;
    end
    CM = [CM;CMscurr];
    CMavgs = [CMavgs;mean(CMscurr)];
    timearr = [timearr;Timeel];
    steps = [steps;stepsize];
    stepsize = stepsize/2
    close all;
end
fprintf("The minimum moment was:")
disp(CMmin)
fprintf('\n"It occured at a stepsize of ')
disp(CMmin_stepsize)
set(gcf,'color','w');
figure(1)
hold on
plot(steps,CMavgs.*3,'-x')
title('Total Control Moment Required')
ylabel('Total Moment Required (N-m-s)')
xlabel('Step Size (Seconds)')
hold off
figure(2)
plot(steps,timearr,'-x')
title('Simulation Time Required')
ylabel('Time (Seconds)')
xlabel('Step Size (Seconds)')



