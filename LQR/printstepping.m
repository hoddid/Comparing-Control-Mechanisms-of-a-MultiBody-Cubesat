function y = printstepping(u)
% This function Prints the current time to the MATLAB console, so the user
% can tell how far along the simulation is. Because the timestep is so
% small, not every timestep is printed, only some of them. This prevents
% clogging up the console. 

if u > 1 %only print if above 1 second
     div10 = mod(round(u,2),2); 
      if div10 == 0 %rounded timestep divisible by 2?
         fprintf('Tsim = %f\n',round(u,2))  %print it out
      end
 end
y = u;
end
