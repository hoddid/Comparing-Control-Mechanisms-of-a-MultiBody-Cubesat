function K = linearize_get_Gain(inputs)
% This Code takes in the current states, linearizes about them, 
% and returns the optimal gain matrix, K.

%the try catch block prints out if an error occurs, and that it happened
%during linearization. This is useful, as MATLAB's builtin error messages
%are often not very helpful about where the error occurs
try
%Separate states required for linearization from all of the model outputs
Ctorques_in = [inputs(1),inputs(2),inputs(3)].'; %Moments
input_states = [quat2eul(inputs(14:17)),inputs(24:26)].'; %Euler Agngles and anguler velocities

%Fetch Gains (chosen by user)
C1 = inputs(end-6:end-4);
C2 = inputs(end-3:end-1);
R = inputs(end);

%get your A and B matricies (linearized matricies which form an equation xdot = Ax+Bu)
[A,B] = GetLinModFtxu(@running_mex_File_getlinmod,input_states,Ctorques_in,inputs);


% Create Q matrix for cost function generation
Q = [C1(1) 0 0 0 0 0;
    0 C1(2) 0 0 0 0;
    0 0 C1(3) 0 0 0;
    0 0 0 C2(1) 0 0;
    0 0 0 0 C2(2) 0;
    0 0 0 0 0 C2(3)];


%Set the LQR function to be coded as extrinsic (required for SIMULINK to function)
coder.extrinsic('lqr'); 

%create Gain matrix via LQR
K = lqr(A,B,Q,R);
catch e %print out the error message, if one occurs.
    fprintf(1,'There was an error during Linearization! The message was:\n%s',e.message);
    K = zeros(3,6);
    
end