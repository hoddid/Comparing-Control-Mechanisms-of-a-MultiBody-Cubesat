function K = linearize_get_Gain(inputs)
%This Code takes in the current states, linearizes about them, 
% and returns the optimal gain matrix, K.
%global C1 C2 R

%create states:
Ctorques_in = [inputs(1),inputs(2),inputs(3)].'; % assume no moments
input_states = [EP2ang_arr(inputs(14:17)),inputs(24:26)].'; % Do linearization about the uler ANGLEs not params.
C1 = inputs(end-6:end-4);
C2 = inputs(end-3:end-1);
R = inputs(end);

%get your A and B matricies (linearized matricies which form an equation xdot = Ax+Bu)
[A,B] = GetLinModFtxu(@running_mex_File_getlinmod,input_states,Ctorques_in,inputs);

%^the function getlinmod gets is supposed to take in the file which makes, 
% the state inputs (position, velocity, etc) to the getmex (Y0) and the
% inputs (moments)

%C = eye(6);
%D = zeros(6,3);
%R = 50;
Q1 = 1;
%C1 = 10;
%C2 = .5;
Q = [C1(1)*Q1 0 0 0 0 0;
    0 C1(2)*Q1 0 0 0 0;
    0 0 C1(3)*Q1 0 0 0;
    0 0 0 C2(1)*Q1 0 0;
    0 0 0 0 C2(2)*Q1 0;
    0 0 0 0 0 C2(3)*Q1];
% Note to future dave: These correspond to the design variables?
% mayve add a constant to some, so it cares more about the position , as
% that doesnt converge fast.

%R = eye(3);
%sym = issymmetric(R) % symm?
%[~,p] = chol(R);% pos def?

coder.extrinsic('lqr'); 
%K = zeros(3,6);
K = lqr(A,B,Q,R);
end