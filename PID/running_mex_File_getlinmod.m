function c = running_mex_File_getlinmod(input_states,Ctorques,other_ins)
%this function exists for getting the linearized model 
%getlinmodfxtu requires the model to be of the format model(states,inputs)
%to linearize it for getting the A and B, matrix, thus this is included. 

%mex helloworld.c
%Assemble input vector using disjointed inputs
%needed because GETLINMODfxtu increases all inputs given,
%and increasing things like inertia to linearize would be inaccurate,
%as only motion states are changing.
Inertias = other_ins(4:12);%,other_ins(5),other_ins(6),other_ins(7)...
EPs = eul2quat(input_states(1:3).');
Vels = input_states(4:6).';
Qs = other_ins(18:20);
Massb = other_ins(13);
U1_3 = other_ins(21:23);
%set Tfinal to be small (slightly larger size than integration step), else linearization will be about huge margin
U7_9_and_times_errs = [other_ins(27:30),other_ins(32)*1.5,other_ins(32:35)];

input = [Ctorques.',Inertias,Massb,EPs,Qs,U1_3,Vels,U7_9_and_times_errs];
results = OSCAR_MEX(input);
clear mex;
%get rough numerical derivative of Qs
Qs_init = input_states(1:3);
%Q1 =  EP2ang(EPS_init);
Qs_result = EP2ang([results(11);results(12);results(13);results(14)]);
%Q2 =  EP2ang(EPS_result);
Qdot = (Qs_result.'-Qs_init)/results(end-3);

%get rough numerical derivative of Rotational velocity
Vels_init = input_states(4:6);
%Q1 =  EP2ang(EPS_init);
Vels_result = [results(21);results(22);results(23)];
%Q2 =  EP2ang(EPS_result);
Vels_dot = (Vels_result-Vels_init)./results(end-3);

Xdot = [Qdot;Vels_dot];
c = Xdot;
%c = [results(10),results(11),results(12),results(13),results(20),results(21),results(22)];
%c = [Q,results(20),results(21),results(22)];
end