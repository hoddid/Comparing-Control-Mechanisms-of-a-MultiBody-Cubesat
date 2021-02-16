function [all_init] = get_initial_inputs()
% This function reads the initial values from the excel sheet via a
% automatically generated MATLAB excel input read file, which was modified
% such that the the file name is taken as an input.
global Fname
[CTORQUE1,CTORQUE2,CTORQUE3,IB11, IB12, IB22, IB23,IB31, IB33,IRW1,IRW2, IRW3, MASSB, EP1, EP2, EP3, EP4, Q1, Q2, Q3, U1, U2, U3, U4, U5, U6, U7, U8, U9, TINITIAL, TFINAL, INTEGSTP, PRINTINT, ABSERR, RELERR] = getInputs(Fname);
all_init = [CTORQUE1, CTORQUE2, CTORQUE3, IB11, IB12, IB22, IB23,IB31, IB33, IRW1, IRW2, IRW3, MASSB, EP1, EP2, EP3, EP4, Q1, Q2, Q3, U1, U2, U3, U4, U5, U6, U7, U8, U9, TINITIAL, TFINAL, INTEGSTP, PRINTINT, ABSERR, RELERR];
end