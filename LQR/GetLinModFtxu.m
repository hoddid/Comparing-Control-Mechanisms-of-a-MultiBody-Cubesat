function  [A,B]=GetLinModFtxu(ftxu,xs,us,others) 
% Hurst, Joshua, 2019. State Space Linearization for Mechatronics 4520
%% [A B]=GetLinModFtxu(ftxu,t,xs,us) 
% ftxu if the function that evaluates xdot
% is must accept inputs t,x,u 
%
% Numerical determination of state space matrices:
% If you have f(x,u) for xdot=f(x,u) for your function then you 
% can numerically determine the values of the partial derivatives.
% To simulate the nonlinear equations you have to put the equations
% in this form (typically) so there is little 
% work involved to numerically verify your expressions 
% for the partial derivatives to get A and B:

% The A matrix is defined as the partial derivative of f with respect to x:
% A = pf/px.  Numerically we can approximate this derivative:
% A ~ [f(x+dx,u)-f(x-dx,u)]/2d
% Similarly
% B ~ [f(x,u+du)-f(x,u-du)]/2d
%
% time argument is typically not used, but this argument is passed in so 
% you can determine A, B with the same function used with the ODE
% solver instead of writing a new one
 
%% numerically determine A:
n=length(xs); % the number of state variables
% a small change in x:
d=1e-7;  % "small" with respect to typically values of x or u.

% determine [pf/px1 pf/px2 - - - pf/pxn] by looping through
% each state variable x1, x2 - - -, xn
A = zeros(6);
for i=1:n % loop through a small change in each state:
    dx=zeros(n,1);
    dx(i,1)=d;
    %disp('A matrix')
    A(:,i)=(ftxu(xs+dx,us,others)-ftxu(xs-dx,us,others))/(2*d);
end
B = zeros(6,3);
%% Numerically determine B:
m=length(us);% the number of inputs
if(m==0) % if ther are no inputs make B zero
    B=zeros(n,1);
else 
    for i=1:m % loop through a small change in each input:
        du=zeros(m,1);
        du(i,1)=d;
        %disp('B matrix')
        B(:,i)=(ftxu(xs,us+du,others)-ftxu(xs,us-du,others))/(2*d);
    end
end


