function Q =  EP2ang_interp(EP)
% Use Builtin functions to convert from Euler Parameters to Euler angles
EP0 = EP(1);
EP1 = EP(2);
EP2 = EP(3);
EP3 = EP(4);
[q1,q2,q3] = quat2angle([EP0 EP1 EP2 EP3]);

Q = [q1;q2;q3];
end