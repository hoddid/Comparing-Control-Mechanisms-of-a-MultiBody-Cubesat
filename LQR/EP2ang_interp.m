function Q =  EP2ang_interp(EP)

EP0 = EP(1);
EP1 = EP(2);
EP2 = EP(3);
EP3 = EP(4);
[q1,q2,q3] = quat2angle([EP0 EP1 EP2 EP3]);

%Correct for overflow (Keeps angles strictly positive, thus ensuring controller always takes 'closest' route)
%if q1 < 0
%    q1 = q1+2*pi;
%end
%if q2 <0
%    q2 = q2+2*pi;
%end
%if q3 <0
%    q3 = q3+2*pi;
%end
    
Q = [q1;q2;q3];
end