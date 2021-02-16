function Q =  EP2ang_arr(EP)
%Change Euler parameters to angles
EP0 = EP(1);
EP1 = EP(2);
EP2 = EP(3);
EP3 = EP(4);

q1 = atan2( (2*(EP0*EP1 + EP2*EP3)),(1-2*(EP1^2 + EP2^2))  );

q2 = asin(2*(EP0*EP2 - EP3*EP1));

q3 = atan2( (2*EP0*EP3 + EP1*EP2),(1-2*(EP2^2 + EP3^2)) );

Q = [q1,q2,q3];
end



