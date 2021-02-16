function Moments = PIDcontroller(EPs,StatesInt,StatesDeriv,Kp,Kd,Kint)
ang = EP2ang(EPs).'; % Do linearization about the uler ANGLEs not params.


U1 = Kp*ang(1);
U2 = Kp*ang(2);
U3 = Kp*ang(3);

Moments = [U1,U2,U3];
end