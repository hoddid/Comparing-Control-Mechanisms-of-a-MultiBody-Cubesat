function angles_int = Angle_integrator(Input)

Angles = Input(1:3);
Timearr = Input(3:end);
if length(Timearr) > length(Angles(1))
     Time_tointegrate = Timearr((end-length(Angles(1))),end);
     int1 = trapz(Time_tointegrate,Angles(1));
     int2 = trapz(Time_tointegrate,Angles(2));
     int3 = trapz(Time_tointegrate,Angles(3));
 
else
    Time_tointegrate = Timearr;
    int1 = trapz(Time_tointegrate,Angles(1));
    int2 = trapz(Time_tointegrate,Angles(2));
    int3 = trapz(Time_tointegrate,Angles(3));
end

angles_int = [int1,int2,int3];

end