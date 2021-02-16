%% Autolev script for motion of Spacecraft (ideally 3U CubeSat OSCaR) using three Reacion Wheel
%% RW1, RW2, RW3 (associated with direction B1>, B2>, B3>, respectively) for attitude control.
%% In this initial treatment the control will be exercised as a reactive control moment applied
%% to each reaction wheel by the carrier B (for Body or Bus). perhaps in future versions we will 
%% treat this as a controls problem where the angular velocity of each reaction wheel with 
%% respect to the carrier B is controlled/prescribed by the control law.

%% For this analysis we will be using Reaction wheels aligned with orthogonal B directions B1>,B2>,B3>. 
%% Each reaction wheel many be independently sized by adjusting its moment of inertial value about
%% is spin axis. Because of its symmetry and the fact that the RW centers of mass and spin axes 
%% with respect to the carrier B, the masses, and moment of initial values are also considered as
%% part of the B. (Talk to me about this if you need clarification).
% 

PAUSE 3 
AUTOZ ON

CONSTANTS B1DIMENSION, B2DIMENSION, B3DIMENSION, RBO1, RBO2, RBO3
%% B1DIMENSION = The B1> dimension (width) of the three Unit (3U) OSCaR CubeSat. 
%%               Nominally this will be 0.10 [m].
%% B2DIMENSION = The B2> dimension (depth) of the three Unit (3U) OSCaR CubeSat. 
%%               Nominally this will be 0.10 [m].
%% B3DIMENSION = The B3> dimension (Length) of the three Unit (3U) OSCaR CubeSat. 
%%               Nominally this will be 0.30 [m].

%% RBO1 = The B1> dimension of position vector from Body B Geometric Center (BGC) to Body B CM (BO). 
%% RBO2 = The B2> dimension of position vector from Body B Geometric Center (BGC) to Body B CM (BO).
%% RBO3 = The B3> dimension of position vector from Body B Geometric Center (BGC) to Body B CM (BO)


VARIABLES Q{3}',EP{4}', U{9}'
% Q1 = X (n1>) Location of Body B geometric Center ( ideal coincident with the Body B Center of Mass)
% Q2 = Y (N2>) Location of Body B geometric Center ( ideal coincident with the Body B Center of Mass
% Q3 = Z (N3>) Location of Body B geometric Center ( ideal coincident with the Body B Center of Mass


% ep{4} EULER PARAMETERs


BODIES B, RW1, RW2, RW3
%% B = CubeSat Body/Chassis/Carrier
%% RWj = j-th Reaction Wheel.

Newtonian N % Define Newtonian Frame to be "N"


POINTS O % Absolute Original 
POINTS BGC % Body B Geometric Center (not generally co-incident with BO) 

% Define names to be used for the respective masses of each body
% Mass Properties Calculations (Remember the mass of the RWs are already
%      considered in the overall mass of OSCaR B)
MASS B=MASSB, RW1=0, RW2=0, RW3=0


INERTIA B, IB11, IB22, IB33, IB12, IB23, IB31
INERTIA RW1, IRW1, 0, 0, 0, 0, 0
INERTIA RW2, 0, IRW2, 0, 0, 0, 0
INERTIA RW3, 0, 0, IRW3, 0, 0, 0


% Important Position Vectors
P_O_BO> = Q1*N1> + Q2*N2> + Q3*N3> % Position of Body/Bus/Carrier CM for analysis and Animake
P_BGC_BO> = RBO1*B1> + RBO2*B2> + RBO3*B3> % position vector of Body B CM (BO) relative to 
   % Body B Geometric Center BGC


% Orientation relations (use Euler Parameters for detemining and tracking orientation)
DIRCOS(N,B,EULER,EP1,EP2,EP3,EP4)
SIMPROT(B,RW1,1,0) % RW1 rotates relative to B about direction B1>. Angle is set to 0 because
                   % RW1 symmetry about B1> axis, and the rotation angle does not need to be 
                   % tracked in this applicaion.(this also simplifies and accelerates program)
SIMPROT(B,RW2,2,0) % RW2 rotates relative to B about direction B2>. Angle is set to 0 because
                   % RW2 symmetry about B2> axis, and the rotation angle does not need to be 
                   % tracked in this applicaion.(this also simplifies and accelerates program)
SIMPROT(B,RW3,3,0) % RW3 rotates relative to B about direction B1>. Angle is set to 0 because
                   % RW3 symmetry about B3> axis, and the rotation angle does not need to be 
                   % tracked in this applicaion.(this also simplifies and accelerates program)

ANIMATE(N,O,B) % generate data for animake animation

%  Angular velocities
W_B_N> = U4*B1> + U5*B2> + U6*B3>
W_RW1_B> = U7*B1>
W_RW2_B> = U8*B2>
W_RW3_B> = U9*B3>


% Velocities
V_O_N> = 0> % Point O is inertially fixed
V_BO_N> = U1*N1> + U2*N2> + U3*N3> % Velocity of CM of the body/Chassis/carrier B
V_RW1O_N> = V_BO_N> % these commands are provided to make AutoLev happy but will 
                    % not ultimately contribute because MRW1 = 0
V_RW2O_N> = V_BO_N>
V_RW3O_N> = V_BO_N>

% Kinematical Differential equations
Q1'=U1
Q2'=U2
Q3'=U3

% Tell AutoLev that the three Control Torques will be specified quantities
CONSTANTS CTORQUE1, CTORQUE2, CTORQUE3 

%Apply Control moments to Reaction Wheels by Carrier
TORQUE(B/RW1,CTORQUE1*B1>) % Control Torque applied to RW1
TORQUE(B/RW2,CTORQUE2*B2>) % Control Torque applied to RW2
TORQUE(B/RW3,CTORQUE3*B3>) % Control Torque applied to RW3



KINDIFFS(N,B,EULER,EP1,EP2,EP3,EP4)

ZERO = FR() + FRSTAR()

KANE()


%UNITSYSTEM KG,METER,SECOND

UNITS [B1DIMENSION, B2DIMENSION, B3DIMENSION, RBO1, RBO2, RBO3]= METRES, [MASSB]= KG, [U1, U2, U3] = m/s, [U4,U5,U6,U7,U8,U9]=Rad/s


OUTPUT T, Q1, Q2, Q3
OUTPUT T, U1, U2, U3
OUTPUT T, U4, U5, U6
OUTPUT T, U7, U8, U9

CODE DYNAMICS() OSCAR001.C,SUBS