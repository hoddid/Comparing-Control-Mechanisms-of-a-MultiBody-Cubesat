/* The name of this program is oscar001.c */
/* Created by AUTOLEV on Tue Nov 17 12:45:53 2020 */
/* Autolev script for motion of Spacecraft (ideally 3U CubeSat OSCaR) using three 
  Reacion Wheel */
/* RW1, RW2, RW3 (associated with direction B1>, B2>, B3>, respectively) for attitude 
  control. */
/* In this initial treatment the control will be exercised as a reactive control 
  moment applied */
/* to each reaction wheel by the carrier B (for Body or Bus). perhaps in future 
  versions we will */
/* treat this as a controls problem where the angular velocity of each reaction 
  wheel with */
/* respect to the carrier B is controlled/prescribed by the control law. */
/* For this analysis we will be using Reaction wheels aligned with orthogonal 
  B directions B1>,B2>,B3>. */
/* Each reaction wheel many be independently sized by adjusting its moment of 
  inertial value about */
/* is spin axis. Because of its symmetry and the fact that the RW centers of mass 
  and spin axes */
/* with respect to the carrier B, the masses, and moment of initial values are 
  also considered as */
/* part of the B. (Talk to me about this if you need clarification). */
/* B1DIMENSION = The B1> dimension (width) of the three Unit (3U) OSCaR CubeSat. */
/* Nominally this will be 0.10 [m]. */
/* B2DIMENSION = The B2> dimension (depth) of the three Unit (3U) OSCaR CubeSat. */
/* Nominally this will be 0.10 [m]. */
/* B3DIMENSION = The B3> dimension (Length) of the three Unit (3U) OSCaR CubeSat. */
/* Nominally this will be 0.30 [m]. */
/* RBO1 = The B1> dimension of position vector from Body B Geometric Center (BGC) 
  to Body B CM (BO). */
/* RBO2 = The B2> dimension of position vector from Body B Geometric Center (BGC) 
  to Body B CM (BO). */
/* RBO3 = The B3> dimension of position vector from Body B Geometric Center (BGC) 
  to Body B CM (BO) */
/* B = CubeSat Body/Chassis/Carrier */
/* RWj = j-th Reaction Wheel. */

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "ctype.h"
#include "stdarg.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#define  _NAN      9.99999999999999E+305

void     eqns1     (double T, double VAR[], double VARp[], char boundary);
void     solve     (int size, double *a[], double b[], double x[]);
void     readf     (FILE *Fp, char message[], double *next, ...);
void     pgets     (FILE *Fp, char line[], double *x);
int      kutta     (void(*eqns)(double, double*, double*, char),
                    int numy, double y[], double *t, double integstp,
                    double abserr, double relerr, char com);
					
void     output    (FILE *Fptr[], double T);
void     writef    (FILE *Fp, char format[], ...);



double   CTORQUE1,CTORQUE2,CTORQUE3,IB11,IB12,IB22,IB23,IB31,IB33,IRW1,IRW2,IRW3,MASSB;
double   EP1,EP2,EP3,EP4,Q1,Q2,Q3,U1,U2,U3,U4,U5,U6,U7,U8,U9;
double   EP1p,EP2p,EP3p,EP4p,Q1p,Q2p,Q3p,U1p,U2p,U3p,U4p,U5p,U6p,U7p,U8p,U9p;
double   Pi,DEGtoRAD,RADtoDEG,z[46],_COEF[9][9],*COEF[9],RHS[9];


double CTORQUE1_mat,CTORQUE2_mat,CTORQUE3_mat,IB11_mat,IB12_mat,IB22_mat,IB23_mat,IB33_mat;
double IRW1_mat,IRW2_mat,IRW3_mat,MASSB_mat;


//Keep all these inputs. They are needed for MATLAB's MEX generation
//PLHS = pointer left hand side = outputs
//PRHS = pointer right hand size = inputs
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
FILE     *Fptr[6];
char     string[256];
int      iprint, printint, iloop;
double   T, tinitial,tfinal,integstp,abserr,relerr,_printint;
double   VAR[16];





//DECLARE MEX variables
//say inputs will be arrays
mxArray *a_in_m, *b_in_m, *c_out_m, *d_out_m;


double *a, *b, *c, *d;
 


//int dimx, dimy, numdims;


//associate inputs with pointers
a_in_m = mxDuplicateArray(prhs[0]);
//b_in_m = mxDuplicateArray(prhs[1]);
//^if i need a second input array


//DEFINE OUTPUTS: (size1,size2,real or imaginary nums)
c_out_m = plhs[0] = mxCreateDoubleMatrix(1,32,mxREAL);
d_out_m = plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

//associate pointers for inputs and putputs
a = mxGetPr(a_in_m);
b = mxGetPr(b_in_m);
c = mxGetPr(c_out_m);
d = mxGetPr(d_out_m);







/* Initialize COEF pointers to beginning of each row */
for(iloop=0;  iloop<9;  iloop++)  COEF[iloop] = &(_COEF[iloop][0]);

/* Open input and output files */

for(iloop=0;  iloop<=5;  iloop++)
  {
  if( !iloop ) strcpy(string, "oscar001.in");
  else sprintf(string, "oscar001.%d", iloop);
  if( (Fptr[iloop] = fopen(string, iloop ? "w" : "r")) == NULL)
    {//mexPrintf("Error: unable to open file %s\n", string); // commented for speed catches a lot;
	}
  }
 
/* Read message from input file */
//readf(Fptr[0],string,NULL);

/* Read values of constants from input file */
//readf(Fptr[0],NULL,&CTORQUE1,&CTORQUE2,&CTORQUE3,&IB11,&IB12,&IB22,&IB23,&IB31,&
//  IB33,&IRW1,&IRW2,&IRW3,&MASSB,NULL);
  
// mexPrintf("element 1 is %f %f %f %f %f %f %f %f %f %f %f %f\n",a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11]);
CTORQUE1 = a[0];
//mexPrintf("ctorque mat is %f, but read val is %f \n",CTORQUE1_mat,CTORQUE1);
CTORQUE2 = a[1];
CTORQUE3 = a[2];
IB11 = a[3];
IB12 = a[4];
IB22 = a[5];
IB23 = a[6];
IB31 = a[7];
IB33 = a[8];
IRW1 = a[9];
IRW2 = a[10];
IRW3 = a[11];
MASSB = a[12];
    
  

/* Read the initial value of each variable from input file */
//readf(Fptr[0],NULL,&EP1,&EP2,&EP3,&EP4,&Q1,&Q2,&Q3,&U1,&U2,&U3,&U4,&U5,&U6,&U7,&
 // U8,&U9,NULL);
  
  
EP1 = a[13];
EP2 = a[14];
EP3 = a[15];
EP4 = a[16];
Q1 = a[17];
Q2 = a[18];
Q3 = a[19];
U1 = a[20];
U2 = a[21];
U3 = a[22];
U4 = a[23];
U5 = a[24];
U6 = a[25];
U7 = a[26];
U8 = a[27];
U9 = a[28];
 
  

/* Read integration parameters from input file */
//readf(Fptr[0],NULL,&tinitial,&tfinal,&integstp,&_printint,&abserr,&relerr,NULL);

tinitial = a[29];
tfinal = a[30];
integstp = a[31];
printint = a[32];
abserr = a[33];
relerr = a[34];



printint = (int)_printint;







/* Write heading(s) to output file(s) */
fprintf(stdout,  " FILE: oscar001.1\n\n *** %s\n\n", string);
fprintf(stdout,  "        T          P_O_BO[1]      P_O_BO[2]      P_O_BO[3]      N_B[1,1]       N_B[1,2]       N_B[1,3]       N_B[2,1]       N_B[2,2]       N_B[2,3]       N_B[3,1]       N_B[3,2]       N_B[3,3]\n"
                 "     (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n" );
fprintf(Fptr[1], " FILE: oscar001.1\n\n *** %s\n\n", string);
fprintf(Fptr[1], "        T          P_O_BO[1]      P_O_BO[2]      P_O_BO[3]      N_B[1,1]       N_B[1,2]       N_B[1,3]       N_B[2,1]       N_B[2,2]       N_B[2,3]       N_B[3,1]       N_B[3,2]       N_B[3,3]\n"
                 "     (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n" );
fprintf(Fptr[2], " FILE: oscar001.2\n\n *** %s\n\n", string);
fprintf(Fptr[2], "        T             Q1             Q2             Q3\n"
                 "     (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n" );
fprintf(Fptr[3], " FILE: oscar001.3\n\n *** %s\n\n", string);
fprintf(Fptr[3], "        T             U1             U2             U3\n"
                 "     (UNITS)         (M/S)          (M/S)          (M/S)\n\n" );
fprintf(Fptr[4], " FILE: oscar001.4\n\n *** %s\n\n", string);
fprintf(Fptr[4], "        T             U4             U5             U6\n"
                 "     (UNITS)        (RAD/S)        (RAD/S)        (RAD/S)\n\n" );
fprintf(Fptr[5], " FILE: oscar001.5\n\n *** %s\n\n", string);
fprintf(Fptr[5], "        T             U7             U8             U9\n"
                 "     (UNITS)        (RAD/S)        (RAD/S)        (RAD/S)\n\n" );

/* Degrees to radians conversion */
  Pi       = 3.14159265358979323846264338327950288419716939937510582;
  DEGtoRAD = 0.01745329251994329576923690768488612713442871888541725;
  RADtoDEG = 57.2957795130823208767981548141051703324054724665643215;

/* Evaluate constants */
  z[40] = IB11 + IRW1;
  z[42] = IB22 + IRW2;
  z[44] = IB33 + IRW3;

/* Initialize time, print counter, variables array for integrator */
  T      = tinitial;
  iprint = 0;
  VAR[0] = EP1;
  VAR[1] = EP2;
  VAR[2] = EP3;
  VAR[3] = EP4;
  VAR[4] = Q1;
  VAR[5] = Q2;
  VAR[6] = Q3;
  VAR[7] = U1;
  VAR[8] = U2;
  VAR[9] = U3;
  VAR[10] = U4;
  VAR[11] = U5;
  VAR[12] = U6;
  VAR[13] = U7;
  VAR[14] = U8;
  VAR[15] = U9;

/* Initalize numerical integrator, with call to eqns1 at T=TINITIAL */
kutta(eqns1, 16, VAR, &T, integstp, abserr, relerr, 0);

/* Numerically integrate; print results */
while(1)
  {
  if( tfinal>=tinitial && T+.01*integstp>=tfinal) iprint=-7;
  if( tfinal<=tinitial && T+.01*integstp<=tfinal) iprint=-7;
  if( iprint <= 0 )
    {
	//THIS IS WHERE TO PUT THE OUTPUT VARIABLES
	// 	
	c[0] = IB11;
	c[1] = IB12;
    c[2] = IB22; 
    c[3] = IB23;
	c[4] = IB31;
	
    c[5] = IB33;
    c[6] = IRW1;
    c[7] = IRW2;
    c[8] = IRW3;
    c[9] = MASSB;
	c[10] = EP1;
	c[11] = EP2;
	c[12] = EP3;
	c[13] = EP4;
	c[14] = Q1;
	c[15] = Q2;
	c[16] = Q3;
	c[17] = U1;
	c[18] = U2;
	c[19] = U3;
	c[20] = U4;
	c[21] = U5;
	c[22] = U6;
	c[23] = U7;
	c[24] = U8;
	c[25] = U9;
	c[26] = tinitial;
	//mexPrintf("tFINAL is %d \n", tfinal);
	//mexPrintf("tinit is %d \n", tinitial);
	//mexPrintf("integstep is %5d \n", integstp);
    c[27] = tfinal;
    c[28] = integstp;
    c[29] = printint;
    c[30] = abserr; 
    c[31] = relerr;
		
		
		
    output(Fptr, T);
    if( iprint == -7 ) break;
    iprint = printint;
    }
  if( !kutta(eqns1, 16, VAR, &T, integstp, abserr, relerr, 1) )
    {
    output(Fptr, T);
    for(iloop=0;  iloop<=5;  iloop++)
      fputs( "\n\tError: Numerical integration failed to converge\n",
              iloop ? Fptr[iloop]:stdout);
    break;
    }
  iprint--;
  }

/* Inform user of input and output filename(s) */
puts( "\n Input is in the file oscar001.in" );
puts( "\n Output is in the file(s) oscar001.i  (i=1, ..., 5)" );
puts( "\n The output quantities and associated files are listed in file "
    "oscar001.dir\n" );
//return 0;
}


/* ................................ EQNS1 ............................. */
void     eqns1        (double T, double VAR[], double VARp[], char boundary)
{

/* Update variables after integration step */
  EP1 = VAR[0];
  EP2 = VAR[1];
  EP3 = VAR[2];
  EP4 = VAR[3];
  Q1 = VAR[4];
  Q2 = VAR[5];
  Q3 = VAR[6];
  U1 = VAR[7];
  U2 = VAR[8];
  U3 = VAR[9];
  U4 = VAR[10];
  U5 = VAR[11];
  U6 = VAR[12];
  U7 = VAR[13];
  U8 = VAR[14];
  U9 = VAR[15];

  Q1p = U1;
  Q2p = U2;
  Q3p = U3;
  EP1p = 0.5*EP2*U6 + 0.5*EP4*U4 - 0.5*EP3*U5;
  EP2p = 0.5*EP3*U4 + 0.5*EP4*U5 - 0.5*EP1*U6;
  EP3p = 0.5*EP1*U5 + 0.5*EP4*U6 - 0.5*EP2*U4;
  EP4p = -0.5*EP1*U4 - 0.5*EP2*U5 - 0.5*EP3*U6;
  z[10] = IB11*U4;
  z[11] = IB12*U5;
  z[12] = IB31*U6;
  z[13] = IB12*U4;
  z[14] = IB22*U5;
  z[15] = IB23*U6;
  z[16] = IB31*U4;
  z[17] = IB23*U5;
  z[18] = IB33*U6;
  z[19] = U4*z[13] + U4*z[14] + U4*z[15] - U5*z[10] - U5*z[11] - U5*z[12];
  z[20] = U6*z[10] + U6*z[11] + U6*z[12] - U4*z[16] - U4*z[17] - U4*z[18];
  z[21] = U5*z[16] + U5*z[17] + U5*z[18] - U6*z[13] - U6*z[14] - U6*z[15];
  z[22] = U4 + U7;
  z[25] = IRW1*z[22];
  z[26] = U5*z[25];
  z[27] = U6*z[25];
  z[28] = U5 + U8;
  z[31] = IRW2*z[28];
  z[32] = U4*z[31];
  z[33] = U6*z[31];
  z[34] = U6 + U9;
  z[37] = IRW3*z[34];
  z[38] = U4*z[37];
  z[39] = U5*z[37];
  z[41] = z[21] + z[39] - z[33];
  z[43] = z[20] + z[27] - z[38];
  z[45] = z[19] + z[32] - z[26];

  COEF[0][0] = -MASSB;
  COEF[0][1] = 0;
  COEF[0][2] = 0;
  COEF[0][3] = 0;
  COEF[0][4] = 0;
  COEF[0][5] = 0;
  COEF[0][6] = 0;
  COEF[0][7] = 0;
  COEF[0][8] = 0;
  COEF[1][0] = 0;
  COEF[1][1] = -MASSB;
  COEF[1][2] = 0;
  COEF[1][3] = 0;
  COEF[1][4] = 0;
  COEF[1][5] = 0;
  COEF[1][6] = 0;
  COEF[1][7] = 0;
  COEF[1][8] = 0;
  COEF[2][0] = 0;
  COEF[2][1] = 0;
  COEF[2][2] = -MASSB;
  COEF[2][3] = 0;
  COEF[2][4] = 0;
  COEF[2][5] = 0;
  COEF[2][6] = 0;
  COEF[2][7] = 0;
  COEF[2][8] = 0;
  COEF[3][0] = 0;
  COEF[3][1] = 0;
  COEF[3][2] = 0;
  COEF[3][3] = -z[40];
  COEF[3][4] = -IB12;
  COEF[3][5] = -IB31;
  COEF[3][6] = -IRW1;
  COEF[3][7] = 0;
  COEF[3][8] = 0;
  COEF[4][0] = 0;
  COEF[4][1] = 0;
  COEF[4][2] = 0;
  COEF[4][3] = -IB12;
  COEF[4][4] = -z[42];
  COEF[4][5] = -IB23;
  COEF[4][6] = 0;
  COEF[4][7] = -IRW2;
  COEF[4][8] = 0;
  COEF[5][0] = 0;
  COEF[5][1] = 0;
  COEF[5][2] = 0;
  COEF[5][3] = -IB31;
  COEF[5][4] = -IB23;
  COEF[5][5] = -z[44];
  COEF[5][6] = 0;
  COEF[5][7] = 0;
  COEF[5][8] = -IRW3;
  COEF[6][0] = 0;
  COEF[6][1] = 0;
  COEF[6][2] = 0;
  COEF[6][3] = -IRW1;
  COEF[6][4] = 0;
  COEF[6][5] = 0;
  COEF[6][6] = -IRW1;
  COEF[6][7] = 0;
  COEF[6][8] = 0;
  COEF[7][0] = 0;
  COEF[7][1] = 0;
  COEF[7][2] = 0;
  COEF[7][3] = 0;
  COEF[7][4] = -IRW2;
  COEF[7][5] = 0;
  COEF[7][6] = 0;
  COEF[7][7] = -IRW2;
  COEF[7][8] = 0;
  COEF[8][0] = 0;
  COEF[8][1] = 0;
  COEF[8][2] = 0;
  COEF[8][3] = 0;
  COEF[8][4] = 0;
  COEF[8][5] = -IRW3;
  COEF[8][6] = 0;
  COEF[8][7] = 0;
  COEF[8][8] = -IRW3;
  RHS[0] = 0;
  RHS[1] = 0;
  RHS[2] = 0;
  RHS[3] = z[41];
  RHS[4] = z[43];
  RHS[5] = z[45];
  RHS[6] = -CTORQUE1;
  RHS[7] = -CTORQUE2;
  RHS[8] = -CTORQUE3;
  solve(9,COEF,RHS,VARp);

/* Update variables after uncoupling equations */
  U1p = VARp[0];
  U2p = VARp[1];
  U3p = VARp[2];
  U4p = VARp[3];
  U5p = VARp[4];
  U6p = VARp[5];
  U7p = VARp[6];
  U8p = VARp[7];
  U9p = VARp[8];

/* Update derivative array prior to integration step */
  VARp[0] = EP1p;
  VARp[1] = EP2p;
  VARp[2] = EP3p;
  VARp[3] = EP4p;
  VARp[4] = Q1p;
  VARp[5] = Q2p;
  VARp[6] = Q3p;
  VARp[7] = U1p;
  VARp[8] = U2p;
  VARp[9] = U3p;
  VARp[10] = U4p;
  VARp[11] = U5p;
  VARp[12] = U6p;
  VARp[13] = U7p;
  VARp[14] = U8p;
  VARp[15] = U9p;

}


/* ................................ OUTPUT ............................. */
void     output (FILE *Fptr[], double T)
{
int      i1;

/* Evaluate output quantities */
  z[1] = 1 - 2*pow(EP2,2) - 2*pow(EP3,2);
  z[2] = 2*EP1*EP2 - 2*EP3*EP4;
  z[3] = 2*EP1*EP3 + 2*EP2*EP4;
  z[4] = 2*EP1*EP2 + 2*EP3*EP4;
  z[5] = 1 - 2*pow(EP1,2) - 2*pow(EP3,2);
  z[6] = 2*EP2*EP3 - 2*EP1*EP4;
  z[7] = 2*EP1*EP3 - 2*EP2*EP4;
  z[8] = 2*EP1*EP4 + 2*EP2*EP3;
  z[9] = 1 - 2*pow(EP1,2) - 2*pow(EP2,2);

/* Write output to screen and to output file(s) */
  writef(stdout, " %- 14.6E", T,Q1,Q2,Q3,z[1],z[2],z[3],z[4],z[5],z[6],z[7],
  z[8],z[9],_NAN);
  writef(Fptr[1]," %- 14.6E", T,Q1,Q2,Q3,z[1],z[2],z[3],z[4],z[5],z[6],z[7],
  z[8],z[9],_NAN);
  writef(Fptr[2]," %- 14.6E", T,Q1,Q2,Q3,_NAN);
  writef(Fptr[3]," %- 14.6E", T,U1,U2,U3,_NAN);
  writef(Fptr[4]," %- 14.6E", T,U4,U5,U6,_NAN);
  writef(Fptr[5]," %- 14.6E", T,U7,U8,U9,_NAN);
}


/*................................... READF ................................*/
void     readf     (FILE *Fp, char message[], double *next, ...) 
{
va_list  args;                                   /* Variable Argument List  */
char     line[256];

if( message )
  {
  pgets(Fp,line,NULL);   pgets(Fp,line,NULL);   pgets(Fp,message,NULL);
  pgets(Fp,line,NULL);   pgets(Fp,line,NULL);
  }
else                                             /* Get number from file    */
  {
  for( va_start(args,next);  next;  next=va_arg(args,double *) )
    pgets(Fp,line,next);
  va_end(args);                        /* Help function make normal return  */
  }
pgets(Fp,line,NULL);                   /* Always get a newline at the end   */
}


/*................................... PGETS ................................*/
void     pgets     (FILE *Fp, char line[], double *x)
{
static   long      line_number=0;
char               *s1;

line_number++;
if( !fgets(line,256,Fp) )
  {
//  mexPrintf("\n Unable to read line %ld of input file."
 //        "\n Premature end of file found while reading input file\n", line_number);
  
  }
if( !x ) return;                                 /* No number expected      */
if( strlen(line) >= 60 )                         /* Expecting number here   */
  {
  *x = strtod(line+59,&s1);                      /* Get double number       */
  while( isspace(*s1) )  s1++;                   /* Skip white space to end */
  if( !*s1 ) return;                             /* Check for valid number  */
  }
//mexPrintf("\n An error occured while reading line %ld of the input file."
  //     "\n The program was expecting to find one double precision number"
    //   "\n beginning with the 60th character in the following line:\n\n%s\n",
      //     line_number, line );

}


/* .................................. WRITEF .............................. */
void     writef       (FILE *Fp, char format[], ...)
{
va_list  args;                                   /* Variable Argument List  */
double   next;                                   /* Current Place in List   */

va_start(args,format);                           /* args start after format */
while( (next=va_arg(args,double)) != _NAN )  fprintf(Fp, format, next);
va_end(args);                                    /* End of variable list    */
fprintf(Fp, "\n");                               /* End with newline        */
}


/*****************************************************************************
C**                                                                         **
C** PURPOSE  Solves a set of first order ordinary differential equations    **
C**          of the form dy(i)/dt = F(t,y(1), ..., y(numeqns) (i = 1,       **
C**          ..., numeqns)                                                  **
C**                                                                         **
C** INPUT                                                                   **
C**    eqns: Subroutine that evaluates dy(i)/dt (i = 1, ..., numeqns), the  **
C**          first derivatives of y(1), ..., y(numeqns) with respect to t   **
C**                                                                         **
C** numeqns: The number of differential equations to be solved              **
C**                                                                         **
C**       y: One-dimensional array whose elements are y(1), ..., y(numeqns) **
C**                                                                         **
C**       t: Independent variable                                           **
C**                                                                         **
C** integstp: Maximum integration stepsize                                  **
C**                                                                         **
C**  abserr: Allowable absolute error in y(i)  (i=1, ..., numeqns)          **
C**                                                                         **
C**  relerr: Allowable relative error in y(i)  (i=1, ..., numeqns)          **
C**                                                                         **
C**     com: When com = 2, the Kutta-Merson algorithm (L. Fox, Numerical    **
C**          Solutions of Ordinary and Partial Differential Equations,      **
C**          Palo Alto: Addison-Wesley, 1962, pp. 24-25) is employed to     **
C**          perform the numerical solution of the differential equations.  **
C**          Accordingly, dy(i)/dt (i = 1, ..., numeqns) are evaluated at   **
C**          every integration boundary, including those at Tinitial,       **
C**          Tfinal, and ones created when integstp is halved to satisfy    **
C**          the requirements imposed by abserr and relerr.  Integration    **
C**          is self-starting at each boundary, and the occurrence, at      **
C**          boundaries, of discontinuities in derivatives does not lead    **
C**          to failure of the integration process.                         **
C**                                                                         **
C**          When com = 1, a modified form of the Kutta-Merson algorithm    **
C**          is employed.  It is nearly 20% faster than the one used when   **
C**          com = 2 because no recalculation of derivatives at inte-       **
C**          gration boundaries between Tinitial and Tfinal takes place.    **
C**          Integration is self-starting at Tinitial and Tfinal only.      **
C**          Integration may fail if any of dy(i)/dt (i = 1, ..., numeqns)  **
C**          is discontinuous between Tinitial and Tfinal.                  **
C**                                                                         **
C**          When com = 0, the function eqns is called and dy(i)/dt         **
C**          (i = 1, ..., numeqns) are evaluated, but no integration        **
C**          is performed.                                                  **
C**                                                                         **
C** OUTPUT                                                                  **
C**          The value of t+integstp is returned in t, and the values of    **
C**          y(i) at t+integstp are returned in y.                          **
C**                                                                         **
C** SOURCE                                                                  **
C**          Copyright 1995 by Paul C. Mitiguy, Thomas R. Kane, David A.    **
C**          Levinson, and David B. Schaechter.  Permission is granted      **
C**          to copy, modify, and distribute this subroutine, provided      **
C**          that this copyright notice appear.                             **
C**                                                                         **
C****************************************************************************/
int      kutta        ( void(*eqns)(double, double*, double*, char),
                        int numy, double y[], double *t,
                        double integstp, double abserr, double relerr, char com )
{
static double  f0[100], f1[100], f2[100], y1[100], y2[100];
static int     numcuts = 20;                     /* Max # cuts of integstp  */
static double  hc = 0;                           /* Last value of stepsize  */
char           entry = 1;                        /* Just entered routine    */
char           stepdouble;                       /* Double the stepsize     */
double         tfinal = *t + integstp;           /* Time at end of full step*/
double         error, test;                      /* Testing error criterion */
double         tt, h, h2, h3, h6, h8;
int            i;

if( !com ) { (*eqns)(*t,y,f0,1);  return 1;}     /* Fill array f0 and return*/
if( numy == 0)  { hc = integstp;  return 1;}     /* Check for initial entry */
if( integstp == 0)  return 0;                    /* Cannot integrate forward*/
if( hc*integstp < 0 ) hc = -hc;                  /* Integrate backward      */        
else if( hc == 0 )    hc = integstp;             /* Maybe initial entry     */
h  = hc;                                         /* Current stepsize        */
tt = *t + h;                                     /* Terminal time this step */
*t = tfinal;                                     /* Return updated t value  */

beginning:
while( tt+h != tt )                              /* Check round-off problems*/
  {
  h2 = h * 0.5;                                            /* Half    of h  */
  h3 = h / 3.0;                                            /* Third   of h  */
  h6 = h / 6.0;                                            /* Sixth   of h  */
  h8 = h * 0.125;                                          /* Eighth  of h  */
  if( com==2 || entry)                                     
 {(*eqns)( tt-h,     y, f0, 1 );   entry=0; }              /* Boundary here */
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h3*f0[i];
  (*eqns)( tt-2*h3, y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h6*(f0[i] + f1[i]);
  (*eqns)( tt-2*h3, y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h8*(f0[i] + 3*f1[i]);
  (*eqns)( tt-h2,   y1, f2, 0 );
  for( i=0;  i<numy;  i++)  y1[i] = y[i] + h2*(f0[i] - 3*f1[i] + 4*f2[i]);
  (*eqns)( tt,      y1, f1, 0 );
  for( i=0;  i<numy;  i++)  y2[i] = y[i] + h6*(f0[i] + 4*f2[i] + f1[i]);
  stepdouble = 1;                                /* Assume need to double   */
  for(i=0;  i<numy;  i++)                        /* Check all equations     */
    {
    error = fabs( y1[i] - y2[i] ) * 0.2;         /* Error in estimate       */
    test  = fabs( y1[i] ) * relerr;              /* For relative error      */
    if( error >= test  &&  error >= abserr )     /* Error criterion not met?*/
      {                                          /* Restart w/ half stepsize*/
      hc = h = h2;                               /* Halve the stepsize      */
      tt -= h2;                                  /* Change terminal time    */
      if( numcuts-- > 0 )  goto beginning;       /* Back to beginning       */
      mexPrintf("\n THE STEPSIZE HAS BEEN HALVED TOO MANY TIMES; T = %15.8E"
             "\n NUMERICAL INTEGRATION FAILED TO CONVERGE.\n", *t=tt-h );
      (*eqns)(*t,y,f0,0);                        /* Fill for error display  */
      return 0;
      }
    if( stepdouble && 64*error>test && 64*error>abserr ) stepdouble = 0;
    }
  for(i=0;  i<numy;  i++)  y[i] = y2[i];         /* Update y for next step  */
  if( stepdouble && fabs(h+h)<=fabs(integstp) && fabs(tt+h+h)<=fabs(tfinal) )
    {hc=(h+=h);  numcuts++;}                     /* Double the stepsize     */
  if( tt == tfinal )                             /* End of integration      */
    { (*eqns)(tfinal,y,f0,2);  return 1;}        /* Derivatives at tfinal   */
  tt += h;                                       /* Update terminal time    */
  /*** If next jump puts tt past or very close to tfinal, adjust h and tt ***/
  if( (h>0 && tt>tfinal-0.1*h) || (h<0 && tt<tfinal-0.1*h) )
    { h = tfinal-(tt-h);  tt = tfinal; }         /* Adjust for last jump    */
  if(com == 1)                                   /* Approx this derivative  */
    for(i=0;  i<numy;  i++)  f0[i] = f1[i];
 }
mexPrintf("\nTHE STEPSIZE OF %15.7E IS TOO SMALL RELATIVE TO THE TERMINAL TIME"
       "\nOF %15.7E.  INTEGRATION HALTED BECAUSE OF NUMERICAL ROUND-OFF.\n"
       "\nTHE STEPSIZE MAY HAVE BEEN CUT TOO MANY TIMES.\n\n", h, *t=tt);
(*eqns)(*t,y,f0,0);                              /* Fill for error display  */
return 0;
}



/**************************************************************************** 
C**                                                                        ** 
C** PURPOSE  The matrix equation a x = b is solved for x, where a is an    ** 
C**          n by n matrix, and x and b are n by 1 matrices.               ** 
C**                                                                        ** 
C** INPUT                                                                  **
C**       N: n                                                             ** 
C**                                                                        ** 
C**       A: an N by N double precision array whose elements are those     **      
C**          of the matrix a                                               ** 
C**                                                                        ** 
C**       B: an N by 1 double precision array whose elements are those     **      
C**          of the matrix b                                               ** 
C**                                                                        ** 
C** OUTPUT                                                                 ** 
C**       X: an N by 1 double precision array whose elements are those     **
C**          of the matrix x                                               ** 
C**                                                                        ** 
C***************************************************************************/ 
void           solve  (int n, double *A[], double B[], double X[])
{
static double  scales[100];                      /* Row scaling factors     */
int            i, j, k, swapi;
double         rowmax, big, size, pivot, ratio, sum, *swapa, swapx;
                                                 
for(i=0;  i<n;  i++)                             /*** Begin decomposition ***/
  {
  for(rowmax=0.0, j=0;  j<n;  j++)               /* Check for zero row      */
    if( rowmax < fabs(A[i][j]) )  rowmax = fabs(A[i][j]);
  if( rowmax == 0.0 )                            /* Zero row found          */
    {mexPrintf("ALL ELEMENTS IN ROW %d OF COEF ARE ZEROS\n", i+1);   }
  scales[i] = 1.0 / rowmax;                      
  X[i] = B[i];                                   /* Leave B matrix unchanged*/
  }
for(k=0;  k<n-1;  k++)                           /* Until 1 column is Left  */
  {
  for(big=0.0, i=k;  i<n;  i++)                  /* Check remaining rows    */
    {
    size = fabs( A[i][k] ) * scales[i];          /* Relative size this col  */
    if( size > big )  {swapi=i;   big=size;}     /* Largest relative entry  */
    }
  if( big == 0.0 )                               /* Zero pivot found        */
    {
    mexPrintf("A PIVOT ELEMENT ENCOUNTERED IN THE DECOMPOSITION OF COEF IS ZERO"
           "\nCOEFFICIENT MATRIX IS SINGULAR\n");
    
    }
  if( swapi != k )                               /* Switch rows of A and X  */
    {
    swapa=A[k];  A[k]=A[swapi];  A[swapi]=swapa; /* Change row pointers     */
    swapx=X[k];  X[k]=X[swapi];  X[swapi]=swapx; /* Change X values         */
    }
  pivot = A[k][k];                               /* Value of pivot          */
  for(i=k+1;  i<n;  i++)                         /* Change lower rows       */
    {
    ratio = A[i][k] / pivot;                     /* Multiplicative factor   */
    X[i] -= ratio * X[k];                        /* Modify elements of X    */
    for(j=k+1;  j<n;  j++)  A[i][j] -= ratio * A[k][j];  
    }
  }
if( A[n-1][n-1] == 0.0 )                         /* Check last pivot        */
  {
  mexPrintf("A PIVOT ELEMENT ENCOUNTERED IN THE DECOMPOSITION OF COEF IS ZERO"
         "\nCOEFFICIENT MATRIX IS SINGULAR\n");

  }
X[n-1] = X[n-1] / A[n-1][n-1];                   /* 1st equation of upper   */
for(i=n-2;  i>=0;  i--)                          /* Solve upper * X = Z     */
  {
  for(sum=0.0, j=i+1;  j<n;  j++)    sum += A[i][j] * X[j];
  X[i] = (X[i] - sum) / A[i][i];
  }
}



