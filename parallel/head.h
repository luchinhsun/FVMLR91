/* 2D LR91 */
/* Operator splitting :
FVM for diffusion term
forward Euler method for ODE */
/* ref: An Advanced Algorithm for Solving Partial Differential Equation in Cardiac Conduction. 1999. */
/* Xiang Zhou, 2018/2/23 */
/* ChinHsun Lu, 2018/03 */

#include <math.h>
#include <stdio.h>
#include <sys/timeb.h>
#include <malloc.h>
#include <string.h>

//-*******mesh parameters for LR91 *******
#define length_x 1.0//cm, the size of the tissue (cm)*(cm)
#define length_y 1.0
#define np 6529//6529;//1665;//433;//117;//34; //number of Vertices
#define nt 12800//12800;//3200;//800;//200;//50; //number of Triangles

#define D 0.001//D: diffusion coefficient cm^2/ms

/* Time Step */
#define dt 0.005 // Time step (ms)
#define cutcount (40 / dt)

/* Total Current and Stimulus */
#define st -80.0 // Constant Stimulus (uA/cm^2)
//double tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
//double it[nt] = { 0 };; // Total current (uA/cm^2)
//int stim_i[nt] = { 0 }; //{0,14,2,13,30,11,29,20,49,10,5, 3,48,6,9,8,27,16,50};// apply the stimulus on these elements 

/* Terms for Solution of Conductance and Reversal Potential */
#define R 8314.0 // Universal Gas Constant (J/kmol*K)
#define frdy 96485.0 // Faraday's Constant (C/mol)
#define temp 310.0 // Temperature (K)

/* Ion Concentrations */
#define nai 18.0 // Intracellular Na Concentration (mM)
#define nao 140.0 // Extracellular Na Concentration (mM)
#define cao 1.8 // Extracellular Ca Concentration (mM)
#define ki 145.0 // Intracellular K Concentration (mM)
#define ko 5.4 // Extracellular K Concentration (mM)

/* Fast Sodium Current (time dependant) */
#define gna 23.0 // Max. Conductance of the Na Channel (mS/uF)
#define ena (((R*temp) / frdy)*log(nao/nai)) // Reversal Potential of Na (mV)

/* Time-dependent potassium current*/
#define gk (0.282*sqrt(ko / 5.4)) // Channel Conductance of Rapidly Activating K Current (mS/uF)
#define ek (((R*temp) / frdy)*log(ko / ki)) // Reversal Potential of Rapidly Activating K Current (mV)

/* Potassium Current (time-independent) */
#define gk1 (0.6047*(sqrt(ko / 5.4))) // Channel Conductance of Time Independant K Current (mS/uF)
#define ek1 (((R*temp) / frdy)*log(ko / ki)) // Reversal Potential of Time Independant K Current (mV)

/* Plateau Potassium Current */
#define gkp 0.0183 // Channel Conductance of Plateau K Current (mS/uF)
#define ekp ek1 // Reversal Potential of Plateau K Current (mV)

/* Background Current */
//

void performance();

/* FVM functions */
double usource(double uu, double x, double y, double t);
int find_t(int array[], int size, int p);
double weight_cal(double xa1, double ya1, double xa2, double ya2, double xb1, double yb1, double xb2, double yb2);
double v_norm(double v[2]);
double v_multiply(double v1[2], double v2[2]);
void v_subtract(double v[2], double v1[2], double v2[2]);
void v_assign(double v[2], double v1, double v2);
double m_det(double tmp[2][2]);

void Allocate();
void Init();
void Send_to_Device();
void dV2();
void stimu(int stim_size);
void Update_V();
void FDM();
void Forward_Euler();
void Send_V();
void Free();
