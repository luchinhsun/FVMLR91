/* 2D LR91 */
/* Operator splitting :
FVM for diffusion term
forward Euler method for ODE */
/* ref: An Advanced Algorithm for Solving Partial Differential Equation in Cardiac Conduction. 1999. */
/* Xiang Zhou, 2018/2/23 */

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/timeb.h>
#include <string.h>

//-*******mesh parameters for LR91 *******
double length_x = 1.0;//cm, the size of the tissue (cm)*(cm)
double length_y = 1.0;
int const np = 6529; //6529;//1665;//433;//117;//34; //number of Vertices
int const nt = 12800;//12800;//3200;//800;//200;//50; //number of Triangles

double D = 0.001;//D: diffusion coefficient cm^2/ms

/* Time Step */
double dt = 0.005; // Time step (ms)
double t = 0.0; // Time (ms)
int cutcount = 40 / dt;

/* Voltage */
double V[nt] = { 0 }; // Initial Voltage (mv)
double dV2[nt] = { 0 }; // second order derivatives of Voltage (mv)
double Vnew[nt] = { 0 };// New Voltage (mV)
double dvdt; // Change in Voltage / Change in Time (mV/ms)
double dvdtnew; // New dv/dt (mV/ms)

/* Total Current and Stimulus */
double st = -80.0; // Constant Stimulus (uA/cm^2)
double tstim; //Time Stimulus is Applied (ms)//Time to begin stimulus
int stimtime = (int)(1 / dt + 0.6);  //Time period during which stimulus is applied (ms)
double it[nt] = { 0 };; // Total current (uA/cm^2)
int stim_i[nt] = { 0 }; //{0,14,2,13,30,11,29,20,49,10,5, 3,48,6,9,8,27,16,50};// apply the stimulus on these elements 
int stim_size = 0;

/* Terms for Solution of Conductance and Reversal Potential */
const double R = 8314; // Universal Gas Constant (J/kmol*K)
const double frdy = 96485; // Faraday's Constant (C/mol)
double temp = 310; // Temperature (K)

/* Ion Concentrations */
double nai = 18;; // Intracellular Na Concentration (mM)
double nao = 140; // Extracellular Na Concentration (mM)
double cai[nt] = { 0 }; // Intracellular Ca Concentration (mM)
double cao = 1.8; // Extracellular Ca Concentration (mM)
double ki = 145; // Intracellular K Concentration (mM)
double ko = 5.4; // Extracellular K Concentration (mM)

/* Fast Sodium Current (time dependant) */
double ina[nt] = { 0 }; // Fast Na Current (uA/uF)
double gna = 23; // Max. Conductance of the Na Channel (mS/uF)
double ena = ((R*temp) / frdy)*log(nao / nai); // Reversal Potential of Na (mV)
//double am; // Na alpha-m rate constant (ms^-1)
//double bm; // Na beta-m rate constant (ms^-1)
//double ah; // Na alpha-h rate constant (ms^-1)
//double bh; // Na beta-h rate constant (ms^-1)
//double aj; // Na alpha-j rate constant (ms^-1)
//double bj; // Na beta-j rate constant (ms^-1)
//double mtau; // Na activation
//double htau; // Na inactivation
//double jtau; // Na inactivation
//double mss; // Na activation
//double hss; // Na inactivation
//double jss; // Na slow inactivation
double m[nt] = { 0 }; // Na activation
double h[nt] = { 0 }; // Na inactivation
double jj[nt] = { 0 }; // Na slow inactivation

/* Current through L-type Ca Channel */
double isi[nt] = { 0 }; // Slow inward current (uA/uF)
double esi[nt] = { 0 }; // Reversal Potential of si (mV)
//double dcai; // Change in myoplasmic Ca concentration (mM)
//double ad; // Ca alpha-d rate constant (ms^-1)
//double bd; // Ca beta-d rate constant (ms^-1)
//double af; // Ca alpha-f rate constant (ms^-1)
//double bf; // Ca beta-f rate constant (ms^-1)

double d[nt] = { 0 }; // Voltage dependant activation gate
double f[nt] = { 0 }; // Voltage dependant inactivation gate
double fca[nt] = { 0 }; // Ca dependant inactivation gate -from LR94
//double dss; // Steady-state value of activation gate d
//double taud; // Time constant of gate d (ms^-1)----mistake ???£¿ms£¿
//double fss; // Steady-state value of inactivation gate f
//double tauf; // Time constant of gate f (ms^-1)

/* Time-dependent potassium current*/
//double prnak = 0.01833;
//ek = ((R*temp) / frdy)*log((ko + prnak*nao) / (ki + prnak*nai));
double gk = 0.282*sqrt(ko / 5.4); // Channel Conductance of Rapidly Activating K Current (mS/uF)
double ek = ((R*temp) / frdy)*log(ko / ki); // Reversal Potential of Rapidly Activating K Current (mV)
double ik[nt] = { 0 }; // Rapidly Activating K Current (uA/uF)
double X[nt] = { 0 }; // Rapidly Activating K time-dependant activation  --gate X in LR91
//double ax; // K alpha-x rate constant (ms^-1)
//double bx; // K beta-x rate constant (ms^-1)
//double xss; // Steady-state value of inactivation gate xr  --gate X in LR91
//double taux; // Time constant of gate xr (ms^-1) --gate X in LR91
//double Xi; // K time-independent inactivation --gate Xi in LR91

/* Potassium Current (time-independent) */
double ik1[nt] = { 0 }; // Time-independent K current (uA/uF)
double gk1 = 0.6047*(sqrt(ko / 5.4)); // Channel Conductance of Time Independant K Current (mS/uF)
double ek1 = ((R*temp) / frdy)*log(ko / ki); // Reversal Potential of Time Independant K Current (mV)
//double ak1; // K alpha-ki rate constant (ms^-1)
//double bk1; // K beta-ki rate constant (ms^-1)
//double K1ss; // Steady-state value of K inactivation gate K1

/* Plateau Potassium Current */
double ikp[nt] = { 0 }; // Plateau K current (uA/uF)
double gkp = 0.0183; // Channel Conductance of Plateau K Current (mS/uF)
double ekp = ek1; // Reversal Potential of Plateau K Current (mV)
//double kp; // K plateau factor

/* Background Current */
double ib[nt] = { 0 }; // Background current (uA/uF)

/* Performance compared */
double Vmax, V_left = 0, V_right = 0, left_peak, right_peak, conduction_t = 0;
double APD90; // Time of 90% Repolarization 
double Vold, v_onset;
int pos_L = 0; //the left side center index of the tissue area to test the AP information
int pos_R = 0;
void performance();

/* Ion Current Functions */
void comp_ina(int i); // Calculates Fast Na Current
void comp_ical(int i); // Calculates Currents through L-Type Ca Channel
void comp_ik(int i); // Calculates Time-dependent K Current
void comp_ik1(int i); // Calculates Time-Independent K Current
void comp_ikp(int i); // Calculates Plateau K Current
void comp_ib(int i); // Calculates Background Current
void comp_it(int i); // Calculates Total Current

/* FVM functions */
double usource(double uu, double x, double y, double t);
int find_t(int array[], int size, int p);
double weight_cal(double xa1, double ya1, double xa2, double ya2, double xb1, double yb1, double xb2, double yb2);
double v_norm(double v[2]);
double v_multiply(double v1[2], double v2[2]);
void v_subtract(double v[2], double v1[2], double v2[2]);
void v_assign(double v[2], double v1, double v2);
double m_det(double tmp[2][2]);

int main(int argc, char* argv[])
{
	int i;
	/******* Mesh processing *******/
	double nodes[np][2] = { 0 };
	int elements[nt][3] = { 0 };
	FILE *fp_nodes;
	fp_nodes = fopen("nodes_c64.dat", "r");
	if (NULL == fp_nodes)
		return 2;
	for (i = 0; i < np; i++){
		fscanf(fp_nodes, "%lf %lf\n", &nodes[i][0], &nodes[i][1]);
		//nodes[i][1] = nodes[i][1] / 10;
		//nodes[i][2] = nodes[i][2] / 10;
	}
	fclose(fp_nodes);

	FILE *fp_elements;
	fp_elements = fopen("0base_elements_c64.dat", "r");
	if (NULL == fp_elements)
		return 2;
	for (i = 0; i < nt; i++)
		fscanf(fp_elements, "%d %d %d\n", &elements[i][0], &elements[i][1], &elements[i][2]);
	fclose(fp_elements);

	int boundary_marker[nt] = { 0 }; //whether a boundary triangle or not
	int numberOfneighbors[nt] = { 0 };
	int adj_triangle[nt][3] = { 0 }; // save the index of  three neighbor triangles
	double centroid[nt][2] = { 0 }; //get the centroid coordinates
	double area[nt] = { 0 }; // each triangle's area
	double Sf[nt][3][2] = { 0 };//surface vector
	double gC[nt][3] = { 0 }; //geometric interpolation factor related to the position of the element face f with respect to the nodes C and F.
	double weight[nt][3] = { 0 }; //weighting factor
	double gradientU_C[nt][2] = { 0 }; //the gradient at the centroid of an triangle with volume Vc
	double face_center[nt][3][2] = { 0 }; //get the center coordinate of triangle's three faces
	double D_CF1[nt][2] = { 0 }, D_CF2[nt][2] = { 0 }, D_CF3[nt][2] = { 0 };
	double ap[nt] = { 0 };
	double x1, x2, x3, y1, y2, y3;
	for (i = 0; i < nt; i++){
		x1 = nodes[elements[i][0]][0];
		x2 = nodes[elements[i][1]][0];
		x3 = nodes[elements[i][2]][0];
		y1 = nodes[elements[i][0]][1];
		y2 = nodes[elements[i][1]][1];
		y3 = nodes[elements[i][2]][1];
		v_assign(centroid[i], (1.0 / 3.0)*(x1 + x2 + x3), (1.0 / 3.0)*(y1 + y2 + y3));
		area[i] = fabs(0.5*(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x2*y1 - x3*y2));
		ap[i] = area[i] / dt;
		if (centroid[i][0]<=0.04){// apply the stimulus on the elements if x coordinate <= 0.2 cm
			stim_i[stim_size] = i;
			stim_size++;
		}
		if ((centroid[i][1] <0.51 && centroid[i][1] >0.49) && centroid[i][0] <0.01){
			pos_L = i;
		}
		if ((centroid[i][1] <0.51 && centroid[i][1] >0.49) && centroid[i][0] >0.99){
			pos_R = i;
		}
	}

	int flag1, flag2, flag3;
	int exist1, exist2, exist3;
	int p1, p2, p3;
	int tt;
	double diff_x, diff_y, Sf_tmp[2] = { 0 }, D_Cf[2] = { 0 };
	for (i = 0; i < nt; i++){
		//first, find the cell center's 3 adjacent triangles
		flag1 = 0;
		flag2 = 0;
		flag3 = 0;
		for (tt = 0; tt < nt; tt++){
			if (tt != i){
				exist1 = 0;
				exist2 = 0;
				exist3 = 0;
				p1 = find_t(elements[i], 3, elements[tt][0]);
				p2 = find_t(elements[i], 3, elements[tt][1]);
				p3 = find_t(elements[i], 3, elements[tt][2]);
				if (p1 != 0){
					if (p1 == 1)
						exist1 = 1;
					else if (p1 == 2)
						exist2 = 1;
					else if (p1 == 3)
						exist3 = 1;
				}
				if (p2 != 0){
					if (p2 == 1)
						exist1 = 1;
					else if (p2 == 2)
						exist2 = 1;
					else if (p2 == 3)
						exist3 = 1;
				}
				if (p3 != 0){
					if (p3 == 1)
						exist1 = 1;
					else if (p3 == 2)
						exist2 = 1;
					else if (p3 == 3)
						exist3 = 1;
				}
				if (!(exist1 && exist2 || exist1 && exist3 || exist2 && exist3)) // reduce computation time
					continue;
				if (exist1 && exist2){
					adj_triangle[i][0] = tt+1;
					diff_x = nodes[elements[i][0]][0] - nodes[elements[i][1]][0];
					diff_y = nodes[elements[i][0]][1] - nodes[elements[i][1]][1];
					face_center[i][0][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][1]][0]);
					face_center[i][0][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][1]][1]);
					v_assign(Sf_tmp, diff_y, -diff_x);
					Sf[i][0][0] = diff_y;
					Sf[i][0][1] = -diff_x;
					v_assign(D_Cf, face_center[i][0][0] - centroid[i][0], face_center[i][0][1] - centroid[i][1]);
					v_assign(D_CF1[i], centroid[tt][0] - centroid[i][0], centroid[tt][1] - centroid[i][1]);
					gC[i][0] = weight_cal(centroid[i][0], centroid[i][1], centroid[tt][0], centroid[tt][1],
						nodes[elements[i][0]][0], nodes[elements[i][0]][1], nodes[elements[i][1]][0], nodes[elements[i][1]][1]);//0.5;
					weight[i][0] = 1 / pow(sqrt(D_CF1[i][0] * D_CF1[i][0] + D_CF1[i][1] * D_CF1[i][1]), 3);
					if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 1){
						Sf[i][0][0] = -diff_y;
						Sf[i][0][1] = diff_x;
					}
					flag1 = 1;
				}
				else if (exist1 && exist3){
					adj_triangle[i][1] = tt+1;
					diff_x = nodes[elements[i][0]][0] - nodes[elements[i][2]][0];
					diff_y = nodes[elements[i][0]][1] - nodes[elements[i][2]][1];
					face_center[i][1][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][2]][0]);
					face_center[i][1][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][2]][1]);
					v_assign(Sf_tmp, diff_y, -diff_x);
					Sf[i][1][0] = diff_y;
					Sf[i][1][1] = -diff_x;
					v_assign(D_Cf, face_center[i][1][0] - centroid[i][0], face_center[i][1][1] - centroid[i][1]);
					v_assign(D_CF2[i], centroid[tt][0] - centroid[i][0], centroid[tt][1] - centroid[i][1]);
					gC[i][1] = weight_cal(centroid[i][0], centroid[i][1], centroid[tt][0], centroid[tt][1],
						nodes[elements[i][0]][0], nodes[elements[i][0]][1], nodes[elements[i][2]][0], nodes[elements[i][2]][1]);//0.5;
					weight[i][1] = 1 / pow(sqrt(D_CF2[i][0] * D_CF2[i][0] + D_CF2[i][1] * D_CF2[i][1]), 3);
					if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
						Sf[i][1][0] = -diff_y;
						Sf[i][1][1] = diff_x;
					}
					flag2 = 1;
				}
				else if (exist2 && exist3){
					adj_triangle[i][2] = tt+1;
					diff_x = nodes[elements[i][1]][0] - nodes[elements[i][2]][0];
					diff_y = nodes[elements[i][1]][1] - nodes[elements[i][2]][1];
					face_center[i][2][0] = 0.5*(nodes[elements[i][1]][0] + nodes[elements[i][2]][0]);
					face_center[i][2][1] = 0.5*(nodes[elements[i][1]][1] + nodes[elements[i][2]][1]);
					v_assign(Sf_tmp, diff_y, -diff_x);
					Sf[i][2][0] = diff_y;
					Sf[i][2][1] = -diff_x;
					v_assign(D_Cf, face_center[i][2][0] - centroid[i][0], face_center[i][2][1] - centroid[i][1]);
					v_assign(D_CF3[i], centroid[tt][0] - centroid[i][0], centroid[tt][1] - centroid[i][1]);
					gC[i][2] = weight_cal(centroid[i][0], centroid[i][1], centroid[tt][0], centroid[tt][1],
						nodes[elements[i][1]][0], nodes[elements[i][1]][1], nodes[elements[i][2]][0], nodes[elements[i][2]][1]);//0.5;
					weight[i][2] = 1 / pow(sqrt(D_CF3[i][0] * D_CF3[i][0] + D_CF3[i][1] * D_CF3[i][1]), 3);
					if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
						Sf[i][2][0] = -diff_y;
						Sf[i][2][1] = diff_x;
					}
					flag3 = 1;
				}
			}
		}
		
		if (!flag1 || !flag2 || !flag3){
			if (!flag1){
				diff_x = nodes[elements[i][0]][0] - nodes[elements[i][1]][0];
				diff_y = nodes[elements[i][0]][1] - nodes[elements[i][1]][1];
				face_center[i][0][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][1]][0]);
				face_center[i][0][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][1]][1]);
				v_assign(Sf_tmp, diff_y, -diff_x);
				Sf[i][0][0] = diff_y;
				Sf[i][0][1] = -diff_x;
				v_assign(D_Cf, face_center[i][0][0] - centroid[i][0], face_center[i][0][1] - centroid[i][1]);
				if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 1){
					Sf[i][0][0] = -diff_y;
					Sf[i][0][1] = diff_x;
				}
			}
			if (!flag2){
				diff_x = nodes[elements[i][0]][0] - nodes[elements[i][2]][0];
				diff_y = nodes[elements[i][0]][1] - nodes[elements[i][2]][1];
				face_center[i][1][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][2]][0]);
				face_center[i][1][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][2]][1]);
				v_assign(Sf_tmp, diff_y, -diff_x);
				Sf[i][1][0] = diff_y;
				Sf[i][1][1] = -diff_x;
				v_assign(D_Cf, face_center[i][1][0] - centroid[i][0], face_center[i][1][1] - centroid[i][1]);
				if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
					Sf[i][1][0] = -diff_y;
					Sf[i][1][1] = diff_x;
				}
			}
			if (!flag3){
				diff_x = nodes[elements[i][1]][0] - nodes[elements[i][2]][0];
				diff_y = nodes[elements[i][1]][1] - nodes[elements[i][2]][1];
				face_center[i][2][0] = 0.5*(nodes[elements[i][1]][0] + nodes[elements[i][2]][0]);
				face_center[i][2][1] = 0.5*(nodes[elements[i][1]][1] + nodes[elements[i][2]][1]);
				v_assign(Sf_tmp, diff_y, -diff_x);
				Sf[i][2][0] = diff_y;
				Sf[i][2][1] = -diff_x;
				v_assign(D_Cf, face_center[i][2][0] - centroid[i][0], face_center[i][2][1] - centroid[i][1]);
				if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
					Sf[i][2][0] = -diff_y;
					Sf[i][2][1] = diff_x;
				}
			}
		}
		//numberOfneighbors[i] = length(find(adj_triangle(i, :)~= 0)); //save the number of neighbor triangles
		//if (numberOfneighbors[i] <3)
		if (adj_triangle[i][0] == 0 || adj_triangle[i][1] == 0 || adj_triangle[i][2] == 0)
			boundary_marker[i] = 1; //whether a boundary triangle or not
	}

	double aC[nt] = { 0 };
	double aF[nt][3] = { 0 };
	double bC[nt][3] = { 0 }; //non - orthogonal flux
	double Tf[nt][3][2] = { 0 }; // component Tf is normal to Sf
	double matrix_G[nt][2][2] = { 0 }; //for Least - Square method in Section 9.3
	int Gauss[nt] = { 0 }; //triangles needed to use Green - Gauss method
	int count = 0;
	double unit_CF[nt][3][2] = { 0 }, Ef[2] = { 0 };
	double tmp[2][2] = { 0 };
	double Ef_value[3] = { 0 };
	double Distance_CF1, Distance_CF2, Distance_CF3;

	for (i = 0; i < nt; i++){
		Distance_CF1 = 0;
		Distance_CF2 = 0;
		Distance_CF3 = 0;
		memset(Ef_value, 0, sizeof(Ef_value));
		if (adj_triangle[i][0] != 0){ // if not boundary face
			v_assign(unit_CF[i][0], D_CF1[i][0] / v_norm(D_CF1[i]), D_CF1[i][1] / v_norm(D_CF1[i]));//the unit vector e in the direction of CF
			v_assign(Sf_tmp, Sf[i][0][0], Sf[i][0][1]);
			v_assign(Ef, (v_multiply(Sf_tmp, Sf_tmp) / v_multiply(unit_CF[i][0], Sf_tmp))*unit_CF[i][0][0],//over - relaxed approach in Section 8.6.4
				(v_multiply(Sf_tmp, Sf_tmp) / v_multiply(unit_CF[i][0], Sf_tmp))*unit_CF[i][0][1]);
			Ef_value[0] = v_norm(Ef);
			v_subtract(Tf[i][0], Sf_tmp, Ef);//decomposition, Sf = Tf + Ef;
			Distance_CF1 = v_norm(D_CF1[i]);

			matrix_G[i][0][0] = matrix_G[i][0][0] + weight[i][0] * D_CF1[i][0] * D_CF1[i][0];
			matrix_G[i][0][1] = matrix_G[i][0][1] + weight[i][0] * D_CF1[i][0] * D_CF1[i][1];
			matrix_G[i][1][0] = matrix_G[i][1][0] + weight[i][0] * D_CF1[i][0] * D_CF1[i][1];
			matrix_G[i][1][1] = matrix_G[i][1][1] + weight[i][0] * D_CF1[i][1] * D_CF1[i][1];
		}
		if (adj_triangle[i][1] != 0){ // if not boundary face
			v_assign(unit_CF[i][1], D_CF2[i][0] / v_norm(D_CF2[i]), D_CF2[i][1] / v_norm(D_CF2[i]));
			v_assign(Sf_tmp, Sf[i][1][0], Sf[i][1][1]);
			v_assign(Ef, (v_multiply(Sf_tmp, Sf_tmp) / v_multiply(unit_CF[i][1], Sf_tmp))*unit_CF[i][1][0],
				(v_multiply(Sf_tmp, Sf_tmp) / v_multiply(unit_CF[i][1], Sf_tmp))*unit_CF[i][1][1]);
			Ef_value[1] = v_norm(Ef);
			v_subtract(Tf[i][1], Sf_tmp, Ef);//decomposition, Sf = Tf + Ef;
			Distance_CF2 = v_norm(D_CF2[i]);

			matrix_G[i][0][0] = matrix_G[i][0][0] + weight[i][1] * D_CF2[i][0] * D_CF2[i][0];
			matrix_G[i][0][1] = matrix_G[i][0][1] + weight[i][1] * D_CF2[i][0] * D_CF2[i][1];
			matrix_G[i][1][0] = matrix_G[i][1][0] + weight[i][1] * D_CF2[i][0] * D_CF2[i][1];
			matrix_G[i][1][1] = matrix_G[i][1][1] + weight[i][1] * D_CF2[i][1] * D_CF2[i][1];
		}
		if (adj_triangle[i][2] != 0){ // if not boundary face
			v_assign(unit_CF[i][2], D_CF3[i][0] / v_norm(D_CF3[i]), D_CF3[i][1] / v_norm(D_CF3[i]));
			v_assign(Sf_tmp, Sf[i][2][0], Sf[i][2][1]);
			v_assign(Ef, (v_multiply(Sf_tmp, Sf_tmp) / v_multiply(unit_CF[i][2], Sf_tmp))*unit_CF[i][2][0],
				(v_multiply(Sf_tmp, Sf_tmp) / v_multiply(unit_CF[i][2], Sf_tmp))*unit_CF[i][2][1]);
			Ef_value[2] = v_norm(Ef);
			v_subtract(Tf[i][2], Sf_tmp, Ef);//decomposition, Sf = Tf + Ef;
			Distance_CF3 = v_norm(D_CF3[i]);

			matrix_G[i][0][0] = matrix_G[i][0][0] + weight[i][2] * D_CF3[i][0] * D_CF3[i][0];
			matrix_G[i][0][1] = matrix_G[i][0][1] + weight[i][2] * D_CF3[i][0] * D_CF3[i][1];
			matrix_G[i][1][0] = matrix_G[i][1][0] + weight[i][2] * D_CF3[i][0] * D_CF3[i][1];
			matrix_G[i][1][1] = matrix_G[i][1][1] + weight[i][2] * D_CF3[i][1] * D_CF3[i][1];
		}

		//aC(i) = Ef_value(i, 1) / Distance_CF1 + Ef_value(i, 2) / Distance_CF2 + Ef_value(i, 3) / Distance_CF3;
		if (Distance_CF1 != 0)
			aF[i][0] = D * Ef_value[0] / Distance_CF1;
		if (Distance_CF2 != 0)
			aF[i][1] = D *  Ef_value[1] / Distance_CF2;
		if (Distance_CF3 != 0)
			aF[i][2] = D *  Ef_value[2] / Distance_CF3;
		aC[i] = aF[i][0] + aF[i][1] + aF[i][2];

		//decide whether singular matrix or not
		tmp[0][0] = matrix_G[i][0][0];
		tmp[0][1] = matrix_G[i][0][1];
		tmp[1][0] = matrix_G[i][1][0];
		tmp[1][1] = matrix_G[i][1][1];
		if (fabs(m_det(tmp)) <= 1e-10){//singular matrix cannot be computed
			Gauss[i] = i+1;
		}
	}
	//-******* Mesh processing ******

	//-* Data File
	FILE *single_AP1;
	FILE *single_AP2;
	FILE *snapshot_AP;
	FILE *fevaluation;
	fevaluation = fopen("fevaluation", "w");
	single_AP1 = fopen("single_AP1", "w");
	single_AP2 = fopen("single_AP2", "w");

	//-* Beginning Ion Concentrations
	nai = 18; // Initial Intracellular Na (mM)
	nao = 140; // Initial Extracellular Na (mM)
	ki = 145; // Initial Intracellular K (mM)
	ko = 5.4; // Initial Extracellular K (mM)
	//cai = 0.0002; // Initial Intracellular Ca (mM)
	cao = 1.8; // Initial Extracellular Ca (mM)

	double Uf[nt][3] = { 0 };//function value at the triangle face center
	for (i = 0; i < nt; i++){
		V[i] = -88.654973; // Initial Voltage (mv)
		m[i] = 0.000838;
		h[i] = 0.993336;
		jj[i] = 0.995484;
		d[i] = 0.000003;
		f[i] = 0.999745;
		X[i] = 0.000129;
		cai[i] = 0.0002; // Initial Intracellular Ca (mM)
	}

	int nstep = 4 / dt;; // snapshot interval 4 ms to save data files
	int index = 0;// filename index 
	char filename[100];

	//time recorder
	struct timeb start, end;
	int diff;
	ftime(&start);

	//----------Big loop for time t--------------------------------------
	double vector_h[2] = { 0 }, gradient_f1[2] = { 0 }, gradient_f1_tmp[2] = { 0 },
		gradient_f2[2] = { 0 }, gradient_f2_tmp[2] = { 0 }, gradient_f3[2] = { 0 },
		gradient_f3_tmp[2] = { 0 };
	double sum_t;

	//	steps = (bcl*beats)/dt; // Time Loop Conditions 
	int ncount;
	for (ncount = 1; ncount <= 600 / dt; ncount++){//simulation time is 600ms
		//-**** save data in file "ap"
		int fileflag = 0;
		for (i = 0; i < nt; i++){
				if (ncount%nstep == 0){//save data every 4 ms
					if (fileflag == 0){
						sprintf(filename, "snapshot_AP%d", index);
						snapshot_AP = fopen(filename, "w");
						fileflag = 1;
						index++;
					}
					fprintf(snapshot_AP, "%g\t", V[i]);
					if (i%10==9){
						fprintf(snapshot_AP, "\n");
					}
				}
		}
		if (fileflag == 1){
			fclose(snapshot_AP);
		}
		fprintf(single_AP1, "%g\n", V[pos_L]);
		fprintf(single_AP2, "%g\n", V[pos_R]);
		performance();
		
		//-*********** step 1 to compute source term ****************
		for (i = 0; i < nt; i++){
			comp_ina(i);
			comp_ical(i);
			comp_ik(i);
			comp_ik1(i);
			comp_ikp(i);
			comp_ib(i);
			comp_it(i);

			dV2[i] = -it[i];
		}
		//stimulation with a plane wave
		if (ncount >= 1 && ncount <= stimtime) { //stimulus is hold with 0.6 ms
			for (i = 0; i < stim_size; i++){
				//stim_i[10 + 1] = { 14, 2, 13, 30, 11, 29, 20, 49, 10, 5 };
				dV2[stim_i[i]] = dV2[stim_i[i]] + (-st);
			}
		}
		for (i = 0; i < nt; i++){
			V[i] = V[i] + dt*dV2[i];
		}
		//-*********** step 1 *******************

		//-*********** step 2 to compute diffusion term
		//compute the gradient at the centroid of volume, using Least - Square method in Section 9.3
		for (i = 0; i < nt; i++){
			memset(vector_h, 0, sizeof(vector_h));
			if (adj_triangle[i][0] != 0){
				vector_h[0] = vector_h[0] + weight[i][0] * D_CF1[i][0] * (V[adj_triangle[i][0]-1] - V[i]);
				vector_h[1] = vector_h[1] + weight[i][0] * D_CF1[i][1] * (V[adj_triangle[i][0]-1] - V[i]);
			}
			if (adj_triangle[i][1] != 0){
				vector_h[0] = vector_h[0] + weight[i][1] * D_CF2[i][0] * (V[adj_triangle[i][1]-1] - V[i]);
				vector_h[1] = vector_h[1] + weight[i][1] * D_CF2[i][1] * (V[adj_triangle[i][1]-1] - V[i]);
			}
			if (adj_triangle[i][2] != 0){
				vector_h[0] = vector_h[0] + weight[i][2] * D_CF3[i][0] * (V[adj_triangle[i][2]-1] - V[i]);
				vector_h[1] = vector_h[1] + weight[i][2] * D_CF3[i][1] * (V[adj_triangle[i][2]-1] - V[i]);
			}
			tmp[0][0] = matrix_G[i][0][0];
			tmp[0][1] = matrix_G[i][0][1];
			tmp[1][0] = matrix_G[i][1][0];
			tmp[1][1] = matrix_G[i][1][1];
			if (fabs(m_det(tmp)) > 1e-10){//singular matrix cannot be computed
				//gradientU_C(i, :) = tmp\vector_h;
				if (tmp[0][0] != 0){// solve Ax=b
					gradientU_C[i][1] = (vector_h[1] - (tmp[1][0] * vector_h[0]) / tmp[0][0]) / (tmp[1][1] - (tmp[1][0] * tmp[0][1]) / tmp[0][0]);
					gradientU_C[i][0] = (vector_h[0] - tmp[0][1] * gradientU_C[i][1]) / tmp[0][0];
				}
				else{
					gradientU_C[i][1] = vector_h[0] / tmp[0][1];
					gradientU_C[i][0] = (vector_h[1] - tmp[1][1] * gradientU_C[i][1]) / tmp[1][0];
				}
				//double r1 = tmp[1][1] * gradientU_C[i][1] + tmp[1][2] * gradientU_C[i][2] - vector_h[1];
				//double r2 = tmp[2][1] * gradientU_C[i][1] + tmp[2][2] * gradientU_C[i][2] - vector_h[2];
			}
		}
		//if matrix tmp is singular, Green - Gauss method in Section 9.2 is used
		for (i = 0; i < nt; i++){
			if (Gauss[i]){
				if (adj_triangle[i][0] != 0)
					Uf[i][0] = (1 - gC[i][0])*V[adj_triangle[i][0]-1] + gC[i][0] * V[i];
				else
					Uf[i][0] = V[i];
				if (adj_triangle[i][1] != 0)
					Uf[i][1] = (1 - gC[i][1])*V[adj_triangle[i][1]-1] + gC[i][1] * V[i];
				else
					Uf[i][1] = V[i];
				if (adj_triangle[i][2] != 0)
					Uf[i][2] = (1 - gC[i][2])*V[adj_triangle[i][2]-1] + gC[i][2] * V[i];
				else
					Uf[i][2] = V[i];
				gradientU_C[i][0] = (1.0 / area[i])*(Uf[i][0] * Sf[i][0][0] + Uf[i][1] * Sf[i][1][0] + Uf[i][2] * Sf[i][2][0]);
				gradientU_C[i][1] = (1.0 / area[i])*(Uf[i][0] * Sf[i][0][1] + Uf[i][1] * Sf[i][1][1] + Uf[i][2] * Sf[i][2][1]);
				//aa = gradient_cal(t, centroid(i, 1), centroid(i, 2));
			}
		}
		//compute gradients at midpoins of three faces, since mid - point integration approximation is used.Section 9.4
		//gradient_fi: the interpolated gradient at the intersecting point of face i
		for (i = 0; i < nt; i++){
			if (adj_triangle[i][0] != 0){
				gradient_f1_tmp[0] = (1 - gC[i][0])*gradientU_C[adj_triangle[i][0]-1][0] + gC[i][0] * gradientU_C[i][0];
				gradient_f1_tmp[1] = (1 - gC[i][0])*gradientU_C[adj_triangle[i][0]-1][1] + gC[i][0] * gradientU_C[i][1];
				//Note: here a tmp variable "gradient_f1_tmp" is needed to avoid across assignment
				gradient_f1[0] = gradient_f1_tmp[0] + ((V[adj_triangle[i][0]-1] - V[i]) / v_norm(D_CF1[i]) - v_multiply(gradient_f1_tmp, unit_CF[i][0]))*unit_CF[i][0][0];
				gradient_f1[1] = gradient_f1_tmp[1] + ((V[adj_triangle[i][0]-1] - V[i]) / v_norm(D_CF1[i]) - v_multiply(gradient_f1_tmp, unit_CF[i][0]))*unit_CF[i][0][1];
			}
			else{
				gradient_f1[0] = 0; //boundary edge
				gradient_f1[1] = 0;
			}

			if (adj_triangle[i][1] != 0){
				gradient_f2_tmp[0] = (1 - gC[i][1])*gradientU_C[adj_triangle[i][1]-1][0] + gC[i][1] * gradientU_C[i][0];
				gradient_f2_tmp[1] = (1 - gC[i][1])*gradientU_C[adj_triangle[i][1]-1][1] + gC[i][1] * gradientU_C[i][1];
				gradient_f2[0] = gradient_f2_tmp[0] + ((V[adj_triangle[i][1]-1] - V[i]) / v_norm(D_CF2[i]) - v_multiply(gradient_f2_tmp, unit_CF[i][1]))*unit_CF[i][1][0];
				gradient_f2[1] = gradient_f2_tmp[1] + ((V[adj_triangle[i][1]-1] - V[i]) / v_norm(D_CF2[i]) - v_multiply(gradient_f2_tmp, unit_CF[i][1]))*unit_CF[i][1][1];
			}
			else{
				gradient_f2[0] = 0; //boundary edge
				gradient_f2[1] = 0;
			}

			if (adj_triangle[i][2] != 0){
				gradient_f3_tmp[0] = (1 - gC[i][2])*gradientU_C[adj_triangle[i][2]-1][0] + gC[i][2] * gradientU_C[i][0];
				gradient_f3_tmp[1] = (1 - gC[i][2])*gradientU_C[adj_triangle[i][2]-1][1] + gC[i][2] * gradientU_C[i][1];
				gradient_f3[0] = gradient_f3_tmp[0] + ((V[adj_triangle[i][2]-1] - V[i]) / v_norm(D_CF3[i]) - v_multiply(gradient_f3_tmp, unit_CF[i][2]))*unit_CF[i][2][0];
				gradient_f3[1] = gradient_f3_tmp[1] + ((V[adj_triangle[i][2]-1] - V[i]) / v_norm(D_CF3[i]) - v_multiply(gradient_f3_tmp, unit_CF[i][2]))*unit_CF[i][2][1];
			}
			else{
				gradient_f3[0] = 0; //boundary edge
				gradient_f3[1] = 0;
			}
			bC[i][0] = D * v_multiply(gradient_f1, Tf[i][0]); //non-orthogonal flux
			bC[i][1] = D * v_multiply(gradient_f2, Tf[i][1]);
			bC[i][2] = D * v_multiply(gradient_f3, Tf[i][2]);
		}
		//flux = 0;
		for (i = 0; i < nt; i++){
			sum_t = 0;
			//double aa1, aa2, aa3;
			if (adj_triangle[i][0] != 0){
				sum_t = sum_t + aF[i][0] * V[adj_triangle[i][0]-1];
			//	aa1 = aF[i][1] * (V[adj_triangle[i][1]] - V[i]);
			}
			if (adj_triangle[i][1] != 0){
				sum_t = sum_t + aF[i][1] * V[adj_triangle[i][1]-1];
			//	aa2 = aF[i][2] * (V[adj_triangle[i][2]] - V[i]);
			}
			if (adj_triangle[i][2] != 0){
				sum_t = sum_t + aF[i][2] * V[adj_triangle[i][2]-1];
			//	aa3 = aF[i][3] * (V[adj_triangle[i][3]] - V[i]);
			}
			//Vnew[i] = (sum_t + (ap[i] - (aF[i][1] + aF[i][2] + aF[i][3]))*V[i] + bC[i][1] + bC[i][2] + bC[i][3]) / ap[i];
			Vnew[i] = V[i] + (sum_t + (-(aF[i][0] + aF[i][1] + aF[i][2]))*V[i] + bC[i][0] + bC[i][1] + bC[i][2]) / ap[i];
			//Unew(i) = (-aC(i)*V(i) + (aF(i, 1)*V(adj_triangle(i, 1)) ...
			// +aF(i, 2)*V(adj_triangle(i, 2)) + aF(i, 3)*V(adj_triangle(i, 3))) + sum(bC(i, :))) / ap(i) + V(i);
		}
		for (i = 0; i < nt; i++){
			V[i] = Vnew[i];
		}
		//-*********** step 2 *******
		
		t = t + dt;

		//-***********trancation 1/2 of the plane wave to generate a spiral wave******
		//if (ncount == cutcount){
		//	for (i = 1; i < nx / 2; i++){
		//		for (j = 1; j < ny; j++){
		//			V[i] = -88.654973; // Initial Voltage (mv)
		//			m[i] = 0.000838;
		//			h[i] = 0.993336;
		//			jj[i] = 0.995484;
		//			d[i] = 0.000003;
		//			f[i] = 0.999745;
		//			X[i] = 0.000129;
		//			cai[i] = 0.0002; // Initial Intracellular Ca (mM)
		//		}
		//	}
		//}
	}

	ftime(&end);
	diff = (1000.0*(end.time - start.time) + (end.millitm - start.millitm));
	conduction_t = (right_peak - left_peak)*0.001; //condution time from left side to right side
	fprintf(fevaluation, "%d ms\n%g\n%g\n%g\n", diff, Vmax, APD90, length_x / conduction_t);
	fclose(fevaluation);
	fclose(single_AP1);
	fclose(single_AP2);
}

//find the position of p in array[]
int find_t(int array[], int size, int p){
	int i;
	for (i = 0; i < size; i++){
		if (p == array[i]){
			return i+1;
		}
	}
	return 0;
}

double weight_cal(double xa1, double ya1, double xa2, double ya2, double xb1, double yb1, double xb2, double yb2){
	//compute the intersecting point of two lines
	double k1, b1, k2, b2, xc, yc;
	double Ff[2], CF[2], gC;
	if (xa1 - xa2 == 0){
		k2 = (yb1 - yb2) / (xb1 - xb2);
		b2 = (xb1*yb2 - xb2*yb1) / (xb1 - xb2);
		xc = xa1;
		yc = k2*xc + b2;
	}
	else if (xb1 - xb2 == 0){
		k1 = (ya1 - ya2) / (xa1 - xa2);
		b1 = (xa1*ya2 - xa2*ya1) / (xa1 - xa2);
		xc = xb1;
		yc = k1*xc + b1;
	}
	else{
		k1 = (ya1 - ya2) / (xa1 - xa2);
		b1 = (xa1*ya2 - xa2*ya1) / (xa1 - xa2);
		k2 = (yb1 - yb2) / (xb1 - xb2);
		b2 = (xb1*yb2 - xb2*yb1) / (xb1 - xb2);
		xc = (b2 - b1) / (k1 - k2);
		yc = k1*xc + b1;
	}
	Ff[0] = xa2 - xc;
	Ff[1] = ya2 - yc;
	CF[0] = xa2 - xa1;
	CF[1] = ya2 - ya1;
	gC = sqrt(Ff[0] * Ff[0] + Ff[1] * Ff[1]) / sqrt(CF[0] * CF[0] + CF[1] * CF[1]);
	return gC;
}

//Calculate the 2-norm of a vector
double v_norm(double v[2]){
	return sqrt(v[0] * v[0] + v[1] * v[1]);
}

//vector multiply 
double v_multiply(double v1[2], double v2[2]){
	return v1[0] * v2[0] + v1[1] * v2[1];
}

//vector subtract
void v_subtract(double v[2], double v1[2], double v2[2]){
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
}

//vector assignment
void v_assign(double v[2], double v1, double v2){
	v[0] = v1;
	v[1] = v2;
}

//matrix det()
double m_det(double tmp[2][2]){
	return tmp[0][0] * tmp[1][1] - tmp[0][1] * tmp[1][0];
}

//performance 
void performance(){
	if (V[pos_L] - V_left > 0){
		left_peak = t; // time reach peak at left side
		V_left = V[pos_L];
	}
	if (V[pos_R] - V_right > 0){
		right_peak = t; // time reach peak at right side
		V_right = V[pos_R];
	}
	if (V[pos_L]>Vmax)
		Vmax = V[pos_L];
	if (V[pos_L] >= (Vmax - 0.9*(Vmax - (-88.654973)))) //-88.654973 is the Resting Membrane Potential, also initial voltage
		APD90 = t; //  Time of 90% Repolarization 
}
//-********************************************************/
//-* Functions that describe the currents begin here */

//Fast sodium current
void comp_ina(int i) {
	/*gate variables can not be shared, should be local due to data racing !!!!!!!!*/
	double am = 0.32*(V[i] + 47.13) / (1 - exp(-0.1*(V[i] + 47.13)));
	double bm = 0.08*exp(-V[i] / 11);
	double ah, bh, aj, bj;
	if (V[i] < -40) {
		ah = 0.135*exp((80 + V[i]) / -6.8);
		bh = 3.56*exp(0.079*V[i]) + 310000 * exp(0.35*V[i]);
		aj = (-127140 * exp(0.2444*V[i]) - 0.00003474*exp(-0.04391*V[i]))*((V[i] + 37.78) / (1 + exp(0.311*(V[i] + 79.23))));
		bj = (0.1212*exp(-0.01052*V[i])) / (1 + exp(-0.1378*(V[i] + 40.14)));
	}
	else {
		ah = 0;
		bh = 1 / (0.13*(1 + exp((V[i] + 10.66) / -11.1)));
		aj = 0;
		bj = (0.3*exp(-0.0000002535*V[i])) / (1 + exp(-0.1*(V[i] + 32)));
	}
	double mtau = 1 / (am + bm);
	double htau = 1 / (ah + bh);
	double jtau = 1 / (aj + bj);

	double mss = am*mtau;
	double hss = ah*htau;
	double jss = aj*jtau;

	m[i] = mss - (mss - m[i])*exp(-dt / mtau);
	h[i] = hss - (hss - h[i])*exp(-dt / htau);
	jj[i] = jss - (jss - jj[i])*exp(-dt / jtau);

	ina[i] = gna*m[i] * m[i] * m[i] * h[i] * jj[i] * (V[i] - ena);
}

//Slow inward current
void comp_ical(int i) {
	esi[i] = 7.7 - 13.0287*log(cai[i]);

	double ad = 0.095*exp(-0.01*(V[i] - 5)) / (1 + exp(-0.072*(V[i] - 5)));
	double bd = 0.07*exp(-0.017*(V[i] + 44)) / (1 + exp(0.05*(V[i] + 44)));
	double af = 0.012*exp(-0.008*(V[i] + 28)) / (1 + exp(0.15*(V[i] + 28)));
	double bf = 0.0065*exp(-0.02*(V[i] + 30)) / (1 + exp(-0.2*(V[i] + 30)));

	double taud = 1 / (ad + bd);
	double tauf = 1 / (af + bf);

	double dss = ad*taud;
	double fss = af*tauf;

	d[i] = dss - (dss - d[i])*exp(-dt / taud);
	f[i] = fss - (fss - f[i])*exp(-dt / tauf);

	isi[i] = 0.09*d[i] * f[i] * (V[i] - esi[i]);

	double dcai = -0.0001*isi[i] + 0.07*(0.0001 - cai[i]);

	cai[i] = cai[i] + dcai*dt;
}

//Time-dependent potassium current
void comp_ik(int i) {
	double ax = 0.0005*exp(0.083*(V[i] + 50)) / (1 + exp(0.057*(V[i] + 50)));
	double bx = 0.0013*exp(-0.06*(V[i] + 20)) / (1 + exp(-0.04*(V[i] + 20)));

	double taux = 1 / (ax + bx);
	double xss = ax*taux;
	X[i] = xss - (xss - X[i])*exp(-dt / taux);

	double Xi;
	if (V[i] > -100) {
		Xi = 2.837*(exp(0.04*(V[i] + 77)) - 1) / ((V[i] + 77)*exp(0.04*(V[i] + 35)));
	}
	else {
		Xi = 1;
	}

	ik[i] = gk*X[i] * Xi*(V[i] - ek);
}


//Time-independent potassium current
void comp_ik1(int i) {
	double ak1 = 1.02 / (1 + exp(0.2385*(V[i] - ek1 - 59.215)));
	double bk1 = (0.49124*exp(0.08032*(V[i] - ek1 + 5.476)) + exp(0.06175*(V[i] - ek1 - 594.31))) / (1 + exp(-0.5143*(V[i] - ek1 + 4.753)));
	double K1ss = ak1 / (ak1 + bk1);
	ik1[i] = gk1*K1ss*(V[i] - ek1);
}


//Plateau potassium current
void comp_ikp(int i) {
	double kp = 1 / (1 + exp((7.488 - V[i]) / 5.98));

	ikp[i] = gkp*kp*(V[i] - ekp);
}

//Background current
void comp_ib(int i) {
	ib[i] = 0.03921*(V[i] + 59.87);
}

//-* Total sum of currents is calculated here, if the time is between stimtime = 0 and stimtime = 0.5 (ms), a stimulus is applied */
void comp_it(int i) {
	//	if (t >= 5 && t<(5 + 0.5)) {
	//		it[i] = st + ina[i] + isi[i] + ik[i] + ik1[i] + ikp[i] + ib[i];
	//	}else {
	it[i] = ina[i] + isi[i] + ik[i] + ik1[i] + ikp[i] + ib[i];
	//	}
}


//-* Values are printed to a file called ap. The voltage and currents can be plotted versus time using graphing software. */
//void prttofile() {
//	if (t>(0) && t<(bcl*beats))
//	{
//		fprintf(ap, "%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
//			t, v, nai, ki, cai, ina, isi, ikr, iki, ikp, ib);
//		//printf("%.5f\t%g\n", t, v);
//		//printf("%.3f\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
//		//	t, v, nai, ki, cai, ina, isi, ikr, iki, ikp, ib);
//	}
//	//nai, ki, cai are the Intracellular Concentration of nai, ki, cai
//}
