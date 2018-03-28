#include "head.h"

double t = 0.0; // Time (ms)

/* Performance compared */
double Vmax = 0, V_left = 0, V_right = 0, left_peak, right_peak, conduction_t = 0;
double APD90; // Time of 90% Repolarization 
double Vold, v_onset;
int pos_L = 0; //the left side center index of the tissue h_area to test the AP information
int pos_R = 0;

int stimtime = (int)(1 / dt + 0.6);  //Time period during which stimulus is applied (ms)
int stim_size = 0;

extern double *h_V;
extern double *h_m;
extern double *h_h;
extern double *h_jj;
extern double *h_d;
extern double *h_f;
extern double *h_X;
extern double *h_cai;

/* FVM variable*/
extern int *h_stim_i;
extern int *h_adj_triangle; // save the index of  three neighbor triangles
extern double *h_weight; //h_weighting factor
extern double *h_D_CF1;
extern double *h_D_CF2;
extern double *h_D_CF3;
extern double *h_matrix_G; //for Least - Square method in Section 9.3
extern double *h_gC; //geometric interpolation factor related to the position of the element face f with respect to the nodes C and F.
extern double *h_area; // each triangle's area
extern double *h_Sf;//surface vector
extern double *h_unit_CF;
extern double *h_Tf; // component h_Tf is normal to h_Sf
extern double *h_aF;

int main(int argC, char* argv[])
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

	Allocate();
	Init();
	
	int boundary_marker[nt] = { 0 }; //whether a boundary triangle or not
	int numberOfneighbors[nt] = { 0 };
	double centroid[nt][2] = { 0 }; //get the centroid coordinates
	double gradientU_C[nt][2] = { 0 }; //the gradient at the centroid of an triangle with volume Vc
	double face_center[nt][3][2] = { 0 }; //get the center coordinate of triangle's three faces
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
		h_area[i] = fabs(0.5*(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x2*y1 - x3*y2));
		ap[i] = h_area[i] / dt;
		if (centroid[i][0]<=0.04){// apply the stimulus on the elements if x coordinate <= 0.2 cm
			h_stim_i[stim_size] = i;
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
					h_adj_triangle[i] = tt+1;
					diff_x = nodes[elements[i][0]][0] - nodes[elements[i][1]][0];
					diff_y = nodes[elements[i][0]][1] - nodes[elements[i][1]][1];
					face_center[i][0][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][1]][0]);
					face_center[i][0][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][1]][1]);
					v_assign(Sf_tmp, diff_y, -diff_x);
					h_Sf[i] = diff_y;//h_Sf[i][0][0]
					h_Sf[i+nt] = -diff_x;//h_Sf[i][0][1]
					v_assign(D_Cf, face_center[i][0][0] - centroid[i][0], face_center[i][0][1] - centroid[i][1]);
					h_D_CF1[i] = centroid[tt][0] - centroid[i][0];
					h_D_CF1[i+nt] = centroid[tt][1] - centroid[i][1];
					h_gC[i] = weight_cal(centroid[i][0], centroid[i][1], centroid[tt][0], centroid[tt][1],
						nodes[elements[i][0]][0], nodes[elements[i][0]][1], nodes[elements[i][1]][0], nodes[elements[i][1]][1]);//0.5;
					h_weight[i] = 1 / pow(sqrt(h_D_CF1[i] * h_D_CF1[i] + h_D_CF1[i+nt] * h_D_CF1[i+nt]), 3);
					if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 1){
						h_Sf[i] = -diff_y;
						h_Sf[i+nt] = diff_x;
					}
					flag1 = 1;
				}
				else if (exist1 && exist3){
					h_adj_triangle[i+nt] = tt+1;
					diff_x = nodes[elements[i][0]][0] - nodes[elements[i][2]][0];
					diff_y = nodes[elements[i][0]][1] - nodes[elements[i][2]][1];
					face_center[i][1][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][2]][0]);
					face_center[i][1][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][2]][1]);
					v_assign(Sf_tmp, diff_y, -diff_x);
					h_Sf[i+2*nt] = diff_y;//h_Sf[i][1][0]
					h_Sf[i+3*nt] = -diff_x;//h_Sf[i][1][1]
					v_assign(D_Cf, face_center[i][1][0] - centroid[i][0], face_center[i][1][1] - centroid[i][1]);
					h_D_CF2[i] = centroid[tt][0] - centroid[i][0];
					h_D_CF2[i+nt] = centroid[tt][1] - centroid[i][1];
					h_gC[i+nt] = weight_cal(centroid[i][0], centroid[i][1], centroid[tt][0], centroid[tt][1],
						nodes[elements[i][0]][0], nodes[elements[i][0]][1], nodes[elements[i][2]][0], nodes[elements[i][2]][1]);//0.5;
					h_weight[i+nt] = 1 / pow(sqrt(h_D_CF2[i] * h_D_CF2[i] + h_D_CF2[i+nt] * h_D_CF2[i+nt]), 3);
					if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
						h_Sf[i+2*nt] = -diff_y;
						h_Sf[i+3*nt] = diff_x;
					}
					flag2 = 1;
				}
				else if (exist2 && exist3){
					h_adj_triangle[i+2*nt] = tt+1;
					diff_x = nodes[elements[i][1]][0] - nodes[elements[i][2]][0];
					diff_y = nodes[elements[i][1]][1] - nodes[elements[i][2]][1];
					face_center[i][2][0] = 0.5*(nodes[elements[i][1]][0] + nodes[elements[i][2]][0]);
					face_center[i][2][1] = 0.5*(nodes[elements[i][1]][1] + nodes[elements[i][2]][1]);
					v_assign(Sf_tmp, diff_y, -diff_x);
					h_Sf[i+4*nt] = diff_y;
					h_Sf[i+5*nt] = -diff_x;
					v_assign(D_Cf, face_center[i][2][0] - centroid[i][0], face_center[i][2][1] - centroid[i][1]);
					h_D_CF3[i] = centroid[tt][0] - centroid[i][0];
					h_D_CF3[i+nt] = centroid[tt][1] - centroid[i][1];
					h_gC[i+2*nt] = weight_cal(centroid[i][0], centroid[i][1], centroid[tt][0], centroid[tt][1],
						nodes[elements[i][1]][0], nodes[elements[i][1]][1], nodes[elements[i][2]][0], nodes[elements[i][2]][1]);//0.5;
					h_weight[i+2*nt] = 1 / pow(sqrt(h_D_CF3[i] * h_D_CF3[i] + h_D_CF3[i+nt] * h_D_CF3[i+nt]), 3);
					if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
						h_Sf[i+4*nt] = -diff_y;
						h_Sf[i+5*nt] = diff_x;
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
				h_Sf[i] = diff_y;
				h_Sf[i+nt] = -diff_x;
				v_assign(D_Cf, face_center[i][0][0] - centroid[i][0], face_center[i][0][1] - centroid[i][1]);
				if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 1){
					h_Sf[i] = -diff_y;
					h_Sf[i+nt] = diff_x;
				}
			}
			if (!flag2){
				diff_x = nodes[elements[i][0]][0] - nodes[elements[i][2]][0];
				diff_y = nodes[elements[i][0]][1] - nodes[elements[i][2]][1];
				face_center[i][1][0] = 0.5*(nodes[elements[i][0]][0] + nodes[elements[i][2]][0]);
				face_center[i][1][1] = 0.5*(nodes[elements[i][0]][1] + nodes[elements[i][2]][1]);
				v_assign(Sf_tmp, diff_y, -diff_x);
				h_Sf[i+2*nt] = diff_y;
				h_Sf[i+3*nt] = -diff_x;
				v_assign(D_Cf, face_center[i][1][0] - centroid[i][0], face_center[i][1][1] - centroid[i][1]);
				if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
					h_Sf[i+2*nt] = -diff_y;
					h_Sf[i+3*nt] = diff_x;
				}
			}
			if (!flag3){
				diff_x = nodes[elements[i][1]][0] - nodes[elements[i][2]][0];
				diff_y = nodes[elements[i][1]][1] - nodes[elements[i][2]][1];
				face_center[i][2][0] = 0.5*(nodes[elements[i][1]][0] + nodes[elements[i][2]][0]);
				face_center[i][2][1] = 0.5*(nodes[elements[i][1]][1] + nodes[elements[i][2]][1]);
				v_assign(Sf_tmp, diff_y, -diff_x);
				h_Sf[i+4*nt] = diff_y;
				h_Sf[i+5*nt] = -diff_x;
				v_assign(D_Cf, face_center[i][2][0] - centroid[i][0], face_center[i][2][1] - centroid[i][1]);
				if ((D_Cf[0] * Sf_tmp[0] + D_Cf[1] * Sf_tmp[1]) < 0){
					h_Sf[i+4*nt] = -diff_y;
					h_Sf[i+5*nt] = diff_x;
				}
			}
		}
		//numberOfneighbors[i] = length(find(h_adj_triangle(i, :)~= 0)); //save the number of neighbor triangles
		//if (numberOfneighbors[i] <3)
		if (h_adj_triangle[i] == 0 || h_adj_triangle[i+nt] == 0 || h_adj_triangle[i+2*nt] == 0)
			boundary_marker[i] = 1; //whether a boundary triangle or not
	}

	double aC[nt] = { 0 };
	double bC[nt][3] = { 0 }; //non - orthogonal flux
	int Gauss[nt] = { 0 }; //triangles needed to use Green - Gauss method
	int count = 0;
	double Ef[2] = { 0 };
	double tmp[2][2] = { 0 };
	double Ef_value[3] = { 0 };
	double Distance_CF1, Distance_CF2, Distance_CF3;

	for (i = 0; i < nt; i++){
		Distance_CF1 = 0;
		Distance_CF2 = 0;
		Distance_CF3 = 0;
		memset(Ef_value, 0, sizeof(Ef_value));
		if (h_adj_triangle[i] != 0){ // if not boundary face
			h_unit_CF[i] = h_D_CF1[i] / sqrt(h_D_CF1[i]*h_D_CF1[i]+h_D_CF1[i+nt]*h_D_CF1[i+nt]);
			h_unit_CF[i+nt] = h_D_CF1[i+nt] / sqrt(h_D_CF1[i]*h_D_CF1[i]+h_D_CF1[i+nt]*h_D_CF1[i+nt]);//the unit vector e in the direction of CF
			v_assign(Sf_tmp, h_Sf[i], h_Sf[i+nt]);
			v_assign(Ef, (v_multiply(Sf_tmp, Sf_tmp) / (h_unit_CF[i] * Sf_tmp[0] + h_unit_CF[i+nt] * Sf_tmp[1]))*h_unit_CF[i],//over - relaxed approach in Section 8.6.4
				(v_multiply(Sf_tmp, Sf_tmp) / (h_unit_CF[i] * Sf_tmp[0] + h_unit_CF[i+nt] * Sf_tmp[1]))*h_unit_CF[i+nt]);
			Ef_value[0] = v_norm(Ef);
			h_Tf[i] = Sf_tmp[0] - Ef[0];//decomposition, h_Sf = h_Tf + Ef;
			h_Tf[i+nt] = Sf_tmp[1] - Ef[1];//decomposition, h_Sf = h_Tf + Ef;
			Distance_CF1 = sqrt(h_D_CF1[i]*h_D_CF1[i]+h_D_CF1[i+nt]*h_D_CF1[i+nt]);

			h_matrix_G[i] = h_matrix_G[i] + h_weight[i] * h_D_CF1[i] * h_D_CF1[i];//h_matrix_G[i][0][0]
			h_matrix_G[i+nt] = h_matrix_G[i+nt] + h_weight[i] * h_D_CF1[i] * h_D_CF1[i+nt];//h_matrix_G[i][0][1]
			h_matrix_G[i+2*nt] = h_matrix_G[i+2*nt] + h_weight[i] * h_D_CF1[i] * h_D_CF1[i+nt];//h_matrix_G[i][1][0]
			h_matrix_G[i+3*nt] = h_matrix_G[i+3*nt] + h_weight[i] * h_D_CF1[i+nt] * h_D_CF1[i+nt];//h_matrix_G[i][1][1]
		}
		if (h_adj_triangle[i+nt] != 0){ // if not boundary face
			h_unit_CF[i+2*nt] = h_D_CF2[i] / sqrt(h_D_CF2[i]*h_D_CF2[i]+h_D_CF2[i+nt]*h_D_CF2[i+nt]);
			h_unit_CF[i+3*nt] = h_D_CF2[i+nt] / sqrt(h_D_CF2[i]*h_D_CF2[i]+h_D_CF2[i+nt]*h_D_CF2[i+nt]);//the unit vector e in the direction of CF
			v_assign(Sf_tmp, h_Sf[i+2*nt], h_Sf[i+3*nt]);
			v_assign(Ef, (v_multiply(Sf_tmp, Sf_tmp) / (h_unit_CF[i+2*nt] * Sf_tmp[0] + h_unit_CF[i+3*nt] * Sf_tmp[1]))*h_unit_CF[i+2*nt],//over - relaxed approach in Section 8.6.4
				(v_multiply(Sf_tmp, Sf_tmp) / (h_unit_CF[i+2*nt] * Sf_tmp[0] + h_unit_CF[i+3*nt] * Sf_tmp[1]))*h_unit_CF[i+3*nt]);
			Ef_value[1] = v_norm(Ef);
			h_Tf[i+2*nt] = Sf_tmp[0] - Ef[0];//decomposition, h_Sf = h_Tf + Ef;
			h_Tf[i+3*nt] = Sf_tmp[1] - Ef[1];//decomposition, h_Sf = h_Tf + Ef;
			Distance_CF2 = sqrt(h_D_CF2[i]*h_D_CF2[i]+h_D_CF2[i+nt]*h_D_CF2[i+nt]);

			h_matrix_G[i] = h_matrix_G[i] + h_weight[i+nt] * h_D_CF2[i] * h_D_CF2[i];
			h_matrix_G[i+nt] = h_matrix_G[i+nt] + h_weight[i+nt] * h_D_CF2[i] * h_D_CF2[i+nt];
			h_matrix_G[i+2*nt] = h_matrix_G[i+2*nt] + h_weight[i+nt] * h_D_CF2[i] * h_D_CF2[i+nt];
			h_matrix_G[i+3*nt] = h_matrix_G[i+3*nt] + h_weight[i+nt] * h_D_CF2[i+nt] * h_D_CF2[i+nt];
		}
		if (h_adj_triangle[i+2*nt] != 0){ // if not boundary face
			h_unit_CF[i+4*nt] = h_D_CF3[i] / sqrt(h_D_CF3[i]*h_D_CF3[i]+h_D_CF3[i+nt]*h_D_CF3[i+nt]);
			h_unit_CF[i+5*nt] = h_D_CF3[i+nt] / sqrt(h_D_CF3[i]*h_D_CF3[i]+h_D_CF3[i+nt]*h_D_CF3[i+nt]);//the unit vector e in the direction of CF
			v_assign(Sf_tmp, h_Sf[i+4*nt], h_Sf[i+5*nt]);
			v_assign(Ef, (v_multiply(Sf_tmp, Sf_tmp) / (h_unit_CF[i+4*nt] * Sf_tmp[0] + h_unit_CF[i+5*nt] * Sf_tmp[1]))*h_unit_CF[i+4*nt],//over - relaxed approach in Section 8.6.4
				(v_multiply(Sf_tmp, Sf_tmp) / (h_unit_CF[i+4*nt] * Sf_tmp[0] + h_unit_CF[i+5*nt] * Sf_tmp[1]))*h_unit_CF[i+5*nt]);
			Ef_value[2] = v_norm(Ef);
			h_Tf[i+4*nt] = Sf_tmp[0] - Ef[0];//decomposition, h_Sf = h_Tf + Ef;
			h_Tf[i+5*nt] = Sf_tmp[1] - Ef[1];//decomposition, h_Sf = h_Tf + Ef;
			Distance_CF3 = sqrt(h_D_CF3[i]*h_D_CF3[i]+h_D_CF3[i+nt]*h_D_CF3[i+nt]);

			h_matrix_G[i] = h_matrix_G[i] + h_weight[i+2*nt] * h_D_CF3[i] * h_D_CF3[i];
			h_matrix_G[i+nt] = h_matrix_G[i+nt] + h_weight[i+2*nt] * h_D_CF3[i] * h_D_CF3[i+nt];
			h_matrix_G[i+2*nt] = h_matrix_G[i+2*nt] + h_weight[i+2*nt] * h_D_CF3[i] * h_D_CF3[i+nt];
			h_matrix_G[i+3*nt] = h_matrix_G[i+3*nt] + h_weight[i+2*nt] * h_D_CF3[i+nt] * h_D_CF3[i+nt];
		}

		//aC(i) = Ef_value(i, 1) / Distance_CF1 + Ef_value(i, 2) / Distance_CF2 + Ef_value(i, 3) / Distance_CF3;
		if (Distance_CF1 != 0)
			h_aF[i] = D * Ef_value[0] / Distance_CF1;
		if (Distance_CF2 != 0)
			h_aF[i+nt] = D *  Ef_value[1] / Distance_CF2;
		if (Distance_CF3 != 0)
			h_aF[i+2*nt] = D *  Ef_value[2] / Distance_CF3;
		aC[i] = h_aF[i] + h_aF[i+nt] + h_aF[i+2*nt];

		//decide whether singular matrix or not
		tmp[0][0] = h_matrix_G[i];
		tmp[0][1] = h_matrix_G[i+nt];
		tmp[1][0] = h_matrix_G[i+2*nt];
		tmp[1][1] = h_matrix_G[i+3*nt];
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

	double Uf[nt][3] = { 0 };//function value at the triangle face center
	for (i = 0; i < nt; i++){
		h_V[i] = -88.654973; // Initial Voltage (mv)
		h_m[i] = 0.000838;
		h_h[i] = 0.993336;
		h_jj[i] = 0.995484;
		h_d[i] = 0.000003;
		h_f[i] = 0.999745;
		h_X[i] = 0.000129;
		h_cai[i] = 0.0002; // Initial Intracellular Ca (mM)
	}

	int nstep = 4 / dt;; // snapshot interval 4 ms to save data files
	int index = 0;// filename index 
	char filename[100];

	//time recorder
	struct timeb start, end;
	int diff;
	ftime(&start);

	Send_to_Device();
	//	steps = (bcl*beats)/dt; // Time Loop Conditions
	int ncount;
	for (ncount = 1; ncount <= 600 / dt; ncount++){//simulation time is 600ms
		//-**** save data in file "ap"
		int fileflag = 0;
		if (ncount%nstep == 0){
			Send_V();
			for (i = 0; i < nt; i++){//save data every 4 ms
					if (fileflag == 0){
						sprintf(filename, "snapshot_AP%d", index);
						snapshot_AP = fopen(filename, "w");
						fileflag = 1;
						index++;
					}
					fprintf(snapshot_AP, "%g\t", h_V[i]);
					if (i%10==9){
						fprintf(snapshot_AP, "\n");
					}
			}
		}
		if (fileflag == 1){
			fclose(snapshot_AP);
		}
		fprintf(single_AP1, "%g\n", h_V[pos_L]);
		fprintf(single_AP2, "%g\n", h_V[pos_R]);
		performance();

		//-*********** step 1 to compute source term ****************
		
		dV2();
		//stimulation with a plane wave
		if (ncount >= 1 && ncount <= stimtime) { //stimulus is hold with 0.6 ms
			stimu(stim_size);
		}
		Update_V();
		
		//-*********** step 1 *******************

		//-*********** step 2 to compute diffusion term

		FDM();
		Forward_Euler();
		
		//-*********** step 2 *******
		t = t + dt;

	}

	ftime(&end);
	free();
	diff = (1000.0*(end.time - start.time) + (end.millitm - start.millitm));
	conduction_t = (right_peak - left_peak)*0.001; //condution time from left side to right side
	fprintf(fevaluation, "%d ms\n%g\n%g\n%g\n", diff, Vmax, APD90, length_x / conduction_t);
	fclose(fevaluation);
	fclose(single_AP1);
	fclose(single_AP2);
}

void Init(){
	int i;
	
	for(i=0;i<nt;i++){
		h_stim_i[i] = 0;
		h_area[i] = 0;
	}
	for(i=0;i<2*nt;i++){
		h_D_CF1[i] = 0;
		h_D_CF2[i] = 0;
		h_D_CF3[i] = 0;
	}
	for(i=0;i<3*nt;i++){	
		h_adj_triangle[i] = 0;
		h_weight[i] = 0;
		h_gC[i] = 0;
		h_aF[i] = 0;
	}
	for(i=0;i<2*2*nt;i++){	
		h_matrix_G[i] = 0;
	}
	for(i=0;i<2*3*nt;i++){	
		h_Sf[i] = 0;
		h_unit_CF[i] = 0;
		h_Tf[i] = 0;
	}
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
	if (h_V[pos_L] - V_left > 0){
		left_peak = t; // time reach peak at left side
		V_left = h_V[pos_L];
	}
	if (h_V[pos_R] - V_right > 0){
		right_peak = t; // time reach peak at right side
		V_right = h_V[pos_R];
	}
	if (h_V[pos_L]>Vmax)
		Vmax = h_V[pos_L];
	if (h_V[pos_L] >= (Vmax - 0.9*(Vmax - (-88.654973)))) //-88.654973 is the Resting Membrane Potential, also initial voltage
		APD90 = t; //  Time of 90% Repolarization 
}

