#include "head.h"

#define tpb 256

extern double *d_it;
extern double *d_V;
extern double *d_dV2;
extern double *d_Vnew;
extern double *d_m;
extern double *d_h;
extern double *d_jj;
extern double *d_d;
extern double *d_f;
extern double *d_X;
extern double *d_cai;

extern double *dcai;

extern int *d_stim_i;
extern int *d_adj_triangle;
extern double *d_weight;
extern double *d_D_CF1;
extern double *d_D_CF2;
extern double *d_D_CF3;
extern double *d_matrix_G;
extern double *d_gC;
extern double *d_area;
extern double *d_Sf;
extern double *d_unit_CF;
extern double *d_Tf;
extern double *d_aF;

extern double *d_gradientU_C;

__device__ void comp_it(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *dcai, double *d_X, double *d_it,  int I, int k) {
	d_it[k] = 0.0;

	//comp_ina
        double am = 0.32*(d_V[k] + 47.13) / (1 - exp(-0.1*(d_V[k] + 47.13)));
        double bm = 0.08*exp(-d_V[k] / 11);
	double ah, bh, aj ,bj;
        if (d_V[k] < -40.0) {
                ah = 0.135*exp((80 + d_V[k]) / -6.8);
                bh = 3.56*exp(0.079*d_V[k]) + 310000 * exp(0.35*d_V[k]);
                aj = (-127140 * exp(0.2444*d_V[k]) - 0.00003474*exp(-0.04391*d_V[k]))*((d_V[k] + 37.78)/(1 + exp(0.311*(d_V[k] + 79.23))));
                bj = (0.1212*exp(-0.01052*d_V[k])) / (1 + exp(-0.1378*(d_V[k] + 40.14)));
        }
        else {
                ah = 0;
                bh = 1 / (0.13*(1 + exp((d_V[k] + 10.66) / -11.1)));
                aj = 0;
                bj = (0.3*exp(-0.0000002535*d_V[k])) / (1 + exp(-0.1*(d_V[k] + 32)));
        }
        double mtau = 1 / (am + bm);
        double htau = 1 / (ah + bh);
	double jtau = 1 / (aj + bj);

        double mss = am*mtau;
        double hss = ah*htau;
        double jss = aj*jtau;

        d_m[k] = mss - (mss - d_m[k])*exp(-dt / mtau);
        d_h[k] = hss - (hss - d_h[k])*exp(-dt / htau);
        d_jj[k] = jss - (jss - d_jj[k])*exp(-dt / jtau);
		
        d_it[k] += gna*d_m[k] * d_m[k] * d_m[k] * d_h[k] * d_jj[k] * (d_V[k] - ena);
	//comp_ical
	__shared__ double esi[tpb];
	__shared__ double isi[tpb];
        esi[I] = 7.7 - 13.0287*log(d_cai[k]);

        double ad = 0.095*exp(-0.01*(d_V[k] - 5)) / (1 + exp(-0.072*(d_V[k] - 5)));
        double bd = 0.07*exp(-0.017*(d_V[k] + 44)) / (1 + exp(0.05*(d_V[k] + 44)));
        double af = 0.012*exp(-0.008*(d_V[k] + 28)) / (1 + exp(0.15*(d_V[k] + 28)));
        double bf = 0.0065*exp(-0.02*(d_V[k] + 30)) / (1 + exp(-0.2*(d_V[k] + 30)));

        double taud = 1 / (ad + bd);
        double tauf = 1 / (af + bf);

        double dss = ad*taud;
        double fss = af*tauf;

        d_d[k] = dss - (dss - d_d[k])*exp(-dt / taud);
        d_f[k] = fss - (fss - d_f[k])*exp(-dt / tauf);

        isi[I] = 0.09*d_d[k] * d_f[k] * (d_V[k] - esi[I]);

        dcai[k] = -0.0001*isi[I] + 0.07*(0.0001 - d_cai[k]);

        d_cai[k] = d_cai[k] + dcai[k]*dt;
	d_it[k] = d_it[k] + isi[I];
	//comp_ik
        double ax = 0.0005*exp(0.083*(d_V[k] + 50)) / (1 + exp(0.057*(d_V[k] + 50)));
        double bx = 0.0013*exp(-0.06*(d_V[k] + 20)) / (1 + exp(-0.04*(d_V[k] + 20)));

        double taux = 1 / (ax + bx);
        double xss = ax*taux;
        d_X[k] = xss - (xss - d_X[k])*exp(-dt / taux);

	double Xi;
        if (d_V[k] > -100) {
                Xi = 2.837*(exp(0.04*(d_V[k] + 77)) - 1)/((d_V[k] + 77)*exp(0.04*(d_V[k] + 35)));
        }
        else {
                Xi = 1;
        }
        d_it[k] += gk*d_X[k] * Xi*(d_V[k] - ek);
	//comp_ik1
        double ak1 = 1.02 / (1 + exp(0.2385*(d_V[k] - ek1 - 59.215)));
        double bk1 = (0.49124*exp(0.08032*(d_V[k] - ek1 + 5.476)) + exp(0.06175*(d_V[k] - ek1 - 594.31))) / (1 + exp(-0.5143*(d_V[k] - ek1 + 4.753)));
        double K1ss = ak1 / (ak1 + bk1);

        d_it[k] += gk1*K1ss*(d_V[k] - ek1);
	//comp_ikp
        double kp = 1 / (1 + exp((7.488 - d_V[k]) / 5.98));

        d_it[k] += gkp*kp*(d_V[k] - ekp);
	//comp_ib
        d_it[k] += 0.03921*(d_V[k] + 59.87);

}

__global__ void comp_dV2(double *d_V, double *d_m, double *d_h, double *d_jj, double *d_d, double *d_f, double *d_cai, double *dcai, double *d_X, double *d_it, double *d_dV2){

	int k = threadIdx.x + blockIdx.x * blockDim.x;
    int I = threadIdx.x;

    if(k<nt){
		comp_it(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, I, k);
		d_dV2[k] = -d_it[k];
	}
}

void dV2(){
	int bpg;

        bpg = (nt+tpb-1)/tpb;
        comp_dV2<<<bpg, tpb>>>(d_V, d_m, d_h, d_jj, d_d, d_f, d_cai, dcai, d_X, d_it, d_dV2);
}

__global__ void plane_waves(double *d_dV2, int *d_stim_i, int stim_size){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<stim_size){
		d_dV2[d_stim_i[k]] = d_dV2[d_stim_i[k]] - st ;
	}
}

void stimu(int stim_size){
	int bpg;
	bpg = (stim_size+tpb-1)/tpb;
	plane_waves<<<bpg, tpb>>>(d_dV2, d_stim_i, stim_size);
	//cudaDeviceSynchronize();
}

__global__ void CUDAUpdate_V(double *d_V, double *d_dV2){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nt){
		d_V[k] = d_V[k] + dt*d_dV2[k];
	}
}

void Update_V(){
	int bpg;
	bpg = (nt+tpb-1)/tpb;
	CUDAUpdate_V<<<bpg, tpb>>>(d_V, d_dV2);
	//cudaDeviceSynchronize();
}

__global__ void diffusion(double *d_V, int *d_adj_triangle, double *d_weight, double *d_D_CF1, double *d_D_CF2, double *d_D_CF3, double *d_matrix_G, double *d_gC, double *d_area, double *d_Sf, double *d_gradientU_C){
	int k = threadIdx.x + blockIdx.x * blockDim.x;
    int I = threadIdx.x;
	
	__shared__ double vector_h[tpb][2];
	//__shared__ double vector_h[tpb];
	//__shared__ double vector_h1[tpb];
	__shared__ double tmp[tpb][2][2];
	__shared__ double Uf[tpb][3];
	
	if(k<nt){
		//compute the gradient at the centroid of volume, using Least - Square method in Section 9.3
		vector_h[I][0] = 0.0;
		vector_h[I][1] = 0.0;
		vector_h[I][0] = vector_h[I][0]	+ (d_weight[k] * d_D_CF1[k] * (d_V[d_adj_triangle[k]-1] - d_V[k]))*(d_adj_triangle[k] != 0) + (d_weight[k+nt] * d_D_CF2[k] * (d_V[d_adj_triangle[k+nt]-1] - d_V[k]))*(d_adj_triangle[k+nt] != 0) + (d_weight[k+2*nt] * d_D_CF3[k] * (d_V[d_adj_triangle[k+2*nt]-1] - d_V[k]))*(d_adj_triangle[k+2*nt] != 0);
		vector_h[I][1] = vector_h[I][1] + (d_weight[k] * d_D_CF1[k+nt] * (d_V[d_adj_triangle[k]-1] - d_V[k]))*(d_adj_triangle[k] != 0) + (d_weight[k+nt] * d_D_CF2[k+nt] * (d_V[d_adj_triangle[k+nt]-1] - d_V[k]))*(d_adj_triangle[k+nt] != 0) + (d_weight[k+2*nt] * d_D_CF3[k+nt] * (d_V[d_adj_triangle[k+2*nt]-1] - d_V[k]))*(d_adj_triangle[k+2*nt] != 0);
		//vector_h[I] = 0.0;
		//vector_h1[I] = 0.0;
		//vector_h[I] = vector_h[I] + (d_weight[k] * d_D_CF1[k] * (d_V[d_adj_triangle[k]-1] - d_V[k]))*(d_adj_triangle[k] != 0) + (d_weight[k+nt] * d_D_CF2[k] * (d_V[d_adj_triangle[k+nt]-1] - d_V[k]))*(d_adj_triangle[k+nt] != 0) + (d_weight[k+2*nt] * d_D_CF3[k] * (d_V[d_adj_triangle[k+2*nt]-1] - d_V[k]))*(d_adj_triangle[k+2*nt] != 0);
		//vector_h1[I] = vector_h1[I] + (d_weight[k] * d_D_CF1[k+nt] * (d_V[d_adj_triangle[k]-1] - d_V[k]))*(d_adj_triangle[k] != 0) + (d_weight[k+nt] * d_D_CF2[k+nt] * (d_V[d_adj_triangle[k+nt]-1] - d_V[k]))*(d_adj_triangle[k+nt] != 0) + (d_weight[k+2*nt] * d_D_CF3[k+nt] * (d_V[d_adj_triangle[k+2*nt]-1] - d_V[k]))*(d_adj_triangle[k+2*nt] != 0);
		
		tmp[I][0][0] = d_matrix_G[k];
		tmp[I][0][1] = d_matrix_G[k+nt];
		tmp[I][1][0] = d_matrix_G[k+2*nt];
		tmp[I][1][1] = d_matrix_G[k+3*nt];
		if(fabs(tmp[I][0][0]*tmp[I][1][1] - tmp[I][1][0]*tmp[I][0][1])>1e-10){
			if (tmp[I][0][0] != 0){// solve Ax=b
					d_gradientU_C[k+nt] = (vector_h[I][1] - (tmp[I][1][0] * vector_h[I][0]) / tmp[I][0][0]) / (tmp[I][1][1] - (tmp[I][1][0] * tmp[I][0][1]) / tmp[I][0][0]);
					d_gradientU_C[k] = (vector_h[I][0] - tmp[I][0][1] * d_gradientU_C[k+nt]) / tmp[I][0][0];
			}else{
					d_gradientU_C[k+nt] = vector_h[I][0] / tmp[I][0][1];
					d_gradientU_C[k] = (vector_h[I][1] - tmp[I][1][1] * d_gradientU_C[k+nt]) / tmp[I][1][0];
			}
		}
		
		//if matrix tmp is singular, Green - Gauss method in Section 9.2 is used
		if (fabs(tmp[I][0][0] * tmp[I][1][1] - tmp[I][0][1] * tmp[I][1][0])<= 1e-10){
			Uf[I][0] = ((1 - d_gC[k])*d_V[d_adj_triangle[k]-1] + d_gC[k] * d_V[k])*(d_adj_triangle[k] != 0) + d_V[k]*(d_adj_triangle[k] == 0);
			Uf[I][1] = ((1 - d_gC[k+nt])*d_V[d_adj_triangle[k+nt]-1] + d_gC[k+nt] * d_V[k])*(d_adj_triangle[k+nt] != 0) + d_V[k]*(d_adj_triangle[k+nt] == 0);
			Uf[I][2] = ((1 - d_gC[k+2*nt])*d_V[d_adj_triangle[k+2*nt]-1] + d_gC[k+2*nt] * d_V[k])*(d_adj_triangle[k+2*nt] != 0) + d_V[k]*(d_adj_triangle[k+2*nt] == 0);
			d_gradientU_C[k] = (1.0 / d_area[k])*(Uf[I][0] * d_Sf[k] + Uf[I][1] * d_Sf[k+2*nt] + Uf[I][2] * d_Sf[k+4*nt]);
			d_gradientU_C[k+nt] = (1.0 / d_area[k])*(Uf[I][0] * d_Sf[k+nt] + Uf[I][1] * d_Sf[k+3*nt] + Uf[I][2] * d_Sf[k+5*nt]);
		}
	}
}

__global__ void diffusion1(double *d_V, int *d_adj_triangle, double *d_D_CF1, double *d_D_CF2, double *d_D_CF3, double *d_gC, double *d_area, double *d_unit_CF, double *d_Tf, double *d_aF, double *d_Vnew, double *d_gradientU_C){	
	int k = threadIdx.x + blockIdx.x * blockDim.x;
    int I = threadIdx.x;

	__shared__ double gradient_f1[tpb][2];
	__shared__ double gradient_f1_tmp[tpb][2];
	__shared__ double gradient_f2[tpb][2];
	__shared__ double gradient_f2_tmp[tpb][2];
	__shared__ double gradient_f3[tpb][2];
	__shared__ double gradient_f3_tmp[tpb][2];
	__shared__ double bC[tpb][3];
	__shared__ double sum_t[tpb];
	
	if(k<nt){	
		//compute gradients at midpoins of three faces, since mid - point integration approximation is used.Section 9.4
		//gradient_fi: the interpolated gradient at the intersecting point of face i
		gradient_f1_tmp[I][0] = (1 - d_gC[k])*d_gradientU_C[d_adj_triangle[k]-1] + d_gC[k] * d_gradientU_C[k];
		gradient_f1_tmp[I][1] = (1 - d_gC[k])*d_gradientU_C[d_adj_triangle[k]-1+nt] + d_gC[k] * d_gradientU_C[k+nt];
		//Note: here a tmp variable "gradient_f1_tmp" is needed to avoid across assignment
		gradient_f1[I][0] = (gradient_f1_tmp[I][0] + ((d_V[d_adj_triangle[k]-1] - d_V[k]) / sqrt(d_D_CF1[k]*d_D_CF1[k]+d_D_CF1[k+nt]*d_D_CF1[k+nt]) - (gradient_f1_tmp[I][0] * d_unit_CF[k] + gradient_f1_tmp[I][1] * d_unit_CF[k+nt]))*d_unit_CF[k])*(d_adj_triangle[k] != 0);
		gradient_f1[I][1] = (gradient_f1_tmp[I][1] + ((d_V[d_adj_triangle[k]-1] - d_V[k]) / sqrt(d_D_CF1[k]*d_D_CF1[k]+d_D_CF1[k+nt]*d_D_CF1[k+nt]) - (gradient_f1_tmp[I][0] * d_unit_CF[k] + gradient_f1_tmp[I][1] * d_unit_CF[k+nt]))*d_unit_CF[k+nt])*(d_adj_triangle[k] != 0);
		bC[I][0] = D * (gradient_f1[I][0] * d_Tf[k] + gradient_f1[I][1] * d_Tf[k+nt]); //non-orthogonal flux
		
		gradient_f2_tmp[I][0] = (1 - d_gC[k+nt])*d_gradientU_C[d_adj_triangle[k+nt]-1] + d_gC[k+nt] * d_gradientU_C[k];
		gradient_f2_tmp[I][1] = (1 - d_gC[k+nt])*d_gradientU_C[d_adj_triangle[k+nt]-1+nt] + d_gC[k+nt] * d_gradientU_C[k+nt];
		gradient_f2[I][0] = (gradient_f2_tmp[I][0] + ((d_V[d_adj_triangle[k+nt]-1] - d_V[k]) / sqrt(d_D_CF2[k]*d_D_CF2[k]+d_D_CF2[k+nt]*d_D_CF2[k+nt]) - (gradient_f2_tmp[I][0] * d_unit_CF[k+2*nt] + gradient_f2_tmp[I][1] * d_unit_CF[k+3*nt]))*d_unit_CF[k+2*nt])*(d_adj_triangle[k+nt] != 0);
		gradient_f2[I][1] = (gradient_f2_tmp[I][1] + ((d_V[d_adj_triangle[k+nt]-1] - d_V[k]) / sqrt(d_D_CF2[k]*d_D_CF2[k]+d_D_CF2[k+nt]*d_D_CF2[k+nt]) - (gradient_f2_tmp[I][0] * d_unit_CF[k+2*nt] + gradient_f2_tmp[I][1] * d_unit_CF[k+3*nt]))*d_unit_CF[k+3*nt])*(d_adj_triangle[k+nt] != 0);
		bC[I][1] = D * (gradient_f2[I][0] * d_Tf[k+2*nt] + gradient_f2[I][1] * d_Tf[k+3*nt]);
		
		gradient_f3_tmp[I][0] = (1 - d_gC[k+2*nt])*d_gradientU_C[d_adj_triangle[k+2*nt]-1] + d_gC[k+2*nt] * d_gradientU_C[k];
		gradient_f3_tmp[I][1] = (1 - d_gC[k+2*nt])*d_gradientU_C[d_adj_triangle[k+2*nt]-1+nt] + d_gC[k+2*nt] * d_gradientU_C[k+nt];
		gradient_f3[I][0] = (gradient_f3_tmp[I][0] + ((d_V[d_adj_triangle[k+2*nt]-1] - d_V[k]) / sqrt(d_D_CF3[k]*d_D_CF3[k]+d_D_CF3[k+nt]*d_D_CF3[k+nt]) - (gradient_f3_tmp[I][0] * d_unit_CF[k+4*nt] + gradient_f3_tmp[I][1] * d_unit_CF[k+5*nt]))*d_unit_CF[k+4*nt])*(d_adj_triangle[k+2*nt] != 0);
		gradient_f3[I][1] = (gradient_f3_tmp[I][1] + ((d_V[d_adj_triangle[k+2*nt]-1] - d_V[k]) / sqrt(d_D_CF3[k]*d_D_CF3[k]+d_D_CF3[k+nt]*d_D_CF3[k+nt]) - (gradient_f3_tmp[I][0] * d_unit_CF[k+4*nt] + gradient_f3_tmp[I][1] * d_unit_CF[k+5*nt]))*d_unit_CF[k+5*nt])*(d_adj_triangle[k+2*nt] != 0);
		bC[I][2] = D * (gradient_f3[I][0] * d_Tf[k+4*nt] + gradient_f3[I][1] * d_Tf[k+5*nt]);

		sum_t[I] = 0.0;
		sum_t[I] = sum_t[I] + (d_aF[k] * d_V[d_adj_triangle[k]-1])*(d_adj_triangle[k] != 0) + (d_aF[k+nt] * d_V[d_adj_triangle[k+nt]-1])*(d_adj_triangle[k+nt] != 0) + (d_aF[k+2*nt] * d_V[d_adj_triangle[k+2*nt]-1])*(d_adj_triangle[k+2*nt] != 0);

		d_Vnew[k] = d_V[k] + (sum_t[I] + (-(d_aF[k] + d_aF[k+nt] + d_aF[k+2*nt]))*d_V[k] + bC[I][0] + bC[I][1] + bC[I][2]) / (d_area[k] / dt);
	}
}


void FDM(){
	int bpg;
	bpg = (nt+tpb-1)/tpb;
	diffusion<<<bpg, tpb>>>(d_V, d_adj_triangle, d_weight, d_D_CF1, d_D_CF2, d_D_CF3, d_matrix_G, d_gC, d_area, d_Sf, d_gradientU_C);
	diffusion1<<<bpg, tpb>>>(d_V, d_adj_triangle, d_D_CF1, d_D_CF2, d_D_CF3, d_gC, d_area, d_unit_CF, d_Tf, d_aF, d_Vnew, d_gradientU_C);
}

__global__ void Euler(double *d_V, double *d_Vnew){
	int k = threadIdx.x + blockIdx.x * blockDim.x;

	if(k<nt){
		d_V[k] = d_Vnew[k];
	}
}

void Forward_Euler(){
	int bpg;
    bpg = (nt+tpb-1)/tpb;
	Euler<<<bpg, tpb>>>(d_V, d_Vnew);
	//cudaDeviceSynchronize();
}
