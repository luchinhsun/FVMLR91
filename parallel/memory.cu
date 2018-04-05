#include "head.h"

double *h_t;
double *d_t;
double *h_V;
double *d_V;
double *d_dV2;
double *d_Vnew;
double *d_it;

double *h_m;
double *d_m;
double *h_h;
double *d_h;
double *h_jj;
double *d_jj;
double *h_d;
double *d_d;
double *h_f;
double *d_f;
double *h_X;
double *d_X;
double *h_cai;
double *d_cai;

double *h_it;

double *dcai;

/* FVM variable*/
int *h_stim_i;
int *h_adj_triangle; // save the index of  three neighbor triangles
double *h_weight; //weighting factor
double *h_D_CF1;
double *h_D_CF2;
double *h_D_CF3;
double *h_matrix_G; //for Least - Square method in Section 9.3
double *h_gC; //geometric interpolation factor related to the position of the element face f with respect to the nodes C and F.
double *h_area; // each triangle's h_area
double *h_Sf;//surface vector
double *h_unit_CF;
double *h_Tf; // component Tf is normal to Sf
double *h_aF;

int *d_stim_i;
int *d_adj_triangle; // save the index of  three neighbor triangles
double *d_weight; //weighting factor
double *d_D_CF1;
double *d_D_CF2;
double *d_D_CF3;
double *d_matrix_G; //for Least - Square method in Section 9.3
double *d_gC; //geometric interpolation factor related to the position of the element face f with respect to the nodes C and F.
double *d_area; // each triangle's h_area
double *d_Sf;//surface vector
double *d_unit_CF;
double *d_Tf; // component Tf is normal to Sf
double *d_aF;

double *d_gradientU_C;

void Allocate(){
	cudaError_t Error;
	size_t size = nt*sizeof(double);

	h_t = (double*)malloc(size);
	Error = cudaMalloc((void**)&d_t, size);
	printf("CUDA error = %s\n",cudaGetErrorString(Error));

	h_V = (double*)malloc(size);
	cudaMalloc((void**)&d_V, size);
	cudaMalloc((void**)&d_dV2, size);
	cudaMalloc((void**)&d_Vnew, size);

	cudaMalloc((void**)&d_it, size);

	h_m = (double*)malloc(size);
	cudaMalloc((void**)&d_m, size);
	h_h = (double*)malloc(size);
        cudaMalloc((void**)&d_h, size);
	h_jj = (double*)malloc(size);
        cudaMalloc((void**)&d_jj, size);
	h_d = (double*)malloc(size);
        cudaMalloc((void**)&d_d, size);
	h_f = (double*)malloc(size);
        cudaMalloc((void**)&d_f, size);
	h_X = (double*)malloc(size);
        cudaMalloc((void**)&d_X, size);
	h_cai = (double*)malloc(size);
        cudaMalloc((void**)&d_cai, size);

	h_it = (double*)malloc(size);

	cudaMalloc((void**)&dcai, size);
	
	h_stim_i = (int *)malloc(nt*sizeof(int));
	h_adj_triangle = (int *)malloc(nt*3*sizeof(int));
	h_weight = (double *)malloc(nt*3*sizeof(double));
	h_D_CF1 = (double *)malloc(nt*2*sizeof(double));
	h_D_CF2 = (double *)malloc(nt*2*sizeof(double));
	h_D_CF3 = (double *)malloc(nt*2*sizeof(double));
	h_matrix_G = (double *)malloc(nt*2*2*sizeof(double));
	h_gC = (double *)malloc(nt*3*sizeof(double));
	h_area = (double *)malloc(nt*sizeof(double));
	h_Sf = (double *)malloc(nt*3*2*sizeof(double));
	h_unit_CF = (double *)malloc(nt*3*2*sizeof(double));
	h_Tf = (double *)malloc(nt*3*2*sizeof(double));
	h_aF = (double *)malloc(nt*3*sizeof(double));
	
	cudaMalloc((void**)&d_stim_i, nt*sizeof(int));
	cudaMalloc((void**)&d_adj_triangle, nt*3*sizeof(int));
	cudaMalloc((void**)&d_weight, nt*3*sizeof(double));
	cudaMalloc((void**)&d_D_CF1, nt*2*sizeof(double));
	cudaMalloc((void**)&d_D_CF2, nt*2*sizeof(double));
	cudaMalloc((void**)&d_D_CF3, nt*2*sizeof(double));
	cudaMalloc((void**)&d_matrix_G, nt*2*2*sizeof(double));
	cudaMalloc((void**)&d_gC, nt*3*sizeof(double));
	cudaMalloc((void**)&d_area, nt*sizeof(double));
	cudaMalloc((void**)&d_Sf, nt*3*2*sizeof(double));
	cudaMalloc((void**)&d_unit_CF, nt*3*2*sizeof(double));
	cudaMalloc((void**)&d_Tf, nt*3*2*sizeof(double));
	cudaMalloc((void**)&d_aF, nt*3*sizeof(double));
	
	cudaMalloc((void**)&d_gradientU_C, nt*2*sizeof(double));	
}

void Free(){

	free(h_t);free(h_V);free(h_m);free(h_h);
	free(h_jj);free(h_d);free(h_f);free(h_X);free(h_cai);
	free(h_it);

	cudaFree(d_t);cudaFree(d_V);cudaFree(d_dV2);cudaFree(d_Vnew);cudaFree(d_it);
	cudaFree(d_m);cudaFree(d_h);cudaFree(d_jj);cudaFree(d_d);
	cudaFree(d_f);cudaFree(d_X);cudaFree(d_cai);

	cudaFree(dcai);

	free(h_stim_i);free(h_adj_triangle);free(h_weight);
	free(h_D_CF1);free(h_D_CF2);free(h_D_CF3);
	free(h_matrix_G);free(h_gC);free(h_area);
	free(h_Sf);free(h_unit_CF);free(h_Tf);free(h_aF);
	
	cudaFree(d_stim_i);cudaFree(d_adj_triangle);cudaFree(d_weight);
	cudaFree(d_D_CF1);cudaFree(d_D_CF2);cudaFree(d_D_CF3);
	cudaFree(d_matrix_G);cudaFree(d_gC);cudaFree(d_area);
	cudaFree(d_Sf);cudaFree(d_unit_CF);cudaFree(d_Tf);cudaFree(d_aF);
	cudaFree(d_gradientU_C);
}

void Send_to_Device(){
        cudaError_t Error;
        size_t size;
        size = nt*sizeof(double);

	Error = cudaMemcpy(d_t, h_t, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_t->d_t) = %s\n",cudaGetErrorString(Error));
        Error = cudaMemcpy(d_V, h_V, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_V->d_V) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_m, h_m, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_m->d_m) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_h, h_h, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_h->d_h) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_jj, h_jj, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_jj->d_jj) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_d, h_d, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_d->d_d) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_f, h_f, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_f->d_f) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_X, h_X, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_X->d_X) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_cai, h_cai, size, cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_cai->d_cai) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_stim_i, h_stim_i, nt*sizeof(int), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_adj_triangle, h_adj_triangle, nt*3*sizeof(int), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_weight, h_weight, nt*3*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_D_CF1, h_D_CF1, nt*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_D_CF2, h_D_CF2, nt*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_D_CF3, h_D_CF3, nt*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_matrix_G, h_matrix_G, nt*2*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_gC, h_gC, nt*3*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_area, h_area, nt*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_Sf, h_Sf, nt*3*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_unit_CF, h_unit_CF, nt*3*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_Tf, h_Tf, nt*3*2*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));
	Error = cudaMemcpy(d_aF, h_aF, nt*3*sizeof(double), cudaMemcpyHostToDevice);
        if (Error != cudaSuccess)
        printf("CUDA error(copy h_stim_i->d_stim_i) = %s\n",cudaGetErrorString(Error));

}

void Send_V(){
        cudaError_t Error;
        size_t size;
        size = nt*sizeof(double);

        Error = cudaMemcpy(h_V, d_V, size, cudaMemcpyDeviceToHost);
        if (Error != cudaSuccess)
        printf("CUDA error(copy d_V->h_V) = %s\n",cudaGetErrorString(Error));
}
