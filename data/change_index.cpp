#include <stdio.h>
int const nt = 50;//12800;//3200;//800;//200;//50; //number of Triangles

int main(){

	int i;
	int elements[nt][3] = { 0 };
	FILE *fp_elements;
        fp_elements = fopen("elements_c4.dat", "r");
        if (NULL == fp_elements)
                return 2;
        for (i = 0; i < nt; i++)
	        fscanf(fp_elements, "%d %d %d\n", &elements[i][0], &elements[i][1], &elements[i][2]);
        fclose(fp_elements);
	int elements_id0[nt][3] = { 0 };

	for (i = 0; i < nt; i++){
		elements_id0[i][0] = elements[i][0]-1;
		elements_id0[i][1] = elements[i][1]-1;
		elements_id0[i][2] = elements[i][2]-1;
	}

	FILE *f;
	f = fopen("0base_elements_c4.dat", "w");
	for (i = 0; i < nt; i++){
		fprintf(f, "%d %d %d\n", elements_id0[i][0], elements_id0[i][1], elements_id0[i][2]);
	}
	fclose(f);
}
