CUDA_LIB :=/usr/local/cuda/lib64 -lcuda -lcudart

All: A B C
	g++ LR91_FVM.o memory.o forward.o -o test.run -L $(CUDA_LIB)
A:
	nvcc memory.cu --ptxas-options=-v -c
B:
	nvcc forward.cu --ptxas-options=-v -c
C:
	g++ LR91_FVM.cpp -c
clean:
	rm *.o snapshot_AP* *.run single_AP1 single_AP2 fevaluation
