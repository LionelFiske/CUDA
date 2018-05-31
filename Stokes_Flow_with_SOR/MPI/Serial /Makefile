
LIBS =  -L/usr/local/Cellar/lapack/3.8.0 \
 -llapack -framework Accelerate  
INCLUDE = -I/usr/local/Cellar/lapack/3.8.0/include 


StokesStandard: StokesStandard.o  
	mpicc -o StokesStandard StokesStandard.o $(LIBS) -g

StokesStandard.o: StokesStandard.c 
	mpicc -c StokesStandard.c -g

clean:
	rm StokesStandard.o StokesStandard StokesStandard.out 
