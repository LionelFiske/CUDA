
LIBS = -llapack -framework Accelerate  
#INCLUDE = -I/usr/local/Cellar/lapack/3.8.0/include 


cudaSor: cudaSor.o  
	nvcc -o cudaSor cudaSor.o $(LIBS) -g 

cudaSOR.o: cudaSor.c 
	nvcc -c cudaSor.c -g 

clean:
	rm cudaSor.o cudaSor cudaSor.out 
