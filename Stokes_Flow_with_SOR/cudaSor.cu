//Lionel Fiske
//Cuda SOR




#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
#include<cuda_runtime.h>

//void KernelU(double* dev_U, double* dev_maxResU, int* dev_N, double* dev_dx, double* dev_mu);


static void HandleError(cudaError_t err, const char* file, int line){

	if(err != cudaSuccess){
		printf("%s in %s at line %d \n ", cudaGetErrorString(err), file, line );
		exit(1);
	}
}

#define HANDLE_ERROR(err) (HandleError(err, __FILE__,__LINE__))


//__global__ void KernelU(double* dev_U,double* dev_P, double* dev_ResU,double* dev_numberOfBlocksU , int* dev_N, double* dev_P_const, double* dev_dx, double* dev_mu, double* dev_omega){

//__shared__ double localU[pointsPerBlock];

/*
//Since we are not interested in the residual I will keep resU resV and resP locally on card
double uC=100;
double uE=100;
double uW=100;
double uN=100;
double uS=100;
double pC=100;
double pW=100;
double P_const=dev_P_const[0];
double omega=*dev_omega;

int i;
int j;
int gi;
int gj;
int rowsPerBlock=dev_N[0]/ dev_numberOfBlocksU[0]; //The local block cN
int dx=dev_dx[0];
int mu=dev_mu[0];
int N=dev_N[0];
int g_tid=(blockIdx.x)*blockDim.x+threadIdx.x;
int tid=threadIdx.x;



for(int color=1; color<3; color++){

//Compute U Residual and update U

//Not in first row or last row 	

	if(blockIdx.x!=0 && blockIdx.x!=*dev_numberOfBlocksU){	
		if(N/2<tid<N/2*(rowsPerBlock) ){

			gi=tid%dev_N; 
			gj=tid/N;

			i=tid%N; //U has N cols
			j=tid/N;

			i=2*i + (j+color)%2; //Color selector 

			//Center points of stencil 
			uC=dev_U[i+N*j];
			pC=dev_P[i+(N-1)*j]; 


			//Handle the x left boundary 
			if(gi!=0 ){
				uW=dev_U[i-1+N*j];
				pW=dev_P[(i-1)+(N-1)*j];
			}

			//left points
			if(gi==0 ){
				uW=uC;
				pW=2*P_const - pC;
			}

			//handle x right interior points
			if(gi!=N-1){
			uE=dev_U[i+1+N*j];
			}

			//handle x right boundary points
			if(gi==N-1){
			uE=uC;
			pC=-dev_P[i-1+(N-1)*j]; 
			}


			//Handle the up interior point
			if(gj!=0 ){
				uN=dev_U[i+N*(j-1)];
			}

			//Handle the y top boundary case
			if(gj==0 ){
				uN=-uC;
			}

			//Handle the y down interior case
			if(gj!=(N-1)-1 ){
				uS=dev_U[i+N*(j+1)];
			}

			//Handle the y down boundary case
			if(gj==(N-1)-1 ){
				uS=-uC;
			}

			// Compute Residual 
				resU= mu*(uE + uW - 2*uC)+ mu*(uN + uS - 2*uC) - dx*(pC-pW);
				dev_U[i+N*j]=dev_U[i+N*j]+ omega*resU;


			}

		}



	if(blockIdx.x==0){	
		if(N/2<tid<N/2*(rowsPerBlock) ){

			gi=tid%dev_N; 
			gj=tid/N;

			i=tid%N; //U has N cols
			j=tid/N;

			i=2*i + (j+color)%2; //Color selector 

			//Center points of stencil 
			uC=dev_U[i+N*j];
			pC=dev_P[i+(N-1)*j]; 


			//Handle the x left boundary 
			if(gi!=0 ){
				uW=dev_U[i-1+N*j];
				pW=dev_P[(i-1)+(N-1)*j];
			}

			//left points
			if(gi==0 ){
				uW=uC;
				pW=2*P_const - pC;
			}

			//handle x right interior points
			if(gi!=N-1){
			uE=dev_U[i+1+N*j];
			}

			//handle x right boundary points
			if(gi==N-1){
			uE=uC;
			pC=-dev_P[i-1+(N-1)*j]; 
			}


			//Handle the up interior point
			if(gj!=0 ){
				uN=dev_U[i+N*(j-1)];
			}

			//Handle the y top boundary case
			if(gj==0 ){
				uN=-uC;
			}

			//Handle the y down interior case
			if(gj!=(N-1)-1 ){
				uS=dev_U[i+N_col*(j+1)];
			}

			//Handle the y down boundary case
			if(gj==(N-1)-1 ){
				uS=-uC;
			}

			// Compute Residual 
				resU= mu*(uE + uW - 2*uC)+ mu*(uN + uS - 2*uC) - dx*(pC-pW);
				dev_U[i+N*j]=dev_U[i+N*j]+ omega*resU;

			//Check if this is the biggest residual 


			}

		}


	if(blockIdx.x==*dev_numberOfBlocksU){	
		if(N/2<tid<N/2*(rowsPerBlock) ){

			gi=tid%dev_N; 
			gj=tid/N;

			i=tid%N; //U has N cols
			j=tid/N;

			i=2*i + (j+color)%2; //Color selector 

			//Center points of stencil 
			uC=dev_U[i+N*j];
			pC=dev_P[i+(N-1)*j]; 


			//Handle the x left boundary 
			if(gi!=0 ){
				uW=dev_U[i-1+N*j];
				pW=dev_P[(i-1)+(N-1)*j];
			}

			//left points
			if(gi==0 ){
				uW=uC;
				pW=2*P_const - pC;
			}

			//handle x right interior points
			if(gi!=N-1){
			uE=dev_U[i+1+N*j];
			}

			//handle x right boundary points
			if(gi==N-1){
			uE=uC;
			pC=-dev_P[i-1+(N-1)*j]; 
			}


			//Handle the up interior point
			if(gj!=0 ){
				uN=dev_U[i+N*(j-1)];
			}

			//Handle the y top boundary case
			if(gj==0 ){
				uN=-uC;
			}

			//Handle the y down interior case
			if(gj!=(N-1)-1 ){
				uS=dev_U[i+N_col*(j+1)];
			}

			//Handle the y down boundary case
			if(gj==(N-1)-1 ){
				uS=-uC;
			}

			// Compute Residual 
				resU= mu*(uE + uW - 2*uC)+ mu*(uN + uS - 2*uC) - dx*(pC-pW);
				dev_U[i+N*j]=dev_U[i+N*j]+ omega*resU;

			//Check if this is the biggest residual 


			}

		}

	}


*/


//}



int main(int argc, char* argv[]){

//Since the error tolerance is 1e-5 I will just use floats instead of doubles to reduce memory transfer times
// and number of blocks

// Get Passed in Params
int	   N=atoi(argv[1]); 
double mu=atof(argv[2]); 
double omega=atof(argv[3]); 
double P_const=atof(argv[4]); 
double tol=atof(argv[5]); 
double dx=1/N ;


//as in the book so in my code
	cudaDeviceProp prop;
	int dev;
	memset(&prop, 0, sizeof(cudaDeviceProp));
	prop.multiProcessorCount =13; 
	HANDLE_ERROR(cudaChooseDevice(&dev, &prop));
	HANDLE_ERROR(cudaSetDevice(dev));
	HANDLE_ERROR(cudaGetDeviceProperties(&prop,dev));
// Now props has the info I need to// my code

//Host Data


	double* u=(double*) malloc((N-1)*N*sizeof(double)); 				
	double* v=(double*) malloc((N-1)*N*sizeof(double)); 	
	double* p=(double*) malloc((N-1)*(N-1)*sizeof(double)); 		
/*
doesnt work?
	memset(u , 0, (N-1)*N*sizeof(double));	
	memset(v , 0, (N-1)*N*sizeof(double));	
	memset(p , 0, (N-1)*(N-1)*sizeof(double));	

*/
	for(int i=0; i<(N*(N-1)); i++ ){
		u[i]=0;
		v[i]=0;
		if(i<(N-1)*(N-1))
		p[i]=0;	
	}





	double currentResidual=100;

//Compute the number of blocks I need.
	int storageOfRow=N*sizeof(double);
	int totalFastMemory=prop.sharedMemPerBlock;
	int rowsPerBlock=(total_fast_memory/storageOfRow);

	//So we want each block to carry some of U V and P and they need additional storage
	// for the rows bove and below so u has N rows but we assume that rowsperblock is 2 less than it actually is 
	//then we have N for U and (N-1) for v and P 


//Handle Device Data 
	
	double* dev_u;	
	double* dev_v; 	
	double* dev_p; 

	double*	dev_N;
	double* dev_mu; 
	double* dev_dx; 
	double* dev_omega; 
	double* dev_P_const; 
	double* dev_tol; 

	double* dev_maxResU; 
	double* dev_maxResV; 
	double* dev_maxResP; 

	int numberOfBlocksU=N/(rowsPerBlock-2) ; 
	int numberOfBlocksV=(N-1)/(rowsPerBlock-2) ; 
	int numberOfBlocksP=(N-1)/(rowsPerBlock-2) ; 

	int* dev_numberOfBlocksU; 
	int* dev_numberOfBlocksV; 
	int* dev_numberOfBlocksP; 



	HANDLE_ERROR(cudaMalloc((double*) &dev_u, (N-1)*N,sizeof(double)));
	HANDLE_ERROR(cudaMalloc((double*) &dev_v, (N-1)*N,sizeof(double)));
	HANDLE_ERROR(cudaMalloc((double*) &dev_p, (N-1)*(N-1),sizeof(double)));


/*
	HANDLE_ERROR(cudaMalloc((void*) &dev_N, sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_mu, sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_dx, sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_omega, sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_P_const, sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_tol, sizeof(double)));

	HANDLE_ERROR(cudaMalloc((void*) &dev_ResU, (N-1)*N*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_ResV, (N-1)*N*sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_ResP,(N-1)*(N-1)*sizeof(double)));

	HANDLE_ERROR(cudaMalloc((void*) &dev_numberOfBlocksU, sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_numberOfBlocksV, sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void*) &dev_numberOfBlocksP, sizeof(int)));


	HANDLE_ERROR(cudaMemcpy( dev_N, &N, sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_mu, &mu, sizeof(double), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_dx, &dx, sizeof(double), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_omega, &omega, sizeof(double), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_P_const, &P_const, sizeof(double), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_tol, &tol, sizeof(double), cudaMemcpyHostToDevice));

	//cudaMemcpy( dev_ResU, &currentResidual, (N-1)*N*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy( dev_ResV, &currentResidual, (N-1)*N*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy( dev_ResP, &currentResidual, (N-1)*(N-1)*sizeof(double), cudaMemcpyHostToDevice);

	HANDLE_ERROR(cudaMemcpy( dev_numberOfBlocksU, &numberOfBlocksU, sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_numberOfBlocksV, &numberOfBlocksV, sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( dev_numberOfBlocksP, &numberOfBlocksP, sizeof(int), cudaMemcpyHostToDevice));




*/


	free(u);
	free(v);
	free(p);

/*
	cudaFree(dev_u);
	cudaFree(dev_v);
	cudaFree(dev_p);




	cudaFree(dev_numberOfBlocksU);
	cudaFree(dev_numberOfBlocksV);
	cudaFree(dev_numberOfBlocksP);

	cudaFree(dev_N);
	cudaFree(dev_mu);
	cudaFree(dev_dx);
	cudaFree(dev_omega);
	cudaFree(dev_P_const);
	cudaFree(dev_omega);
	cudaFree(dev_tol);
	cudaFree(dev_maxResU);
	cudaFree(dev_maxResV);
	cudaFree(dev_maxResP);
*/

}