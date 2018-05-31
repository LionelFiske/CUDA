//Lionel Fiske
//ESAPPM 444 - intro to hating lapack
// This code is working using serial jacobi with the same scheme as assigned 
// use omega less than .2 and 50k its for convergence. 



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Accelerate/Accelerate.h>
#include <mpi.h>


void extractBlack(double* M, double* M_B , int N_row , int N_col );
void insertBlack(double* M, double* M_B , int N_row , int N_col);
void extractRed(double* M, double* M_R , int N_row , int N_col);
void insertRed(double* M_R, double* M , int N_row , int N_col);
void printMat(double* M, int rows, int cols);

void residualU(double* resU,double* U,double* P , int N_row, int N_col, double mu, double dx, double P_const);
void residualV(double* resV,double* V,double* P , int N_row, int N_col, double mu, double dx, double P_const);
void residualP(double* resP, double* U, double* V, double* P , int N_row, int N_col, double mu, double dx, double P_const);
//void residualP(double* U,double* V,double* P , int N, int M);

void updateRule(double* M, double* resM , int start_row, int stop_row ,int N_col , int N_row, double omega);


int main(int argc, const char * argv[]) {


//MPI setup
	MPI_Init(&argc , (char***) &argv) ;
	int rank;
	int numCores;
	MPI_Comm_size(MPI_COMM_WORLD,&numCores); //How many chunks
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  //Get your rank

if(rank==0){
	printf("numcores= %d \n", numCores );
}
 
//Scalar Declarations

	double dx; //I am assuming dx and dy are same 
	double L;
	double mu;
	//double lambda;
	double omega;
	int N;
	int P_const;
	int number_steps;
	int wFreq;
	int T_end=2;

//I am assuming only rank 0 has IO priv. so I am B_Casting my argins. 
	if(rank==0){
		N=atoi(argv[1]); 
		//number_steps=atoi(argv[2]); 
		mu=atof(argv[2]); 
		omega=atof(argv[3]); 
		P_const=atof(argv[4]); 
	}
		MPI_Bcast( &N, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Bcast( &number_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast( &mu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast( &omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		int localN=N/numCores ;

// Matrix Declarations
	/*
	double* u; 		
	double* uBlack; 		
	double* uRed;

	double* v; 		
	double* vBlack;
	double* vRed;

	double* p; 		
	double* pBlack;
	double* pRed;
	*/
/////////////////////////////////////////////////

//Allocations
	L=1;
	dx=(double)L/N;
	//dt=(double) T_end/number_steps;
	//lambda=(double) dt/(dx);
	wFreq=number_steps/10;

printf("dx=%lf, and N=%d \n", dx ,N);
//Notes
//u \in \R^{NxN-1}
//u_black \in \R^{Nx floor( (N-1)/2 ) }
//u_red \in \R^{Nx ceil( (N-1)/2 ) }


//v \in \R^{N-1 x N}
//v_black \in \R^{N-1 x ( (N)/2 ) }
//v_red \in \R^{N-1 x ( (N)/2 ) }


//p \in \R^{N+1 xN-1}
//p_black \in \R^{N+1 x floor( (N-1)/2 ) }
//p_red \in \R^{N+1 x ceil( (N-1)/2 ) }


	//N=10;
	double chN1 = ceil( (double) (N)/2  ); //This is always N I am dumb
	double fhN1 = floor( (double) (N-1)/2  );

	//printf("ciel is %lf\n", chN1 );
	//printf("floor is %lf\n", fhN1 );

	double* u=calloc((N-1)*N,sizeof(double)); 		
	double* uBlack=calloc( (int) fhN1 *N,sizeof(double)); 			
	double* uRed=calloc((N/2)*N,sizeof(double)); 
	double* resU=calloc((N-1)*N,sizeof(double)); 		

	double* v=calloc((N-1)*N,sizeof(double)); 		
	double* vBlack=calloc((N-1)*N/2,sizeof(double)); 	//When conpressing data to send I get half the cols (assuming N is even)
	double* vRed=calloc((N-1)*N/2,sizeof(double)); 	
	double* resV=calloc((N-1)*N,sizeof(double)); 	

	double* p=calloc((N-1)*(N-1),sizeof(double)); 		
	double* pBlack=calloc(((N-1)*(N-1) +1)/2,sizeof(double));
	double* pRed=calloc(((N-1)*(N-1) +1)/2,sizeof(double));
	double* resP=calloc((N-1)*(N-1),sizeof(double)); 	



	double* u_exact=calloc((N-1)*N,sizeof(double)); 

	double* p_exact=calloc((N-1)*(N-1),sizeof(double)); 		


	//double* recieveBuffer=calloc((N-1)*N/2,sizeof(double));


// Save location
	 FILE *fileidU ;
	 FILE *fileidV ;
	 FILE *fileidP ;
	if(rank==0){
     fileidU = fopen( "StokesU.out" , "w" ) ; 	
     fileidV = fopen( "StokesV.out" , "w" ) ; 	    
     fileidP = fopen( "StokesP.out" , "w" ) ; 	
}




for(int i=0; i<50000; i++){
//printMat(u,N-1,N);

residualU( resU, u, p , N-1, N, mu, dx, P_const);
updateRule(u,resU ,0, N-1 ,N-1, N,  omega);

//printMat(resU,N-1,N);

residualV( resV, v, p , N, N-1, mu, dx, P_const);
updateRule(v,resV ,0, N-1 ,N-1 , N,  omega);

//printf("v \n" );
//printMat(v,N,N-1);

residualP( resP, u, v, p , N-1, N-1, mu, dx, P_const);
updateRule(p,resP ,0, N-1 ,N-1 , N-1,  omega);

//printf("p \n" );
//printMat(p,N-1,N-1);


}




//////////////////////////////////////////////////////
// Save final result.
//MPI_Barrier( MPI_COMM_WORLD);

	if(rank==0){
	printf("last step \n");

	fwrite(u, sizeof( double ) , N*(N-1) , fileidU ) ; 
	fwrite(v, sizeof( double ) , N*(N-1) , fileidV ) ; 
	fwrite(p, sizeof( double ) , (N-1)*(N-1) , fileidP ) ; 

    fclose(fileidU);
    fclose(fileidV);
    fclose(fileidP);
	}	


////////////////////////////////////////////////////
//No memory leaks here. 

	free(u_exact);


	free(u);
	free(uRed);
	free(uBlack);

	free(v);
	free(vRed);
	free(vBlack);

	free(p);
	free(pRed);
	free(pBlack);
	MPI_Finalize();


}

void printMat(double* M, int rows, int cols){

//Prints matrices 

for(int j=0; j<rows; j++ ){

	for(int i=0; i<cols; i++ ){

		printf("%.4lf ", M[i+cols*j] );

	}
	printf("\n ");
	}

		printf("\n");
}



void residualU(double* resU,double* U, double* P , int N_row, int N_col, double mu, double dx, double P_const  ){


// Set these to random numbers 
	double uC=-420;
	double uE=666;
	double uW=420;
	double uN=-666;
	double uS=420;
	double pC=666;
	double pW=420;

	int N_col_p=N_col -1; 
	


		for(int j=0; j<N_row; j++){
			for(int i=0; i<(N_col) ; (i++) ){




				uC=U[i+N_col*j];
				pC=P[i+N_col_p*j];


				//Handle the x Boundary cases
				if(i!=0 ){
					uW=U[i-1+N_col*j];
					pW=P[(i-1)+N_col_p*j];
				}

				if(i==0 ){
					uW=uC;
					pW=2*P_const - pC;
				}

				if(i!=N_col-1){
					uE=U[i+1+N_col*j];
				}

				if(i==N_col-1){
					uE=uC;
					pC=-P[i-1+N_col_p*j]; 
				}


				//Handle the y Boundary cases
				if(j!=0 ){
					uN=U[i+N_col*(j-1)];
				}

				if(j==0 ){
					uN=-uC;
				}

				if(j!=N_row-1 ){
					uS=U[i+N_col*(j+1)];
				}

				if(j==N_row-1 ){
					uS=-uC;
				}

				///////////////////////// Compute Residual 

				resU[i+N_col*j]= mu*(uE + uW - 2*uC)+ mu*(uN + uS - 2*uC) - dx*(pC-pW);

				
				/////////////////////////
		}
	}

}





void residualV(double* resV,double* V,double* P , int N_row, int N_col, double mu, double dx, double P_const){
	
	double vC=-420;
	double vE=666;
	double vW=420;
	double vN=-666;
	double vS=420;
	double pC=666;
	double pN=420;
	//double pN=420;

	int N_col_p=N_col; 
	


		for(int j=0; j<N_row; j++){
			for(int i=0; i<(N_col) ; (i++) ){
//printf("%lf", V[i+N_col*j]);

				vC=V[i+N_col*j];
				pC=P[i+N_col_p*j];


				//Handle the x Boundary cases
				if(i!=0 ){
					vW=V[i-1+N_col*j];
				}

				if(i==0 ){
					vW=vC;
		
				}

				if(i!=N_col-1){
					vE=V[i+1+N_col*j];

				}

				if(i==N_col-1){
					vE=vC;
				}


				//Handle the y Boundary cases
				if(j!=0 ){
					vN=V[i+N_col*(j-1)];
					pN=P[i+N_col_p*(j-1)];

				}

				if(j!=N_row-1 ){
					vS=V[i+N_col*(j+1)];
				}

				///If V on y boundary we force it to remain 0
				
				if(j==N_row-1 ){
						 vC=0;
						 vE=0;
						 vW=0;
						 vN=0;
						 vS=0;
						 pC=0;
						 pN=0;
				}


				if(j==0 ){
						 vC=0;
						 vE=0;
						 vW=0;
						 vN=0;
						 vS=0;
						 pC=0;
						 pN=0;
				}

				///////////////////////// Compute Residual 
				resV[i+N_col*j]= mu*(vE + vW - 2*vC)  + mu*(vN + vS - 2*vC) - dx*(pC-pN);
				/////////////////////////
			}
		}
	}



//residualP( resP, u, v, p , N+1, N-1, mu, dx, P_const, color);

void residualP(double* resP, double* U, double* V, double* P , int N_row, int N_col, double mu, double dx, double P_const){

	double vC=-420;
	double vS=666;
	double uC=420;
	double uE=-666;


	int N_col_u=N_col+1; 
	


		for(int j=0; j<N_row; j++){
			for(int i=0; i<(N_col) ; (i++) ){




				vC=V[i+N_col*j];
				vS=V[i+N_col*(j+1)];

				uC=U[i+N_col_u*j];
				uE=U[i+1+N_col_u*j];


				///////////////////////// Compute Residual 

				resP[i+N_col*j]=-(uE-uC) - (vS-vC);

				/////////////////////////
			}
		}

}







void updateRule(double* M, double* resM , int start_row, int stop_row ,int N_row , int N_col, double omega){


// if N_col%2 ==0 then the number of red or black in each row does not change 
// if N_col%2 ==1 then the number of red or black alternates between 0 and 1. 	


		for(int j=0; j<N_row; j++){
			for(int i=0; i<(N_col) ; i++){


				M[i+(N_col)*j]= M[i+(N_col)*j] + omega* resM[i+(N_col)*j];

			}
		}

}













