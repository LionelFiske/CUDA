//Lionel Fiske
//ESAPPM 444 - intro to hating lapack
//
// In this code I solve the Stokes flow using SOR/
//call as  mpirun -np <blah> Stokes N mu omega p_const
//Super well commented this time, I promise. 
//Code assumes N is even.



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

void residualU(double* resU,double* U,double* P ,int start_row, int stop_row, int N_row, int N_col, double mu, double dx, double P_const,int color);
void residualV(double* resV,double* V,double* P ,int start_row, int stop_row, int N_row, int N_col, double mu, double dx, double P_const, int color);
void residualP(double* resP, double* U, double* V, double* P ,int start_row, int stop_row, int N_row, int N_col, double mu, double dx, double P_const, int color);
//void residualP(double* U,double* V,double* P , int N, int M);

void updateRule(double* M, double* resM , int start_row, int stop_row ,int N_col , int N_row, double omega, int color);


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
	double P_const;
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
		MPI_Bcast( &P_const, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

	//N=10;
	//double chN1 = ceil( (double) (N)/2  ); //This is always N I am dumb
	//double fhN1 = floor( (double) (N-1)/2  );

	//printf("ciel is %lf\n", chN1 );
	//printf("floor is %lf\n", fhN1 );

	double* u=calloc((N-1)*N,sizeof(double)); 		
	double* uSendBuffer=calloc((N/2) *(N-1),sizeof(double)); 	
	double* uBuffer=calloc( (N-1) *N/2,sizeof(double)); 			
	double* resU=calloc((N-1)*N,sizeof(double)); 		

	double* v=calloc((N-1)*N,sizeof(double)); 	
	double* vSendBuffer=calloc((N)*N/2,sizeof(double)); 
	double* vBuffer=calloc((N)*N/2,sizeof(double)); 		//When conpressing data to send I get half the cols (assuming N is even)
	double* resV=calloc((N-1)*N,sizeof(double)); 	

	double* p=calloc((N-1)*(N-1),sizeof(double)); 		
	double* pSendBuffer=calloc(((N-1)*(N) )/2,sizeof(double));
	double* pBuffer=calloc(((N-1)*(N) )/2,sizeof(double));
	double* resP=calloc((N-1)*(N-1),sizeof(double)); 	


	int* sendCountsU=malloc(numCores*sizeof(int));
	int* sendCountsV=malloc(numCores*sizeof(int));
	int* sendCountsP=malloc(numCores*sizeof(int));
	int* offsetsU=malloc(numCores*sizeof(int));
	int* offsetsV=malloc(numCores*sizeof(int));
	int* offsetsP=malloc(numCores*sizeof(int));



// Set up MPI all gather variables

	//U stuff
	offsetsU[0]=0;
	sendCountsU[numCores-1]=(localN-1)*N/2;
	for(int i=0; i<numCores-1; i++){
		sendCountsU[i]=localN*N/2;

		if(i!=0)
		offsetsU[i]=offsetsU[i-1]+localN*N/2; 

	}
		offsetsU[numCores-1]=offsetsU[numCores-2]+localN*N/2; 
		sendCountsU[numCores-1]=(localN-1)*N/2;

	//V stuff
	offsetsV[0]=0;
	for(int i=0; i<numCores; i++){
		sendCountsV[i]=localN*(N)/2;

		if(i!=0)
		offsetsV[i]=offsetsV[i-1]+localN*(N)/2; 
	}


		//P stuff
	offsetsP[0]=0;
	//sendCountsP[numCores-1]=(localN-1)*(N)/2;
	for(int i=0; i<numCores-1; i++){
		sendCountsP[i]=(localN)*(N)/2; //need a localn+1 and am not syre why 

	if(i!=0)
		offsetsP[i]=offsetsP[i-1]+(localN)*(N)/2; 

	}

	offsetsP[numCores-1]=offsetsP[numCores-2]+(localN)*(N)/2;
	sendCountsP[numCores-1]=(localN-1)*(N)/2;


	if(rank==1){
printf("localN = %d\n",localN );

printf("offsetsU = %d\n",offsetsU[0] );
printf("offsetsV = %d\n",offsetsV[0] );
printf("offsetsP = %d\n",offsetsP[0] );

printf("offsetsU = %d\n",offsetsU[1] );
printf("offsetsV = %d\n",offsetsV[1] );
printf("offsetsP = %d\n",offsetsP[1] );

}

//printf("sendCountsU = %d\n",sendCountsU[0] );
//printf("sendCountsV = %d\n",sendCountsV[0] );
//printf("sendCountsP = %d\n",sendCountsP[0] );

//printf("array size 1= %d", (int) fhN1 *N);
//printf("array size 1= %d",  (N-1)*N/2);
//printf("array size 1= %d", (N-1)*(N) /2);

// Save location
	 FILE *fileidU ;
	 FILE *fileidV ;
	 FILE *fileidP ;
	if(rank==0){
     fileidU = fopen( "StokesU.out" , "w" ) ; 	
     fileidV = fopen( "StokesV.out" , "w" ) ; 	    
     fileidP = fopen( "StokesP.out" , "w" ) ; 	
}

int color;



for(int i=0; i<10000; i++){	


color=1;
//Update Black

		if(rank==numCores-1){
residualU( resU, u, p ,rank*(localN) ,(rank+1)*(localN)-1, N-1, N, mu, dx, P_const, color);
updateRule(u,resU ,rank*(localN), (rank+1)*(localN)-1  ,N-1, N,  omega,  color);
		}

		if(rank!=numCores-1){
residualU( resU, u, p ,rank*(localN) ,(rank+1)*(localN), N-1, N, mu, dx, P_const, color);
updateRule(u,resU ,rank*(localN), (rank+1)*(localN) ,N-1, N,  omega,  color);
		}

residualV( resV, v, p ,rank*(localN), (rank+1)*(localN), N, N-1, mu, dx, P_const, color);
updateRule(v,resV ,rank*(localN), (rank+1)*(localN) ,N , N-1,  omega,  color);


	if(rank==numCores-1){
residualP( resP, u, v, p ,rank*(localN) ,(rank+1)*(localN)-1, N-1, N-1, mu, dx, P_const, color);
updateRule(p,resP ,rank*(localN) ,(rank+1)*(localN)-1 ,N-1 , N-1,  omega,  color);

//printf("sendcounts= %d \n" ,sendCountsP[1]);
//printf("stop= %d \n" ,(rank+1)*(localN)-1 );
//printMat(p,N-1,N-1);

	}

	if(rank!=numCores-1){
residualP( resP, u, v, p ,rank*(localN) ,(rank+1)*(localN), N-1, N-1, mu, dx, P_const, color);
updateRule(p,resP ,rank*(localN) ,(rank+1)*(localN) ,N-1 , N-1,  omega,  color);

//printf("core 0 start = %d, stop = %d \n", rank*(localN), (rank+1)*(localN));

	}

//Extract the black points 

extractBlack(u, uSendBuffer , N-1 ,N );
extractBlack(v, vSendBuffer , N ,N-1 );
extractBlack(p, pSendBuffer , N-1 ,N-1 );

//Consolidate the black points 
MPI_Barrier( MPI_COMM_WORLD);
if(rank==0){
 
 //printf("Counts= %d,  %d \n", sendCountsP[0], sendCountsP[1]);

}

MPI_Barrier( MPI_COMM_WORLD);

if(rank==1){
 //printf("counts= %d,  %d \n", sendCountsP[0], sendCountsP[1]);
//	printMat(pSendBuffer,N-1,N/2);

}


MPI_Allgatherv(&uSendBuffer[offsetsU[rank]],sendCountsU[rank], MPI_DOUBLE,uBuffer,sendCountsU, offsetsU, MPI_DOUBLE, MPI_COMM_WORLD);
MPI_Allgatherv(&vSendBuffer[offsetsV[rank]],sendCountsV[rank], MPI_DOUBLE,vBuffer,sendCountsV, offsetsV, MPI_DOUBLE, MPI_COMM_WORLD);
MPI_Allgatherv(&pSendBuffer[offsetsP[rank]],sendCountsP[rank], MPI_DOUBLE,pBuffer,sendCountsP, offsetsP, MPI_DOUBLE, MPI_COMM_WORLD);


//Insert the black points into arrays

if(rank==0){
//printf(" after black buff /n ");
//printMat(pBuffer,N-1,N/2);
}



insertBlack(uBuffer, u , N-1 ,N );
insertBlack(vBuffer, v, N ,N-1 );
insertBlack(pBuffer, p, N-1 ,N-1 );
if(rank==0){
//	printMat(pBuffer,N-1,N/2);
//	printMat(p,N-1,N-1);
}

color=2;

//Update Red


		if(rank==numCores-1){
residualU( resU, u, p ,rank*(localN) ,(rank+1)*(localN)-1, N-1, N, mu, dx, P_const, color);
updateRule(u,resU ,rank*(localN), (rank+1)*(localN)-1  ,N-1, N,  omega,  color);
//printf("core 1 start = %d, stop = %d \n", rank*(localN), (rank+1)*(localN)-1);

//printMat(resU,N-1,N);
		}

		if(rank!=numCores-1){
residualU( resU, u, p ,rank*(localN) ,(rank+1)*(localN), N-1, N, mu, dx, P_const, color);
updateRule(u,resU ,rank*(localN), (rank+1)*(localN) ,N-1, N,  omega,  color);
		}

residualV( resV, v, p ,rank*(localN), (rank+1)*(localN), N, N-1, mu, dx, P_const, color);
updateRule(v,resV ,rank*(localN), (rank+1)*(localN) ,N , N-1,  omega,  color);


	if(rank==numCores-1){

residualP( resP, u, v, p ,rank*(localN) ,(rank+1)*(localN)-1, N-1, N-1, mu, dx, P_const, color);
updateRule(p,resP ,rank*(localN) ,(rank+1)*(localN) -1 ,N-1 , N-1,  omega,  color);

	}
MPI_Barrier(MPI_COMM_WORLD);
	if(rank!=numCores-1){
residualP( resP, u, v, p ,rank*(localN) ,(rank+1)*(localN), N-1, N-1, mu, dx, P_const, color);
updateRule(p,resP ,rank*(localN) ,(rank+1)*(localN) ,N-1 , N-1,  omega,  color);

//printf("repPres on 1 ");
//printMat(resP,N-1,N-1);

	}


//Extract the red points 
extractRed(u, uSendBuffer , N-1 ,N );
extractRed(v, vSendBuffer , N ,N-1 );
extractRed(p, pSendBuffer , N-1 ,N-1 );

//if(rank==1){
//printf("pSendBuffer \n ");
//printMat(pSendBuffer,N-1,N/2);
//}


//Consolidate the red points 

MPI_Allgatherv(&uSendBuffer[offsetsU[rank]],sendCountsU[rank], MPI_DOUBLE,uBuffer,sendCountsU, offsetsU, MPI_DOUBLE, MPI_COMM_WORLD);
MPI_Allgatherv(&vSendBuffer[offsetsV[rank]],sendCountsV[rank], MPI_DOUBLE,vBuffer,sendCountsV, offsetsV, MPI_DOUBLE, MPI_COMM_WORLD);
MPI_Allgatherv(&pSendBuffer[offsetsP[rank]],sendCountsP[rank], MPI_DOUBLE,pBuffer,sendCountsP, offsetsP, MPI_DOUBLE, MPI_COMM_WORLD);


//Insert the red points into arrays

//if(rank==0){
//printf("pBuffer ");
//printMat(pBuffer,N-1,N/2);
//}

insertRed(uBuffer, u, N-1 ,N );
insertRed(vBuffer, v , N ,N-1 );
insertRed(pBuffer, p , N-1 ,N-1 );



}

//////////////////////////////////////////////////////
// Save final result.
MPI_Barrier( MPI_COMM_WORLD);
//MPI_Barrier(MPI)
	if(rank==0){
	//printf("last step \n");

	fwrite(u, sizeof( double ) , N*(N-1) , fileidU ) ; 
	fwrite(v, sizeof( double ) , N*(N-1) , fileidV ) ; 
	fwrite(p, sizeof( double ) , (N-1)*(N-1) , fileidP ) ; 

    fclose(fileidU);
    fclose(fileidV);
    fclose(fileidP);
	}	


////////////////////////////////////////////////////
//No memory leaks here. 

	//free(u_exact);


	free(u);
	free(uSendBuffer);
	free(resU);

	free(v);
	free(vSendBuffer);
	free(resV);


	free(p);
	free(pSendBuffer);
	free(resP);

	free(uBuffer);
	free(vBuffer);
	free(pBuffer);

	free(sendCountsU);
	free(sendCountsV);
	free(sendCountsP);

	free(offsetsU);
	free(offsetsV);
	free(offsetsP);


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


void extractBlack(double* M, double* M_B , int N_row , int N_col ){

	if(N_col%2 ==0){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<N_col/2; i++){
				M_B[i+(N_col/2) *j] = M[2*i + (1-j%2) +N_col*j];
			}
		}
	}
	
	if(N_col%2 ==1){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<( (N_col+1 )/2  - (j+1)%2); i++){
				M_B[i+((N_col+1)/2)*j] = M[2*i + (1-j%2) +N_col*j];
			}
		}
	}


}

void insertBlack(double* M_B, double* M , int N_row, int N_col){
	
	if(N_col%2 ==0){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<N_col/2; i++){
				M[(2*i) +(1-j%2) +(N_col*j)]=M_B[i+(N_col/2)*j];
			}
		}	
	}

	if(N_col%2 ==1){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<( (N_col+1 )/2  - (j+1)%2); i++){

				M[2*i + (1-j%2) +N_col*j]=M_B[i+((N_col+1)/2) *j];

			}
		}
	}

}




void extractRed(double* M, double* M_R ,  int N_row , int N_col){


	if(N_col%2 ==0){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<N_col/2; i++){
				M_R[i+(N_col/2) *j] = M[2*i + (j%2) +N_col*j];
			}
		}
	}
	
	if(N_col%2 ==1){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<( (N_col+1 )/2  - (j)%2); i++){
				//printf("index= %d \n" , i + ((N_col+1)/2) *j);
				//printf("jndex= %d \n" , i+((N_col+1)/2) *j);
				M_R[i + ((N_col+1)/2) *j] = M[2*i + (j%2) +N_col*j];
				//printf("val= %lf \n" ,M[2*i + (j%2) +N_col*j]);
			}
		}
	}



}	



//insertRed(uBlack, u , 6 , 5);

void insertRed(double* M_R, double* M ,  int N_row , int N_col){

	if(N_col%2 ==0){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<N_col/2; i++){
			M[2*i + (j%2) +N_col*j]=M_R[i+(N_col/2) *j];
			//printf(" FUck lslslsl jndex= %d \n" , i+((N_col+1)/2) *j);
			}
		}
	}
	
	if(N_col%2 ==1){
		for(int j=0; j<N_row; j++){
			for(int i=0; i<( (N_col+1 )/2  - (j)%2); i++){
			//printf("index= %d \n" , i + ((N_col+1)/2) *j);
				//printf(" FUck lslslsl jndex= %d \n" , i+((N_col+1)/2) *j);
				M[2*i + (j%2) +N_col*j]=M_R[i + ((N_col+1)/2) *j];
				//printf("val= %lf \n" ,M[2*i + (j%2) +N_col*j]);
			}
		}
	}

}

void residualU(double* resU,double* U, double* P , int start_row, int stop_row , int N_row, int N_col, double mu, double dx, double P_const, int color  ){
	
	//Computes the residual of U and stores it in ResU. 
	//color codes 1-red 2-black  

	int k; //I have done this in a dumb way but wont fix it unless I have to 

	//set to random large numbers to ensure everything is set correctly 
	double uC=-420;
	double uE=666;
	double uW=420;
	double uN=-666;
	double uS=420;
	double pC=666;
	double pW=420;

	int N_col_p=N_col -1; 



//Red black is a pain in the ass because for even numbers of cols there are the same number of red and blacks in each row
// but for odd cols it alternates between n+1/2 and n-1 /2. I made two cases to accound for this.
	if(N_col%2 ==0){


		for(int j=start_row; j<stop_row; j++){
			for(int i=0; i<N_col/2; (i++) ){

			k=i; //K holds actual value of i while i is mapped to its red or black counter part
			i=2*i + (j+color)%2; 

			//printf("j=%d , index= %d\n", j, i );

			//Very happy with this organizational scheme here. You should have seen the mess I had for my first go the program. 

			//Center points of stencil 
				uC=U[i+N_col*j];
				pC=P[i+N_col_p*j]; 


				//Handle the x left boundary 
				if(i!=0 ){
					uW=U[i-1+N_col*j];
					pW=P[(i-1)+N_col_p*j];
				}

				//left points
				if(i==0 ){
					uW=uC;
					pW=2*P_const - pC;
				}

				//handle x right interior points
				if(i!=N_col-1){
					uE=U[i+1+N_col*j];
				}

				//handle x right boundary points
				if(i==N_col-1){
					uE=uC;
					pC=-P[i-1+N_col_p*j]; 
				}


				//Handle the up interior point
				if(j!=0 ){
					uN=U[i+N_col*(j-1)];
				}

				//Handle the y top boundary case
				if(j==0 ){
					uN=-uC;
				}

				//Handle the y down interior case
				if(j!=N_row-1 ){
					uS=U[i+N_col*(j+1)];
				}

				//Handle the y down boundary case
				if(j==N_row-1 ){
					uS=-uC;
				}


				// Compute Residual 

				resU[i+N_col*j]= mu*(uE + uW - 2*uC)+ mu*(uN + uS - 2*uC) - dx*(pC-pW);
							//printf("Val=%lf \n", resU[i+N_col*j] );
			
			//reset and increment
			i=k; 

			}
		}

	}



	if(N_col%2 ==1){



		for(int j=start_row; j<stop_row; j++){
			for(int i=0; i<(N_col+1 )/2 - (j+color)%2 ; (i++) ){

			k=i;


			i=2*i + (j+color)%2; 

						//Center points of stencel 
				uC=U[i+N_col*j];
				pC=P[i+N_col_p*j]; 


				//Handle the x left boundary 
				if(i!=0 ){
					uW=U[i-1+N_col*j];
					pW=P[(i-1)+N_col_p*j];
				}

				//left points
				if(i==0 ){
					uW=uC;
					pW=2*P_const - pC;
				}

				//handle x right interior points
				if(i!=N_col-1){
					uE=U[i+1+N_col*j];
				}

				//handle x right boundary points
				if(i==N_col-1){
					uE=uC;
					pC=-P[i-1+N_col_p*j]; 
				}


				//Handle the up interior point
				if(j!=0 ){
					uN=U[i+N_col*(j-1)];
				}

				//Handle the y top boundary case
				if(j==0 ){
					uN=-uC;
				}

				//Handle the y down interior case
				if(j!=N_row-1 ){
					uS=U[i+N_col*(j+1)];
				}

				//Handle the y down boundary case
				if(j==N_row-1 ){
					uS=-uC;
				}


				// Compute Residual 

				resU[i+N_col*j]= mu*(uE + uW - 2*uC)+ mu*(uN + uS - 2*uC) - dx*(pC-pW);

			
			//reset and increment
			i=k; 


			}
		}

	}

}



void residualV(double* resV,double* V,double* P , int start_row, int stop_row , int N_row, int N_col, double mu, double dx, double P_const, int color){
	


	//Computes the residual of V and stores it in ResV. 
	//color codes 1-red 2-black  


	int k;
	int N_col_p=N_col; //To accound for the MAC grid sizes
	double vC=-420;
	double vE=666;
	double vW=420;
	double vN=-666;
	double vS=420;
	double pC=666;
	double pN=420; 




if(N_col%2 ==0){

		for(int j=start_row; j<stop_row; j++){
			for(int i=0; i<N_col/2; (i++) ){

			k=i;
			i=2*i + (j+color)%2; 
			//Same idea is resU but with fewer comments

				vC=V[i+N_col*j];
				pC=P[i+N_col_p*j];


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


				if(j!=0 ){
					vN=V[i+N_col*(j-1)];
					pN=P[i+N_col_p*(j-1)];

				}

				if(j!=N_row-1 ){
					vS=V[i+N_col*(j+1)];
				}

				//If V on y boundary we force it to remain 0
				
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


				// Compute Residual 
				resV[i+N_col*j]= mu*(vE + vW - 2*vC)  + mu*(vN + vS - 2*vC) - dx*(pC-pN);


			i=k;


			}
		}
	}



if(N_col%2 ==1){


		for(int j=start_row; j<stop_row; j++){
			for(int i=0; i<((N_col+1)/2 - (j+color)%2); (i++) ){

			k=i;
			i=2*i + (j+color)%2; 
			//Same idea is resU but with fewer comments

				vC=V[i+N_col*j];
				pC=P[i+N_col_p*j];


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


				if(j!=0 ){
					vN=V[i+N_col*(j-1)];
					pN=P[i+N_col_p*(j-1)];

				}

				if(j!=N_row-1 ){
					vS=V[i+N_col*(j+1)];
				}

				//If V on y boundary we force it to remain 0
				
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


				// Compute Residual 
				resV[i+N_col*j]= mu*(vE + vW - 2*vC)  + mu*(vN + vS - 2*vC) - dx*(pC-pN);


			i=k;

			}
		
		}
	}




}


//residualP( resP, u, v, p , N+1, N-1, mu, dx, P_const, color);

void residualP(double* resP, double* U, double* V, double* P ,  int start_row, int stop_row , int N_row, int N_col, double mu, double dx, double P_const, int color){

	int k;
	double vC=-420;
	double vS=666;
	double uC=420;
	double uE=-666;


	int N_col_u=N_col+1; 


//The p update is dead simple.

// if N_col%2 ==0 then the number of red or black in each row does not change 
// if N_col%2 ==1 then the number of red or black alternates between 0 and 1. 	

	if(N_col%2 ==0){
		//printf("err even pressure");

		for(int j=start_row; j<stop_row; j++){

			for(int i=0; i<N_col/2; (i++) ){

			k=i;

			i=2*i + (j+color)%2; 


				vC=V[i+N_col*j];
				vS=V[i+N_col*(j+1)];

				uC=U[i+N_col_u*j];
				uE=U[i+1+N_col_u*j];


				// Compute Residual 

				resP[i+N_col*j]=-(uE-uC) - (vS-vC);


			i=k;
			}

		}

	}

	if(N_col%2 ==1){
		for(int j=start_row; j<stop_row; j++){

			for(int i=0; i<((N_col+1 )/2   ); (i++) ){

			k=i;

			i=2*i + (j+color)%2; 


				vC=V[i+N_col*j];
				vS=V[i+N_col*(j+1)];

				uC=U[i+N_col_u*j];
				uE=U[i+1+N_col_u*j];


				// Compute Residual 

				resP[i+N_col*j]=-(uE-uC) - (vS-vC);


			i=k;

			}
	

		}

	}



}



void updateRule(double* M, double* resM , int start_row, int stop_row ,int N_row , int N_col, double omega, int color){


// if N_col%2 ==0 then the number of red or black in each row does not change 
// if N_col%2 ==1 then the number of red or black alternates between 0 and 1. 	

int k;


	if(N_col%2 ==0){
		//Make this with a start and stop 
		for(int j=start_row; j<stop_row; j++){
			for(int i=0; i<(N_col/2) ; i++){
				k=i;

				i=2*i + (j+color)%2; 

				M[i+(N_col)*j]= M[i+(N_col)*j] + omega* resM[i+(N_col)*j];

				i=k;
			}
		}

	}


	if(N_col%2 ==1){
		//Make this with a start and stop 
		for(int j=start_row; j<stop_row; j++){
			for(int i=0; i<((N_col+1)/2 - (j+color)%2) ; i++){

				k=i;
				//printf("update rule j= %d\n", j );

				i=2*i + (j+color)%2;  //j+color %2 for red will go 0 1 0 1 and for black will go 1 0 1 0 1
				//I feel really proud of this indexing. 

				M[i+(N_col)*j]= M[i+(N_col)*j] + omega* resM[i+(N_col)*j];

				i=k;
			}
		}

	}

}













