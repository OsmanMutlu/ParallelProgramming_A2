/* 
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory 
 * and reimplementation by Scott B. Baden, UCSD
 * 
 * Modified and  restructured by Didem Unat, Koc University
 *
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
using namespace std;


// Utilities
// 

// Timer
// Make successive calls and take a difference to get the elapsed time.
static const double kMicro = 1.0e-6;
double getTime()
{
    struct timeval TV;
    struct timezone TZ;

    const int RC = gettimeofday(&TV, &TZ);
    if(RC == -1) {
            cerr << "ERROR: Bad call to gettimeofday" << endl;
            return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()

// Allocate a 2D array
double **alloc2D(int m,int n){
   double **E;
   int nx=n, ny=m;
   E = (double**)malloc(sizeof(double*)*ny + sizeof(double)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++)
     E[j] = (double*)(E+ny) + j*nx;
   return(E);
}

// double *to1D(int m, int n, double **Mat)
// {
//   double *array[m * n];
//   for (int i=0; i<m; i++) {
//     for (int j=0; j<n; j++) {
//       array[m * i + j] = Mat[i][j];
//     }
//   }
//   return(array);
// }

// double **to2D(int m, int n, double *array)
// {
//   double **Mat = alloc2D(m, n);
//   for (int i=0; i<m; i++) {
//     for (int j=0; j<n; j++) {
//       Mat[i][j] = array[m * i + j];
//     }
//   }
//   return(Mat);
// }

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
 double stats(double **E, int m, int n, double *_mx){
     double mx = -1;
     double l2norm = 0;
     int i, j;
     for (j=1; j<=m; j++)
       for (i=1; i<=n; i++) {
	   l2norm += E[j][i]*E[j][i];
	   if (E[j][i] > mx)
	       mx = E[j][i];
      }
     *_mx = mx;
     l2norm /= (double) ((m)*(n));
     l2norm = sqrt(l2norm);
     return l2norm;
 }

// External functions
extern "C" {
    void splot(double **E, double T, int niter, int m, int n);
}
void cmdLine(int argc, char *argv[], double& T, int& n, int& px, int& py, int& plot_freq, int& no_comm, int&num_threads);


void simulate (double** E,  double** E_prev,double** R,
	       const double alpha, const int n, const int m, const double kk,
	       const double dt, const double a, const double epsilon,
	       const double M1,const double  M2, const double b)
{
  int i, j; 
    /* 
     * Copy data from boundary of the computational box 
     * to the padding region, set up for differencing
     * on the boundary of the computational box
     * Using mirror boundaries
     */

    // Solve for the excitation, the PDE
    for (j=1; j<=m; j++){
      for (i=1; i<=n; i++) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
      }
    }
    
    /* 
     * Solve the ODE, advancing excitation and recovery to the
     *     next timtestep
     */
    for (j=1; j<=m; j++){
      for (i=1; i<=n; i++)
	E[j][i] = E[j][i] -dt*(kk* E[j][i]*(E[j][i] - a)*(E[j][i]-1)+ E[j][i] *R[j][i]);
    }
    
    for (j=1; j<=m; j++){
      for (i=1; i<=n; i++)
	R[j][i] = R[j][i] + dt*(epsilon+M1* R[j][i]/( E[j][i]+M2))*(-R[j][i]-kk* E[j][i]*(E[j][i]-b-1));
    }
    
}

// Main program
int main (int argc, char** argv)
{
  /*
   *  Solution arrays
   *   E is the "Excitation" variable, a voltage
   *   R is the "Recovery" variable
   *   E_prev is the Excitation variable for the previous timestep,
   *      and is used in time integration
   */
  int P, myrank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // status and request will probably be arrays
  MPI_Status status;
  MPI_Request request[2];

  double T=1000.0;
  int m=200,n=200;
  int plot_freq = 0;
  int px = 1, py = 1;
  int no_comm = 0;
  int num_threads=1;

  double **E, **R, **E_prev;

  if(myrank == 0) {
    cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);
  }
  // Various constants - these definitions shouldn't change
  const double a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;

  int sendbuff[6];
  
  sendbuff[0] = n;
  sendbuff[1] = plot_freq;
  sendbuff[2] = px;
  sendbuff[3] = py;
  sendbuff[4] = no_comm;
  sendbuff[5] = num_threads;

  MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(sendbuff, 6, MPI_INT, 0, MPI_COMM_WORLD);
  n = sendbuff[0];
  plot_freq = sendbuff[1];
  px = sendbuff[2];
  py = sendbuff[3];
  no_comm = sendbuff[4];
  num_threads = sendbuff[5];
  m = n;
  // Allocate contiguous memory for solution arrays
  // The computational box is defined on [1:m+1,1:n+1]
  // We pad the arrays in order to facilitate differencing on the 
  // boundaries of the computation box
  int grid_size = n/P;

  // int sizes[2]    = {m+2, n+2}; // gA size
  // int subsizes[2] = {grid_size+2, n+2}; // lA size
  // int starts[2]   = {0,0}; // where this one starts
  // MPI_Datatype type, subarrtype;
  // MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
  // MPI_Type_create_resized(type, 0, (grid_size+2) * sizeof(double), &subarrtype);
  // MPI_Type_commit(&subarrtype);

  double **my_E, **my_R, **my_E_prev, *ghost1, *ghost2;

  my_E = alloc2D(grid_size+2,n+2);
  my_E_prev = alloc2D(grid_size+2,n+2);
  my_R = alloc2D(grid_size+2,n+2);
  
  // For 2D geometry
  // my_E = alloc2D(grid_size,grid_size);
  // my_E_prev = alloc2D(grid_size+2,grid_size+2);
  // my_R = alloc2D(grid_size,grid_size);

  // We do this for every process to avoid communication cost.
  // Other option would be initializing these for only master and then scattering them before while loop began.
  // if(myrank == 0) {
  E = alloc2D(m+2,n+2);
  E_prev = alloc2D(m+2,n+2);
  R = alloc2D(m+2,n+2);

  int i,j;
  // Initialization
  for (j=1; j<=m; j++)
    for (i=1; i<=n; i++)
      E_prev[j][i] = R[j][i] = 0;
  
  for (j=1; j<=m; j++)
    for (i=n/2+1; i<=n; i++)
      E_prev[j][i] = 1.0;

  for (j=m/2+1; j<=m; j++)
    for (i=1; i<=n; i++)
      R[j][i] = 1.0;
  // }

  for (int p = 0; p<P; p++) {
    for (int i = 0; i<grid_size+2; i++) {
      for (int j = 0; j<n+2; j++) {
	if (myrank == p) {
	  my_E[p*grid_size + i][j] = E[p*grid_size + i][j];
	  my_E_prev[p*grid_size + i][j] = E_prev[p*grid_size + i][j];
	  my_R[p*grid_size + i][j] = R[p*grid_size + i][j];
	}
      }
    }
  }

  // if it counts one pointer as n+2
  int sendcounts[P];
  int displs[P];
  for (int i = 0; i < P; i++) {
    sendcounts[i] = (grid_size+2) * (n+2);
  }
  int disp = 0;
  for (int i=0; i < P; i++) {
    displs[i] = disp;
    disp += grid_size;
  }
  
  double dx = 1.0/n;

  // For time integration, these values shouldn't change
  double rp= kk*(b+1)*(b+1)/4;
  double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
  double dtr=1/(epsilon+((M1/M2)*rp));
  double dt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
  double alpha = d*dt/(dx*dx);

  // Start the timer
  double t0 = getTime();

  if (myrank == 0) {
    cout << "Grid Size       : " << n << endl;
    cout << "Duration of Sim : " << T << endl;
    cout << "Time step dt    : " << dt << endl;
    cout << "Process geometry: " << px << " x " << py << endl;
    if (no_comm)
      cout << "Communication   : DISABLED" << endl;
      
    cout << endl;
  }

  // Simulated time is different from the integer timestep number
  // Simulated time
  double t = 0.0;
  // Integer timestep number
  int niter=0;

  while (t<T) {
    
    t += dt;
    niter++;

    simulate(my_E, my_E_prev, my_R, alpha, n, m, kk, dt, a, epsilon, M1, M2, b);

    // Always tag = 1 for receiving and tag = 2 for sending
    // ghost1 is the message coming from north, and ghost2 is the message coming from south
    if (P == 1) { // If there is only one process
      // Mirror upper bound
      for (i=1; i<=n; i++)
	my_E[0][i] = my_E[2][i];
      // Mirror lower bound
      for (i=1; i<=n; i++)
	E[m+1][i] = E[m-1][i];
    }
    else if (myrank == 0) {
      // Mirror upper bound
      for (i=1; i<=n; i++)
	my_E[0][i] = my_E[2][i];
      // Recieve from south neighbour
      MPI_Irecv(ghost2, n+2, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, request[0]); // REQUEST????
      MPI_Send(my_E[grid_size], n+2, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
      my_E[grid_size+1] = ghost2;
    }
    else if (myrank == P-1) {
      // Mirror lower bound
      for (i=1; i<=n; i++)
	E[m+1][i] = E[m-1][i];
      // Recieve from north neighbour
      MPI_Irecv(ghost1, n+2, MPI_DOUBLE, P-2, 1, MPI_COMM_WORLD, request[0]); // REQUEST????
      MPI_Send(my_E[1], n+2, MPI_DOUBLE, P-2, 2, MPI_COMM_WORLD);
      my_E[0] = ghost1;
    }
    else {
      MPI_Irecv(ghost1, n+2, MPI_DOUBLE, myrank-1, 1, MPI_COMM_WORLD, request[0]);
      MPI_Irecv(ghost2, n+2, MPI_DOUBLE, myrank+1, 1, MPI_COMM_WORLD, request[1]);
      MPI_Send(my_E[1], n+2, MPI_DOUBLE, myrank-1, 2, MPI_COMM_WORLD);
      MPI_Send(my_E[grid_size], n+2, MPI_DOUBLE, myrank+1, 2, MPI_COMM_WORLD);
      my_E[grid_size+1] = ghost2;
      my_E[0] = ghost1;
    }

    MPI_Waitall(2, request, &status);

    // Mirror leftmost bound
    for (j=1; j<=m; j++)
      my_E[j][0] = my_E[j][2];
    // Mirror rightmost bound
    for (j=1; j<=m; j++)
      my_E[j][n+1] = my_E[j][n-1];

    //Update E_prev
    my_E_prev = my_E;

    if (plot_freq){
      MPI_Gatherv(my_E, (grid_size+2) * (n+2), MPI_DOUBLE,
      		  E, sendcounts, displs, MPI_DOUBLE,
      		  0, MPI_COMM_WORLD);

      int k = (int)(t/plot_freq);
      if ((t - k * plot_freq) < dt){
    	splot(E,t,niter,m+2,n+2);
      }
    }

  }//end of while loop

  MPI_Gatherv(my_E, (grid_size+2) * (n+2), MPI_DOUBLE,
	      E, sendcounts, displs, MPI_DOUBLE,
	      0, MPI_COMM_WORLD);

  if(myrank == 0) {
    double time_elapsed = getTime() - t0;

    double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed ;
    double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;

    cout << "Number of Iterations        : " << niter << endl;
    cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
    cout << "Sustained Gflops Rate       : " << Gflops << endl; 
    cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 

    double mx;
    double l2norm = stats(E,m,n,&mx);
    cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;

    if (plot_freq){
      cout << "\n\nEnter any input to close the program and the plot..." << endl;
      getchar();
    }
  
    free (E);
    free (E_prev);
    free (R);
  }

  free (my_E);
  free (my_E_prev);
  free (my_R);
  free (ghost1);
  free (ghost2);
  free (my_R);
  MPI_Finalize();
  
  return 0;
}
