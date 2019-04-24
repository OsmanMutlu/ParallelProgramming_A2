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

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem

 double stats(double **E, int m, int n, double *_mx){
     double mx = -1;
     double l2norm = 0;
     int i, j;
     for (j=0; j<m; j++)
       for (i=0; i<n; i++) {
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

  MPI_Status status[4*P];
  MPI_Request request[4*P];
  for (int i=0; i<4*P; i++)
    request[i] = MPI_REQUEST_NULL;

  double T=1000.0;
  int m=200,n=200;
  int plot_freq = 0;
  int px = 1, py = 1;
  int no_comm = 0;
  int num_threads=1;

  double **E, **R, **E_prev;

  if (myrank == 0) {
    cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);

    if (P != px*py) {
      cout << "Your px*py must be equal to the total amount of processors" << endl;
      return(-1);
    }

    if ((n % px != 0) || (n % py != 0)) {
      cout << px << py << n << m << endl;
      cout << "Please make sure that your n is divisible by px and your m is divisible by py" << endl;
      return(-1);
    }

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
  int grid_size_x = n/px;
  int grid_size_y = m/py;

  double **my_E, **my_R, **my_E_prev, *ghost1, *ghost2, ghost3[grid_size_y], ghost4[grid_size_y], send3[grid_size_y], send4[grid_size_y];

  ghost1 = (double*)malloc((grid_size_x+2) * sizeof(double));
  ghost2 = (double*)malloc((grid_size_x+2) * sizeof(double));

  my_E = alloc2D(grid_size_y+2,grid_size_x+2);
  my_E_prev = alloc2D(grid_size_y+2,grid_size_x+2);
  my_R = alloc2D(grid_size_y+2,grid_size_x+2);

  // We do this for every process to avoid communication cost.
  // Other option would be initializing these for only master and then scattering them before while loop began.
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

  // Assign each my_E
  for (int a = 0; a<py; a++) {
    for (int b = 0; b<px; b++) {
      for (int i = 0; i<grid_size_y+2; i++) {
	for (int j = 0; j<grid_size_x+2; j++) {
	  if (myrank == a*px + b) {
	    my_E[i][j] = E[a*grid_size_y + i][b*grid_size_x + j];
	    my_E_prev[i][j] = E_prev[a*grid_size_y + i][b*grid_size_x + j];
	    my_R[i][j] = R[a*grid_size_y + i][b*grid_size_x + j];
	  }
	}
      }
    }
  }

  double dx = 1.0/n;

  // For time integration, these values shouldn't change
  double rp= kk*(b+1)*(b+1)/4;
  double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
  double dtr=1/(epsilon+((M1/M2)*rp));
  double dt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
  double alpha = d*dt/(dx*dx);

  // Start the timer
  double t0;

  if (myrank == 0) {
    t0 = getTime();
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

    // North neighbour -> ghost1
    if (myrank / px == 0) { // They are in first row
      // Mirror upper bound
      for (int i=1; i<=grid_size_x; i++)
	ghost1[i] = my_E_prev[2][i];
    }
    else {
      // We can send like this because it is contigious
      MPI_Irecv(ghost1, grid_size_x+2, MPI_DOUBLE, (myrank / px - 1)*px + (myrank % px), 1, MPI_COMM_WORLD, &request[4*myrank + 0]);
      MPI_Send(my_E_prev[1], grid_size_x+2, MPI_DOUBLE, (myrank / px - 1)*px + (myrank % px), 1, MPI_COMM_WORLD);
    }

    // South neighbour -> ghost2
    if (myrank / px == py - 1) { // They are in last row
      // Mirror lower bound
      for (int i=1; i<=grid_size_x; i++)
        ghost2[i] = my_E_prev[grid_size_y-1][i];
    }
    else {
      // We can send like this because it is contigious
      MPI_Irecv(ghost2, grid_size_x+2, MPI_DOUBLE, (myrank / px + 1)*px + (myrank % px), 1, MPI_COMM_WORLD, &request[4*myrank + 1]);
      MPI_Send(my_E_prev[grid_size_y], grid_size_x+2, MPI_DOUBLE, (myrank / px + 1)*px + (myrank % px), 1, MPI_COMM_WORLD);
    }

    // Left neighbour -> ghost3
    if (myrank % px == 0) { // They are in the leftmost column
      // Mirror leftmost bound
      for (int j=1; j<=grid_size_y; j++)
        ghost3[j-1] = my_E_prev[j][2];
    }
    else {
      MPI_Irecv(ghost3, grid_size_y, MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &request[4*myrank + 2]);
      for (int i=1; i<=grid_size_y; i++)
    	send3[i-1] = my_E_prev[i][1];
      MPI_Send(send3, grid_size_y, MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD);
    }

    // Right neighbour -> ghost4
    if (myrank % px == px - 1) { // They are in the rightmost column
      // Mirror rightmost bound
      for (int j=1; j<=grid_size_y; j++)
    	ghost4[j-1] = my_E_prev[j][grid_size_x-1];
    }
    else {
      MPI_Irecv(ghost4, grid_size_y, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &request[4*myrank + 3]);
      for (int i=1; i<=grid_size_y; i++)
    	send4[i-1] = my_E_prev[i][grid_size_x];
      MPI_Send(send4, grid_size_y, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD);
    }

    MPI_Waitall(4*P, request, status);

    for (int i=1; i<=grid_size_x; i++) { // North
      my_E_prev[0][i] = ghost1[i];
    }
    for (int i=1; i<=grid_size_x; i++) { // South
      my_E_prev[grid_size_y+1][i] = ghost2[i];
    }

    for (int i=1; i<=grid_size_y; i++) { // Left
      my_E_prev[i][0] = ghost3[i-1];
    }
    for (int i=1; i<=grid_size_y; i++) { // Right
      my_E_prev[i][grid_size_x+1] = ghost4[i-1];
    }

    simulate(my_E, my_E_prev, my_R, alpha, grid_size_x, grid_size_y, kk, dt, a, epsilon, M1, M2, b);

    double **tmp = my_E; my_E = my_E_prev; my_E_prev = tmp;

    MPI_Barrier(MPI_COMM_WORLD);

  }//end of while loop

  double **E2,**my_E2;
  E2 = alloc2D(m,n); // only allocating this for rank 0 causes segmentation fault!

  if (P != 1) {

    my_E2 = alloc2D(grid_size_y, grid_size_x);
    
    for (int i=1; i<=grid_size_y; i++)
      for (int j=1; j<=grid_size_x; j++)
    	my_E2[i-1][j-1] = my_E_prev[i][j];
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&my_E2[0][0], grid_size_y*grid_size_x, MPI_DOUBLE, &E2[0][0], grid_size_y*grid_size_x, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  }
  else {
    for (int i=1; i<=m; i++)
      for (int j=1; j<=n; j++)
  	E2[i-1][j-1] = my_E_prev[i][j];
  }

  if (myrank == 0) {

    double time_elapsed = getTime() - t0;

    double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed ;
    double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;

    cout << "Number of Iterations        : " << niter << endl;
    cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
    cout << "Sustained Gflops Rate       : " << Gflops << endl; 
    cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 

    double mx;
    double l2norm = stats(E2,m,n,&mx);
    cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;
  
  }

  MPI_Barrier(MPI_COMM_WORLD); // To make sure that processors other than rank 0 don't come here and start freeing stuff up before we calculate l2norm
  free (E2);
  free (my_E2);
  free (E);
  free (E_prev);
  //free (R); // Throws an error for some reason
  free (my_E);
  free (my_E_prev);
  free (my_R);
  free (ghost1);
  free (ghost2);
  MPI_Finalize();
  
  return 0;
}
