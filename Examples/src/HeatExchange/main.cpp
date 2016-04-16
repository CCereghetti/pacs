#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include <string>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
/*SPIEGAZIONE: 
- per la Challenge 1.2 runnare il file run.sh, che contiene due file con i parametri di imput differenti: variano
sia il nome di uscita del file, che la norma da utilizzare. In questo modo posso confrontare il numero di iterazioni
*/

/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  const double& toler=param.toler;   // Tolerance for stopping criterion
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  const auto& outputname=param.outputname; //Name of the output file
  const auto& outputwhere=param.outputwhere; //Where to print the results (0=screen, 1=file, 2=both)
  const auto& norm=param.norm;	//Say witch norm to use: (0: L², 1: H¹)
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
 /*--------------------------------------
  
  ////CHALLENGE 1.3////
  //costruisco i tre vettori a,b,c (che sono le 3 diagonali della matrice), che userò per applicare Thomas
  //a è la diagonale, b la sottodiagonale, c la sovradiagonale
  vector<double> a(M,2+h*h*act), b(M-1,-1), c(M-1,-1);
  vector<double> alfa(M), beta(M-1);
  alfa[0]=a[0];
  for (int i=1; i<M,i++) {
  	beta[i-1]=b[i-1]/alfa[i-1];
	alfa[i]=a[i]-beta[i]*c[i-1];  
  }
  //costruisco il vettore della soluzione b
  vector<double> f(M,0.);
  f[0]=To-Te;
  //Risolvo Ly=f
  vector<double> y(M);
  y[0]=f[0];
  for (int i=1; i<M; i++)
  	y[i]=f[i]-beta[i]*y[i-1];
  	
  //Risolvo Ux=y
  vector<double> x(M);
  x[M-1]=y[M-1]/alfa[M-1];
  for(int i=M-1; i>=0; i--)
  	x[i]=(y[i])
  	
 -------------------------------------------------*/
 
  // Solution vector
  vector<double> theta(M+1);
  vector<double> dif(M+1);			//vettore contenente la differenza tra x_new-x_old
  // Gauss Siedel is initialised with a linear variation
  // of T
  
  for(unsigned int m=0;m <= M;++m) {
     theta[m]=(1.-m*h)*(To-Te)/Te;
     dif[m]=0;						//inizializzo il vettore delle differenze a zero
     }
  
  // Gauss-Seidel
  // epsilon=||x^{k+1}-x^{k}||
  // Stopping criteria epsilon<=toler
  
  int iter=0;
  double xnew, epsilon;
  cout<<"I compute with the "<<norm<<" norm"<<endl;  
     do
       { epsilon=0.;

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   dif[m]=xnew-theta[m];
	   theta[m] = xnew;
         }
     
     //compute the L² norm
      epsilon+=h/2*(dif[1]*dif[1]);
      for (int m=2; m<M; m++) {
      	epsilon+=h/2*(dif[m]*dif[m]+dif[m-1]*dif[m-1]);     //computo l'integrale col metodo dei trapezi
     }
	 //Last row
	 xnew = theta[M-1]; 
	 dif[M]=xnew-theta[M];
	 theta[M]=  xnew;
	 epsilon+=h/2*(dif[M]*dif[M]+dif[M-1]*dif[M-1]);
	
	//Implement the H¹ norm, if needed
	if(norm==1) {
	epsilon+=1/(2*h)*(dif[1]*dif[1]);
	for(int m=2; m<M+1; m++) {
		epsilon+=1/(2*h)*(dif[m]*dif[m]+dif[m-1]*dif[m-1]);
	}
	}
	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }

 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 
	 ofstream f;
     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);
	
	switch(outputwhere){
	case 0:
	cout<<"Plot on the screen"<<endl;
	    for(int m = 0; m<= M; m++)
       {
         	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
       }
        // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<endl;
	 break;
	 case 1:
     	cout<<"Wiew result in a file named: "<<outputname<<endl;
     	f.open(outputname);
     	for(int m = 0; m<= M; m++)
       {
	 	// \t writes a tab 
         f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
       }
       f.close();
       break;
   case 2:  
     cout<<"Wiew result in a file named: "<<outputname<<" and screen"<<endl;
     f.open(outputname);
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
         	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
       }
      // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<endl;
       f.close();
       break;
       }
     return status;
}
