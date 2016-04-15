#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include <string>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
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

  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  
  ////CHALLENGE 1.3////
  
  //costruisco i tre vettori a,b,c (che sono le 3 diagonali della matrice), che userò per applicare Thomas
  //a è la diagonale, b la sottodiagonale, c la sovradiagonale
  vector<double> a(M,2.+h*h*act), b(M-1,-1.), c(M-1,-1.);
  vector<double> alfa(M), beta(M-1);
  a[M-1]=1;
  alfa[0]=a[0];
  for (int i=1; i<M; i++) {
  	beta[i-1]=b[i-1]/alfa[i-1];
	alfa[i]=a[i]-beta[i-1]*c[i-1];  
  }
  //costruisco il vettore della soluzione b
  vector<double> f(M);
  f[0]=To-Te;
  //Risolvo Ly=f
  vector<double> y(M);
  y[0]=f[0];
  for (int i=1; i<M; i++)
  	y[i]=f[i]-beta[i-1]*y[i-1];
  	
  //Risolvo U*theta=y
  vector<double> theta(M+1);
  theta[M]=y[M-1]/alfa[M-1];		//il vettore theta sarà il mio vettore soluzione
  theta[0]=To-Te;		
  for(int i=M-1; i>0; i--)
  	theta[i]=(y[i-1]-c[i-1]*theta[i+1])/alfa[i-1];
  	

  // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 
	 ofstream ff;
     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);
	
	switch(outputwhere){
	case 0:
	cout<<"Plot on the screen"<<endl;
	    for(int m = 0; m<= M; m++)
       {
         	 // An example of use of tie and tuples
	 		std::tie(coor[m],sol[m],exact[m])=
	  		std::make_tuple(m*h*L,Te+theta[m],thetaa[m]);
       }
        // Using temporary files (another nice use of tie)
     	gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       	"w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       	"w l title 'uex'"<<endl;
	break;
	case 1:
     	cout<<"Wiew result in a file named: "<<outputname<<endl;
     	ff.open(outputname);
     	for(int m = 0; m<= M; m++)
       {
	 	// \t writes a tab 
         ff<<m*h*L<<"\t"<<Te+theta[m]<<"\t"<<thetaa[m]<<endl;
       }
       ff.close();
       break;
   case 2:  
     cout<<"Wiew result in a file named: "<<outputname<<" and screen"<<endl;
     ff.open(outputname);
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         ff<<m*h*L<<"\t"<<Te+theta[m]<<"\t"<<thetaa[m]<<endl;
         	 // An example of use of tie and tuples!
         
		std::tie(coor[m],sol[m],exact[m])=
	    std::make_tuple(m*h*L,Te+theta[m],thetaa[m]);
       }
      // Using temporary files (another nice use of tie)
    	gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<endl;
       ff.close();
   	break;
   }
   
 
     return status;
}
