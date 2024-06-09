//Library------------------------------------------------------------------------------------------
#include <iostream>	// cout,cin
#include <cmath>	// math tools
#include <fstream>	// make file
#include <cstdio>       // printf, NULL 
#include <cstdlib>      // srand, rand 
#include <ctime>        // time 
using namespace std;
//Random-------------------------------------------------------------------------------------------
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum) {
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
long iseed = -36;
void Randomize() {
  iseed = -time(NULL);  
}
inline double Random() { return ran2(&iseed); }
inline int Random(int N) { return int(ran2(&iseed)*N); }
//END_RANDOM
//-------------------------------------------------------------------------------------------------


const int N = 64; 			// lenght of lattice
const int LMask = N-1;			// for boundary condition
int s[N][N];	 		        // Spins
int MCstep;				// Monte carlo step
double J = 1.0;				// ferromagnetic factor
double h = 0.0; 			//External Field
double T = 4.0;				//Temperature
//INITIAL------------------------------------------------------------------------------------------
void initial(){		// in this function we put initial value of spin lattice 
        double w=0;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
                        w = Random();
			s[i][j] = ((w>=0 && w<=0.5)?-1:1);
                  //      s[i][j] = 1;
		}
	} 
} 
//-------------------------------------------------------------------------------------------------

// Rotate i if cross the boundary vs the periodic boundary condition
inline int Period(int i) {
    return i & LMask;
}

//Calculate Delta E if s(i,j) flipped
double Delta_E_Flip(int i, int j) {
    float sum = 0.0;
    sum = s[i][Period(j+1)] + s[i][Period(j-1)] + s[Period(i-1)][j] + s[Period(i+1)][j];
    return 2*s[i][j]* (h + J * sum);
}

void Single_monte_carlo_step(){
	for(int c=0;c<N*N;c++){	
		int i = N * Random();
		int j = N * Random();
		double dE = Delta_E_Flip(i,j);	
	 	if ( dE < 0 ) {    // If the total energy decreases after spin-flip, accept the infinitesimal step
         		   s[i][j] = -s[i][j];
       		} else if ( Random() < exp(-dE/T) ) {
           					 s[i][j] = -s[i][j];
        	}
	}
}
//-------------------------------------------------------------------------------------------------
int main(){
    
ofstream out("absdata.txt");		// make txt file 
ofstream fout("nonabsdata.txt");	// make txt file

	system("clear");		// clear screen of the terminal
        cout<<"------------------------------------------\n";	
        cout<<"|Ising Model,Provided by Amir Masominia  |\n";
        cout<<"------------------------------------------\n";
        cout<<"|Please Insert Value of Temperature:      \t";
        cin>>T;
        cout<<"------------------------------------------\n";
        cout<<"|Please Insert Number of Monte Carlo Steps:\t";
	cin>> MCstep;
        cout<<"------------------------------------------\n";

	initial(); // Putting value of spin +1 or -1
	
        
        for(int i=0;i<MCstep;i++){
			Single_monte_carlo_step();
                        if(i%1000 == 0 ){ system("clear");
					  cout<<"-----------------------------------------------\n";
					  cout<<"|Round\t"<<((i/1000)+1)<<" is running. Temprature: "<<T<<".          |\n";
                                          cout<<"-----------------------------------------------\n";
                        }	// this is just for knowledge of user 
			double sum=0.0;
			for(int l=0;l<N;l++){
				for(int m=0;m<N;m++){
					sum = sum + s[l][m];
				}
			} // calculate of the order parameter
			fout<< i << '\t' << (sum/(N*N)) << '\n'; // (sum/(N*N)) >> normalize of the order parameter
			out<< i << '\t' << abs(sum/(N*N)) << '\n';
	}
      
	fout.close();	//close file
        out.close();	//close file
	
	system("clear");
        cout<<"------------------------------------------\n";
        cout<<"Program successfully finished.           |\n";
        cout<<"------------------------------------------\n";
        return 0;
}  
