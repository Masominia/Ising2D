
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>      /* printf, NULL */
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
using namespace std;
//-------------------------------------------------------------------------------------------------
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


const int N = 64; 
const int LMask  = N-1;
const int MCstep = 15000;
int s[N][N];	 		        // Spins
double E; 	                 	// Energy
double J  = 1.0;
double h  = 0.0;                        // External Field
double T  = 0.7;                        // Temperature
double M  = 0.0;			// Magnetization-Variable
double F_T= 4.0;			// Final Value of Temperature for stopping Temp loop
double I_T= 0.001;			// Interval value of Temp
const int value=3300;
double V [MCstep][value];// I'll put value of mangnezitation of each temperature in this array.
string name;
string png_filename;
char z,x,v,b;
int q=48,w=48,e=65,a=65;
//INITIAL------------------------------------------------------------------------------------------
void initial(){				// In This loop we put each spin in up mode.
        double w=0;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
                        s[i][j] = 1;
		}
	}
}
//-------------------------------------------------------------------------------------------------

// Rotate i if cross the boundary vs the periodic boundary condition
inline int Period(int i) {
    return i & LMask;
}

//Calculate Delta E if s(i,j) flipped
float Delta_E_Flip(int i, int j) {
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
       		}else if ( Random() < exp(-dE/T) ) {
           					 s[i][j] = -s[i][j];
        	}
	}
}
//-------------------------------------------------------------------------------------------------
int main(){
   	 
	ofstream out("Magnezitation.txt");	// I define a txt file for storing value of magnezitation
	ofstream fout("Var-Magnezitation.txt");// I define a txt file for storing value of Variance.
	ofstream put;

	initial(); 			// Putting value of spin +1 or -1
	int count = 0;		  	// this counter because of changing column in V.


	while(T<=F_T){
            system("clear");
            cout<<"------------------------------------------\n";
            cout<<"|Ising Model,Provided by amir MASOMINIA  |\n";
            cout<<"------------------------------------------\n";
	    cout<<"|Final Value of Temperature is: "<<F_T<<".\n";
	    cout<<"|Inteval Value of Temperature : "<<I_T<<".\n";
	    cout<<"------------------------------------------\n";
	    cout<<"|Temperature is :    "<<T<<".                |\n";
	    cout<<"------------------------------------------\n"; 
 	    cout<<"Value of step temperature :  "<< value <<"          |\n";
            cout<<"------------------------------------------\n";
//-------------------------------------------------------------------------------------------------
	if(q>57){w++; q = 48;}
	if(w>57){e++; w = 48;}
	if(e>90){a++; e = 65;}//change 	
	if(a>90){cout<< "Counter of file is over\n";	break;}
	name = char(a);
	name+= char(e);
	name+= char(w);
	name+= char(q);
	
	
	png_filename = char(a);
	png_filename+= char(e);
	png_filename+= char(w);
	png_filename+= char(q);
	

	png_filename += ".png";	
	name += ".txt";
	put.open(name.c_str(),ios::out);
	put << "set term png\n";
  	put << "set output \"" << png_filename << "\"\n";
  	put << "set xrange [ 0 : " << N << " ]\n";
  	put << "set yrange [ 0 : " << N << " ]\n";
  	put << "set nokey\n";
  	put << "set title \" Ising Model, Temperature : "<< T << " \"  \n";
  	put << "unset tics\n";	
	
	for(int l=0;l<N;l++){
		for(int m=0;m<N;m++){
			(s[l][m]<0? 
		put << "set object rectangle from " << m << "," << l << " to " << m+1 << "," << l+1 << " fc rgb 'blue'\n"
		:
		put << "set object rectangle from " << m << "," << l << " to " << m+1 << "," << l+1 << " fc rgb 'red'\n"
		); 	
			
		}
	}
	put << "plot 1\n";
  	put << "quit\n";
	put.close();
	q++;
	string wname;
	wname = "gnuplot " + name;
	system(wname.c_str());
	
//-------------------------------------------------------------------------------------------------	
	    for(int i=0;i<MCstep;i++){

                            Single_monte_carlo_step(); 					// change mode of each spin.


                            if(i%1000 == 0 ){ cout<<"|Round\t"<<((i/1000)+1)<<"\tis running.              |\n";} // debugging
                            double sum=0.0;						// Put zero value for next loop.
                         
			    for(int l=0;l<N;l++){
                                    for(int m=0;m<N;m++){
                                            sum = sum + s[l][m];
                                    } 					// End of <for> loop, counter m
                            }						// End of <for> loop, counter l 

			    
			    V[i][count] =   (sum+0.0)/(N*N);		// Now I putt value of magnezitation in V;
                            M = M + (sum/(N*N));			// I have summation of all value of the magnezitation for finding mean mag.


            } 								//End of <For> loop


            M = (M+0.0) / MCstep;					// Calculating mean Magnezitation of each single metropolice loop.


	    // for calculating Variance of M, first we must find difference M and V. then calculate squre of this difference.
	    for(int i=0;i<MCstep;i++) V[i][count] = pow((M - V[i][count]),2.0); 


	    // Now we have squre of diiference of the magnezitation of each temperature, now we must find average of the squre of difference.
	    double sum = 0.0;  // sum is a variable for putting summation of V. I put zero for first time because of other rounds.
	    for(int i=0;i<MCstep;i++) sum = sum + V[i][count];
	    fout<< T << '\t' << ((sum+0.0)/MCstep) << '\n'; 		// Finally I store these numbers in a file.
           

	    out<< T <<'\t' << M <<'\n'; // print in code.txt file for plotting Magnetization-Temperature Diagram. 


            M = 0.0;							// putting zero value in M for next Metropolice loop.
            T = T + I_T;						// Update Temperature for next while loop.

	    count++;							// update counter of column of the V for next loop
        }								// End of <while> loop

	
	out.close();
	fout.close();
//-------------------------------------------------------------------------------------------------

	system("clear");
	cout<<"------------------------------------------\n";
        cout<<"Now, Program is cleaning .txt file.      |\n";
        cout<<"------------------------------------------\n";

	q=48,w=48,e=65,a=65;
		
	for(int i=0; i<value; i++){
		if(q>57){w++; q = 48;}
		if(w>57){e++; w = 48;}
		if(e>90){a++; e = 65;}	
		if(a>90){cout<<"The counter of file is over\n";	break;}	
	
		name = char(a);	
		name+= char(e);
		name+= char(w);
		name+= char(q);
		name+=  ".txt";
		q++;
		name = "rm -r -f " + name;
		system(name.c_str());
		
	}

//-------------------------------------------------------------------------------------------------
	system("clear");
	cout<<"------------------------------------------\n";
        cout<<"Program successfully finished.           |\n";
        cout<<"------------------------------------------\n";
	
        return 0;
}  
