#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 5277
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-8
#define RNMX (1.0-EPS)

# include <stdio.h>
# include <iostream>
# include <string>
# include <math.h>

# include <vector>
# include <iostream>
# include <fstream>

using namespace std;

class RandomWrap
{
  
private:
	long idum;
	int j;
	long k;
	long idum2;
	long iy;
	long iv[NTAB];
	double temp;
	bool seeded;
  
public:

	RandomWrap()
	{
		seeded = false;
	}

	void set_seed(long seed)
	{
		if (!seeded)
		{
			idum2 = 123456789;
			iy = 0;
			idum = seed;
			if (idum < 1) idum=1; 
			//Be sure to prevent idum = 0.
			idum2=(idum);
			for (j=NTAB+7;j>=0;j--) 
			{		
				//Load the shuffle table (after 8 warm-ups).
				k=(idum)/IQ1;
				idum=IA1*(idum-k*IQ1)-k*IR1;
				if (idum < 0) idum += IM1;
				if (j < NTAB) iv[j] = idum;
			}
			iy=iv[0];
			seeded = true;
		}
	}  

	double d01()
	{
		k=(idum)/IQ1; 
		//Start here when not initializing.
		idum=IA1*(idum-k*IQ1)-k*IR1; 
		//Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method. 
		if (idum < 0) idum += IM1;
		k=idum2/IQ2;
		idum2=IA2*(idum2-k*IQ2)-k*IR2; 
		//Compute idum2=(IA2*idum) % IM2 likewise.
		if (idum2 < 0) idum2 += IM2;
		j=iy/NDIV; 
		//Will be in the range 0..NTAB-1.
		iy=iv[j]-idum2; 
		//Here idum is shuffled, idum and idum2 are combined to generate output. 
		iv[j] = idum;
		if (iy < 1) iy += IMM1;
		if ((temp=AM*iy) > RNMX) 
		{
			//cout << "random call = " << "RNMAX" << "\n";
			return RNMX; //Because users don't expect endpoint values.
		}
		else 
		{
			return temp;
		}
	}  

};

