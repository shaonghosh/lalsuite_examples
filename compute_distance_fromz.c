#include <stdio.h>
#include <math.h>

int main() {
	long double I, z=0, S=0, redshift=1.0, dz=1e-9;
	long double c=299792458.0, H0=67.9, OmegaM=0.3065, OmegaL=0.6935;
	while(z <= redshift) {
		I = 1./pow((OmegaM*pow((1 + z), 3) + OmegaL), 0.5);
		S += I*dz;
		z += dz;
	}
	printf("Value of luminosity distance = %Lf", (c/1000)*(1+z)*S/H0);
}




