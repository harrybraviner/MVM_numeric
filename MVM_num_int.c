#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* The purpose of this code is to numerically integrate equations (2.18), (2.19) and (2.20) from hep-th/1011.3502v1 with initial data representing a small displacement from AdS_d.
This is very similar to the numerical work found near the end of hep-th/0808.1725v2, but we work with the massive vector model and general z */

// Physics constants
const double d = 4;		// Number of spacetime dimensions
const double l_UV = 1;		// The length-scale of the target spacetime in the UV (l in the e.o.m.s)
double z = 2;			// The dynamical exponent of the target spacetime in the UV
double a_0 = 0.01;		// The 'size' of the perturbation
const double K = 1;		// Coefficient in front of r^2*dt^2 in AdS_d. This is purely a gauge choice.

// Programming constants
/* Note: work with rho = log r, since this allows for a much greater range */
double rho_init = 1;			// Start-point for the numerical integration
const unsigned long N_rho = 50000;	// The number of rho steps
double delta_rho = 0.01;			// Step-length
char *outfileName = "output.dat";	// File to which the output will be written



static double m_0, Lambda, l_IR; // Declare these to be global now - will calculate them inside of main()

static unsigned long i;		// Dummy counter variable
// Function declarations for derivatives
double A_prime(double rho, double A, double B, double C, double C_prime);
double B_prime(double rho, double A, double B, double C, double C_prime);
double C_prime_prime(double rho, double A, double B, double C, double C_prime);

int main(int argc, char **argv){
	
	if(argc>1){
		// Set the value of z to the first input
		z = strtod(argv[1],NULL);
	}
	if(argc>2){
		// Set the value of a_0 to the second input
		a_0 = strtod(argv[2],NULL);
	}
	if(argc>3){
		rho_init = strtod(argv[3],NULL);
	}
	if(argc>4){
		delta_rho = strtod(argv[4],NULL);
	}
	printf("Setting m_0^2, Lambda to values for z = %f\n",z);
	printf("\"Size\" of the perturbation set to a_0 = %f\n",a_0);
	printf("rho_init = %f\n",rho_init);
	printf("delta_rho = %f\n",delta_rho);

	FILE *outfile = fopen(outfileName, "w");
	if(outfile == NULL){
		fprintf(stderr, "Unable to open %s in write mode. Quitting.\n",outfileName);
		exit(0);
	}

	/* Calculate various physical quantities based on the parameters supplied */
	m_0 = sqrt((d-2)*z)/l_UV;	// Set the mass of the vector field to the necessary value to attain the UV target
	Lambda = -(z*z + (d-3)*z + (d-2)*(d-2))/(2*l_UV*l_UV);	// Cosmological constant
	l_IR = l_UV*sqrt((d-1)*(d-2)/(z*z+(d-3)*z+(d-2)*(d-2)));	// This is used to set the initial data

	/* Note that we do not work directly with the fields used by CM, as these will become very large at large r. Instead we use A = ln(|f_1|)/rho, B = ln(W)/rho and C = ln(alpha)/rho. All derivatives (including "prime" in this code) are with respect to rho */
	double A, B, C, C_prime;
	double kappa_est,alpha_est,c,d;	// variables for curve-fitting
	long est_step = 60;
#if 0
	/* Initialise A,B to the appropriate values for AdS_d */
	printf("Setting the spacetime at rho = %f to AdS\n",rho_init);
	A = log(K)/rho_init + 2;	// Could possibly generalise this to another Lifshitz space in the IR by replacing the 2 here
	B = 2*log(l_UV/l_IR)/rho_init + 2;
	/* set the gauge field to obey the fall-off condition found earlier, alpha ~ a_0*r^beta */
	a_0 = 0.01;
	double beta = -(d-3)/2.0 + sqrt((d-3)*(d-3)/4.0 + (d-1)*z/(z*z+(d-3)*z+(d-2)*(d-2)));
	C = log(a_0)/rho_init + beta;
	C_prime = -log(a_0)/(rho_init*rho_init);
#endif
#if 1
	/* Set up the spacetime near rho=rho_init to look like a Lifshitz spacetime of length-scale l_UV */
	printf("Setting the spacetime at rho = %f to Li_%f\n",rho_init,z);
	A = log(K)/rho_init + 2*z;
	B = 2;
	C = (0.5*log(2*(z-1)/z)+log(l_UV))/rho_init + z;
	C_prime = -(0.5*log(2*(z-1)/z)+log(l_UV))/(rho_init*rho_init);
#endif

	kappa_est=0;
	alpha_est=0;
	c=0,d=0;

	/* Write column headers to the output file, then write the inital values */
	fprintf(outfile,"# rho\tln|f|/rho\tln(W)/rho\tln(alpha)/rho\td(ln(alpha)/rho)/drho\tkappa estimate\td\talpha estimate\tc\n");
	fprintf(outfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",rho_init,A,B,C,C_prime,kappa_est,d,alpha_est,c);

	// variables used internally by the RK4 loop 
	double Ak_1,Ak_2,Ak_3,Ak_4;
	double Bk_1,Bk_2,Bk_3,Bk_4;
	double Ck_1,Ck_2,Ck_3,Ck_4;
	double C_primek_1,C_primek_2,C_primek_3,C_primek_4;
	double B1,B2;
	for(i=0;i<N_rho;i++){
		/* Use RK4 method to numerically integrate */
		/* k_1 */
		Ak_1 = A_prime(rho_init + i*delta_rho, A, B, C, C_prime);
		Bk_1 = B_prime(rho_init + i*delta_rho, A, B, C, C_prime);
		Ck_1 = C_prime;
		C_primek_1 = C_prime_prime(rho_init + i*delta_rho, A, B, C, C_prime);
		/* k_2 */
		Ak_2 = A_prime(rho_init + (i+0.5)*delta_rho, A + 0.5*delta_rho*Ak_1, B + 0.5*delta_rho*Bk_1, C + 0.5*delta_rho*Ck_1, C_prime + 0.5*delta_rho*C_primek_1);
		Bk_2 = B_prime(rho_init + (i+0.5)*delta_rho, A + 0.5*delta_rho*Ak_1, B + 0.5*delta_rho*Bk_1, C + 0.5*delta_rho*Ck_1, C_prime + 0.5*delta_rho*C_primek_1);
		Ck_2 = C_prime + 0.5*delta_rho*C_primek_1;
		C_primek_2 = C_prime_prime(rho_init + (i+0.5)*delta_rho, A + 0.5*delta_rho*Ak_1, B + 0.5*delta_rho*Bk_1, C + 0.5*delta_rho*Ck_1, C_prime + 0.5*delta_rho*C_primek_1);
		/* k_3 */
		Ak_3 = A_prime(rho_init + (i+0.5)*delta_rho, A + 0.5*delta_rho*Ak_2, B + 0.5*delta_rho*Bk_2, C + 0.5*delta_rho*Ck_2, C_prime + 0.5*delta_rho*C_primek_2);
		Bk_3 = B_prime(rho_init + (i+0.5)*delta_rho, A + 0.5*delta_rho*Ak_2, B + 0.5*delta_rho*Bk_2, C + 0.5*delta_rho*Ck_2, C_prime + 0.5*delta_rho*C_primek_2);
		Ck_3 = C_prime + 0.5*delta_rho*C_primek_2;
		C_primek_3 = C_prime_prime(rho_init + (i+0.5)*delta_rho, A + 0.5*delta_rho*Ak_2, B + 0.5*delta_rho*Bk_2, C + 0.5*delta_rho*Ck_2, C_prime + 0.5*delta_rho*C_primek_2);
		/* k_4 */
		Ak_4 = A_prime(rho_init + (i+1)*delta_rho, A + delta_rho*Ak_3, B + delta_rho*Bk_3, C + delta_rho*Ck_3, C_prime + delta_rho*C_primek_3);
		Bk_4 = B_prime(rho_init + (i+1)*delta_rho, A + delta_rho*Ak_3, B + delta_rho*Bk_3, C + delta_rho*Ck_3, C_prime + delta_rho*C_primek_3);
		Ck_4 = C_prime + delta_rho*C_primek_3;
		C_primek_4 = C_prime_prime(rho_init + (i+1)*delta_rho, A + delta_rho*Ak_3, B + delta_rho*Bk_3, C + delta_rho*Ck_3, C_prime + delta_rho*C_primek_3);
		/* Calculate the next values */
		A = A + (delta_rho/6)*(Ak_1 + 2*Ak_2 + 2*Ak_3 + Ak_4);
		B = B + (delta_rho/6)*(Bk_1 + 2*Bk_2 + 2*Bk_3 + Bk_4);
		C = C + (delta_rho/6)*(Ck_1 + 2*Ck_2 + 2*Ck_3 + Ck_4);
		C_prime = C_prime + (delta_rho/6)*(C_primek_1 + 2*C_primek_2 + 2*C_primek_3 + C_primek_4);
		/* If we have enough data to do so, estimate the value of z and kappa */
		if(i%est_step==0){
			B2=B1;
			B1 = B;
			if(i>=2*est_step){
				kappa_est = log((rho_init+(i+1)*delta_rho)*(B1-2)/((rho_init+(i-est_step+1)*delta_rho)*(B2-2)))/(delta_rho*est_step);
				alpha_est = log((B1-2)/(B2-2))/(est_step*delta_rho);
				d = (rho_init+delta_rho*(i+1))*(B-2)/exp(kappa_est*(rho_init+delta_rho*(i+1)));
				c = (B-2)/exp(kappa_est*(rho_init+delta_rho*(i+1)));
			}
		}
		fprintf(outfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",rho_init+delta_rho*(i+1),A,B,C,C_prime,kappa_est,d,alpha_est,c);
	}

	fclose(outfile);

	return 0;
}

// Function definitions for derivatives
double A_prime(double rho, double A, double B, double C, double C_prime){
	double deriv = - A/rho - (d-3)/rho - (2*l_UV*l_UV*Lambda/(d-2))*exp(rho*(2-B))/rho
	+ (m_0*m_0/(2*(d-2)))*exp(rho*(2-B+2*C-A))/rho
	- (1/(2*(d-2)*l_UV*l_UV))*(rho*C_prime*C_prime + 2*C*C_prime + C*C/rho)*exp(rho*(2*C-A));
	return deriv;
}

double B_prime(double rho, double A, double B, double C, double C_prime){
	double deriv = - B/rho - (d-3)/rho - (2*l_UV*l_UV*Lambda/(d-2))*exp(rho*(2-B))/rho
	- (m_0*m_0/(2*(d-2)))*exp(rho*(2-B+2*C-A))/rho
	- (1/(2*(d-2)*l_UV*l_UV))*(rho*C_prime*C_prime + 2*C*C_prime + C*C/rho)*exp(rho*(2*C-A));
	return deriv;
}

double C_prime_prime(double rho, double A, double B, double C, double C_prime){
	double deriv = (3-d - 2/rho + (m_0*m_0/(2*(d-2)))*exp(rho*(2-B+2*C-A)))*(C_prime + C/rho)
	- (rho*C_prime*C_prime + 2*C*C_prime + C*C/rho)
	+ 2*C/(rho*rho) + (l_UV*l_UV*m_0*m_0)*exp(rho*(2-B))/rho;
	return deriv;
}
