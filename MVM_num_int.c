#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Bfitpower
//#define Bfitlog
#define Ffitpower

//#define AdS_IR
#define Li_IR

/* The purpose of this code is to numerically integrate equations (2.18), (2.19) and (2.20) from hep-th/1011.3502v1 with initial data representing a small displacement from AdS_d or a Lifshitz spacetime
This is very similar to the numerical work found near the end of hep-th/0808.1725v2, but we work with the massive vector model and general z */

// Physics constants
const double d = 6;		// Number of spacetime dimensions
const double L = 1.0;		// The length-scale appearing in the ansatz, and hence the equations of motion.
double z = 2;			// The dynamical exponent of the Lishitz spacetime to which \Lambda, m_0 are appropriate
double a_0 = -0.00000001;		// The 'size' of the perturbation
const double K = 1;		// Coefficient in front of r^2*dt^2 in AdS_d. This is purely a gauge choice.

// Programming constants
/* Note: work with rho = log r, since this allows for a much greater range */
double rho_init = 1;			// Start-point for the numerical integration
const unsigned long N_rho = 50000;	// The number of rho steps
double delta_rho = 0.001;			// Step-length
char *outfileName = "output.dat";	// File to which the output will be written
long est_step = 100;	// number of steps to take between points for curve-fitting estimates



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

	/* Calculate various physical quantities based on the parameters supplied
	   Specifically, we want to choose \Lambda, m_0 such that they both support
	   a Lifshitz spacetime of dynamical exponent z, and whatever spacetime is
	   in the IR with length-scale 1 */
#ifdef AdS_IR
	Lambda = -0.5*(d*d - 3*d + 2);
	m_0 = sqrt(-(2*Lambda*z*(d-2))/(z*z + (d-3)*z + (d-2)*(d-2)));
#endif
#ifdef Li_IR
	Lambda = -0.5*(z*z + (d-3)*z + (d-2)*(d-2));
	m_0 = sqrt(z*(d-2));
#endif

	/* Note that we do not work directly with the fields used by CM, as these will become very large at large r. Instead we use A = ln(|f_1|)/rho, B = ln(W)/rho and C = ln(alpha)/rho. All derivatives (including "prime" in this code) are with respect to rho */
	double A, B, C, C_prime;
#ifdef AdS_IR
	/* Set up the spacetime at rho_init as an AdS space of length-scale 1 */
	printf("Setting the spacetime at rho = %f to AdS\n",rho_init);
	A = (log(K)-2*log(L))/rho_init + 2;
	B = 2*log(L)/rho_init + 2;
	/* Add a perturbation of size a_0 along the irrelevant direction */
	C = log(a_0)/rho_init + (0.5*(1-d) + 0.5*sqrt(4*m_0*m_0 + (d-3)*(d-3)));
	C_prime = -log(a_0)/(rho_init*rho_init);
#endif

#ifdef Li_IR
	/* Set up the spacetime near rho=rho_init to look like a Lifshitz spacetime of length-scale L */
	printf("Setting the spacetime at rho = %f to Li_%f\n",rho_init,z);
	A = log(K)/rho_init + 2*z;
	B = 2*log(L)/rho_init + 2;
	C = z + log(2*(d-2)*(z-1))/(2*rho_init) - log(m_0)/rho_init;
	C_prime = -(0.5*log(2*(d-2)*(z-1)) - log(m_0))/(rho_init*rho_init);
	/* Add a perturbation of size a_0 along the irrelevant direction */
	A += +a_0*sqrt(2*z)*(z+d-3)*(z+d-2)*(z+d-2+sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2)))*exp(0.5*(2-d-z+sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2)))*log(rho_init))/(K*rho_init);
	B += a_0*sqrt(2*z)*(z-1)*(z+d-2)*(z-3*d+6*sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2)))*exp(0.5*(2-d-z+sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2)))*log(rho_init))/rho_init;
	C += a_0*4*L*z*sqrt(z-1)*(z+d-2)*(z+d-3)*exp(0.5*(2-d-z+sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2)))*log(rho_init))/(L*sqrt(2*(z-1)/z)*rho_init);
	C_prime += (0.5*(2-d-z+sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2))) - 1)*(a_0*4*L*z*sqrt(z-1)*(z+d-2)*(z+d-3)*exp(0.5*(2-d-z+sqrt(9*z*z+(4-6*d)*z+(d+6)*(d-2)))*log(rho_init))/(L*sqrt(2*(z-1)/z)*rho_init))/rho_init;
#endif

	/* Variable definitions for curve fitting */
#ifdef Bfitpower
	// Variable definitions for attempting to fit the curve to W ~ a r^2 + ab r^{\kappa+2}
	double kappa,a,b;
	double B1=0,B2=0,B3=0;
	double rho1=0,rho2=0,rho3=0;
	// Open a suitable output file, and write headers
	FILE *Bfitpower_file = fopen("Bfitpower.dat","w");
	fprintf(Bfitpower_file,"rho\tkappa\ta\tb\n");
	printf("Attempting a power-law fit to the subleading contributions to W. This is being written to Bfitpower.dat\n");
#endif

#ifdef Bfitlog
	// Variable definitions for attempting to fit the curve to W ~ a r^2 + ab r^{\kappa+2}log(r)
	double kappa_B_log,a_B_log,b_B_log;
	double pB1_log=0,pB2_log=0,pB3_log=0;
	// Open a suitable output file, and write headers
	FILE *Bfitlog_file = fopen("Bfitlog.dat","w");
	fprintf(Bfitlog_file,"rho\tkappa\ta\tb\n");
	printf("Attempting a power-law*log fit to the subleading contributions to W. This is being written to Bfitlog.dat\n");
#endif

#ifdef Ffitpower
	// Variable definitions for attempting to fi the curve F ~ a_F r^{2 z_F} + a_F b_F r^{2 z_F + \kappa_F}
	double z_F,kappa_F,a_F,b_F;
	double pA1, pA2, pA3, pA4;
	// Open a file to write output to, and write headers
	FILE *Ffitpower_file = fopen("Ffitpower.dat","w");
	fprintf(Ffitpower_file,"rho\tkappa_F\ta_F\tb_F\tz_F\n");
	printf("Attempting to estimate a power-law scaling for F, fit to the subleading contributions to a power-law. This is being written to Ffitpower.dat\n");
#endif

	/* Write column headers to the output file, then write the inital values */
	fprintf(outfile,"# rho\tln|f|/rho\tln(W)/rho\tln(alpha)/rho\td(ln(alpha)/rho)/drho\n");
	fprintf(outfile,"%f\t%f\t%f\t%f\t%f\n",rho_init,A,B,C,C_prime);

	// variables used internally by the RK4 loop 
	double Ak_1,Ak_2,Ak_3,Ak_4;
	double Bk_1,Bk_2,Bk_3,Bk_4;
	double Ck_1,Ck_2,Ck_3,Ck_4;
	double C_primek_1,C_primek_2,C_primek_3,C_primek_4;
	double rho;
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
		rho = rho_init+delta_rho*(i+1);
		// Write solution points to a file
		fprintf(outfile,"%f\t%f\t%f\t%f\t%f\n",rho,A,B,C,C_prime);

		/* Try to fit to an appropriate curve */
#ifdef Bfitpower
		if(i%est_step==0){
			B3=B2,B2=B1,B1=B;
			rho3=rho2,rho2=rho1,rho1=(rho_init+delta_rho*(i+1));
			kappa = log((rho1*B1-rho2*B2-2*est_step*delta_rho)/(rho2*B2-rho3*B3-2*est_step*delta_rho))/(est_step*delta_rho);
			b = (rho1*(B1-2) - rho2*(B2-2))/(exp(kappa*rho1)-exp(kappa*rho2));
			a = exp(rho1*(B1-2)-b*exp(kappa*rho1));
			fprintf(Bfitpower_file,"%f\t%f\t%f\t%f\n",rho1,kappa,a,b);
		}
#endif

#ifdef Bfitlog
		// FIXME: Do not yet have a working method to extract these numbers
#endif

#ifdef Ffitpower
		if(i%est_step==0){
			pA4=pA3,pA3=pA2,pA2=pA1,pA1=(rho_init+delta_rho*(i+1))*A; 
			z_F = (pA2*pA2 + pA3*pA3 - pA1*pA3 - pA2*pA3 + pA1*pA4 - pA2*pA4)/(2*est_step*delta_rho*(-pA1 + 3*pA2 - 3*pA3 + pA4));
			kappa_F = log((pA1 - pA2 - 2*z_F*est_step*delta_rho)/(pA2 - pA3 - 2*z_F*est_step*delta_rho))/(est_step*delta_rho);
			b_F = (pA1 - pA2 - 2*z_F*est_step*delta_rho)/(exp(kappa_F*(rho_init+delta_rho*(i+1)))*(1-exp(-kappa_F*est_step*delta_rho)));
			a_F = exp(pA1-(rho_init+delta_rho*(i+1))*2*z_F - b_F*exp(kappa_F*(rho_init+delta_rho*(i+1))));
			fprintf(Ffitpower_file,"%f\t%f\t%f\t%f\t%f\n",(rho_init+delta_rho*(i+1)),kappa_F,a_F,b_F,z_F);
		}
#endif
	}

	fclose(outfile);
#ifdef Bfitpower
	fclose(Bfitpower_file);
#endif
#ifdef Ffitpower
	fclose(Ffitpower_file);
#endif

	return 0;
}

// Function definitions for derivatives
double A_prime(double rho, double A, double B, double C, double C_prime){
	double deriv = - A/rho - (d-3)/rho - (2*L*L*Lambda/(d-2))*exp(rho*(2-B))/rho
	+ (m_0*m_0/(2*(d-2)))*exp(rho*(2-B+2*C-A))/rho
	- (1/(2*(d-2)*L*L))*(rho*C_prime*C_prime + 2*C*C_prime + C*C/rho)*exp(rho*(2*C-A));
	return deriv;
}

double B_prime(double rho, double A, double B, double C, double C_prime){
	double deriv = - B/rho - (d-3)/rho - (2*L*L*Lambda/(d-2))*exp(rho*(2-B))/rho
	- (m_0*m_0/(2*(d-2)))*exp(rho*(2-B+2*C-A))/rho
	- (1/(2*(d-2)*L*L))*(rho*C_prime*C_prime + 2*C*C_prime + C*C/rho)*exp(rho*(2*C-A));
	return deriv;
}

double C_prime_prime(double rho, double A, double B, double C, double C_prime){
	double deriv = (3-d - 2/rho + (m_0*m_0/(2*(d-2)))*exp(rho*(2-B+2*C-A)))*(C_prime + C/rho)
	- (rho*C_prime*C_prime + 2*C*C_prime + C*C/rho)
	+ 2*C/(rho*rho) + (L*L*m_0*m_0)*exp(rho*(2-B))/rho;
	return deriv;
}
