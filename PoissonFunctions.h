#define _POLIFITGSL_H
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//semente PRNG. 
#define SEED time(NULL)
//para nunca repetir mesma hist. termica: time(NULL)


/*-----------------------------------------------------------------------------------------------------------------*/
/*----------P. Random Number Generator by Parisi & Rapuano (Lucas Nicolao UFSC)------------------------------------*/
/*--------- Será que ele tem tempo de repetição muito curto, e é melhor usar o gerador do Numerical Recipes?-------*/
/*-----------------------------------------------------------------------------------------------------------------*/
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)
#define PI 3.141592654

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

unsigned rand4init(void) {
  unsigned long long y;
  
  y = (myrand*16807LL);
  myrand = (y&0x7fffffff) + (y>>31);	  	  	                  
  if (myrand&0x80000000)
    myrand = (myrand&0x7fffffff) + 1;
  return myrand;
}

void Init_Random(void) {
  unsigned i;
  
  ip=128;
  ip1=ip-24;
  ip2=ip-55;
  ip3=ip-61;
  
  for (i=ip3; i<ip; i++)
    ira[i] = rand4init();
}
/*---------------------------------------------------------------------*/

/*-------------------------------------------------
obs is the list size, on the format (dx,dy)
degree is the degree of polynomial
*store have the fit coefficients, in crescent order
function got in https://rosettacode.org/wiki/Polynomial_regression#C
--------------------------------------------------*/
bool polynomialfit(int obs, int degree, double *dx, double *dy, double *store) /* n, p */   
{
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;
 
  int i, j;
 
  X = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
 
  for(i=0; i < obs; i++) {
    for(j=0; j < degree; j++) {
      gsl_matrix_set(X, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
  }
 
  ws = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);
 
  /* store result ... */
  for(i=0; i < degree; i++)
  {
    store[i] = gsl_vector_get(c, i);
  }
 
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  return true; /* we do not "analyse" the result (cov matrix mainly)
		  to know if the fit is "good" */
}

//------------------------------------------------ Generate a random variable according to rho(y) from fit

void rejeita(double coeff[],double ymax,double max, double meanx,int n,double xy[]) {
	double y,a,b,zmax_left, zmax_right;	
	double x,function=0;
	double valor=50*max;
	int i;
	unsigned rand4init(void);

	//Parameter for f_y = a |x| + b
	a = 0.009;		
	b=0.8;		
	zmax_right = a/2*max*max+b*max;
//	zmax_left = a/2*min*min+b*fabs(min);

	do {
		function = 0.;		

// Here the x is uniformly sorted, which gives little statistics in the tails. I'll use a linear function a|x|+b to the sort.
		x = (FRANDOM-0.5)*(2*(max))+meanx; printf("x=%lf\n");    /*I'm assuming that interval is symmetric, max = -min*/
		                              /*Estou fazendo uma distribuição uniforme - vale a pena? Acho que tinha que avaliar melhor nas caudas*/
		
		//x = (FRANDOM - (max+m)/2)*2*(max+m);
		//if(FRANDOM > 0.5) {x = -b/a + 1./a*sqrt(b*b+2*a*(FRANDOM*zmax_right));}
		//else {x= (-b/a + 1./a*sqrt(b*b+2*a*(FRANDOM*zmax_left))); x = -x;}
		//printf("x = %lf\n",x);

//--------------------

		y = FRANDOM*ymax;

		for(i=0; i <= n; i++)
		{
			function = function + coeff[i]*pow(x,i);
		}
		
		if (y < pow(10.,function)){
			if(y>pow(10.,function)) {printf("y = %lf, func = %lf\n",y,pow(10.,function));}
			valor = x;
			xy[0]=x;
			xy[1]=y;
			}
		} while (valor == 50*max);

		if(xy[0] < 0) {printf("x negativo YESSS \n");}
}
//--------------------------------------
void rejeita_gsl(gsl_rng *r,double coeff[],double ymax,double max, double meanx,int n,double xy[]) {
	double y;	
	double x,function=0;
	double valor=50*max;
	int i;
	unsigned rand4init(void);
	
//	printf("ymax = %lf\n",ymax);

	do {
		function = 0.;		

//		x = (FRANDOM-0.5)*(2*(max))+meanx; printf("x=%lf\n");    /*I'm assuming that interval is symmetric, max = -min*/
		x = gsl_ran_flat(r,-1.*max,max) + meanx; //printf("x=%lf\n",x);		

//--------------------

		y = gsl_ran_flat(r,0.,1.)*ymax;
//		y = FRANDOM*ymax;

		for(i=0; i <= n; i++)
		{
			function = function + coeff[i]*pow(x,i);
		}

//			if(y>pow(10.,function)) {printf("y = %lf, func = %lf\n",y,pow(10.,function));}
		
		if (y < pow(10.,function)){
			valor = x;
			xy[0]=x;
			xy[1]=y;
			}
		} while (valor == 50*max);
}


//------------------------------------ To prepare the K(x) list (discrete sorting)
double fat(int x)
{
	double f=1.;
	int i;

	for(i=1; i<= x; i++) {f=f*x;}
	return f;
}

int prepare_pdf(int size, double deltaN, double x[], double cdf[])
{

	int i,limit=0;
	double sum = 0;
	double delta, first;


	for(i=0;i<= size; i++)
	{ 
		x[i] = pow(3./2., deltaN*2./9.-i/3.);
//		printf("x[%d] = %lf\n",i,x[i]);
		delta = pow(deltaN*2.*log(3./2.),i)*exp(-deltaN*2.*log(3./2.))/fat(i);
		if(i == 0) {first = delta;}
//		printf("x[%d] = %lf, delta = %e\n",i,x[i],delta);
		sum = sum + delta;
		cdf[i] = sum;
		if (delta < 0.000000001*first) 
		{
			limit = i;			
			break;
		}
	}

	
	for(i=0;i <= limit; i++)
	{
		cdf[i] = cdf[i]/sum;
//		printf("%2.10lf %2.10lf\n",x[i],cdf[i]);
	}

	return limit;  //limit is the length of this K(x) list
}

double sorting(int limite, double x[], double cdf[])  /*generate K(x) distributed in discrete variables*/
{
	unsigned rand4init(void);
	double eps=3.;
	int i,index=0;
	double dif,rando;
	
	rando = FRANDOM;

	for (i=0;i<= limite; i++)
	{
		dif = cdf[i]-rando;
		//printf("dif = %lf\n",dif);
		if(dif < eps && dif > 0) 
		{
			eps = dif;
			index = i;
		}
	}
	
//	printf("indice[min] = %d\n",index);
	return x[index];
}

double poisson_gsl_sorting(gsl_rng *r, float deltaN){
	double mu = 2./3.;
	double c = 2.*log(1.5);
	double a = 1.5;
	int var; 

	var = gsl_ran_poisson(r,c*(double)deltaN);
	return pow(a,deltaN*mu/3. - var/3.);
}

double gauss_sorting(gsl_rng *r,float deltaN){ /*generate G(x) distrinuted in gaussian pdf*/
	double var;
	double smed = (double)deltaN * 0.1;
	var = smed + gsl_ran_gaussian_ziggurat(r,sqrt(2.*smed/log(2)));
	return pow(2,-fabs(var));
//	return var;
}

int fitting(FILE *fin, int degree, int *N, double *xmax, double *ymax, double *coeff, double *meanu){
	FILE *fit = NULL;
	FILE *datatest = NULL;
	double *u = NULL, *v = NULL, *tempcoeff = NULL, aux=0.;
	double temp[2];
	int i;

	//Doing the fit
	i=0;
	while(fscanf(fin,"%lf\t%lf\n",&temp[0],&temp[1]) != EOF) {if(temp[1] > 0) {i++;}}
	fit = fopen("./fit.gp","w");
	datatest = fopen("./datatest.dat","w");

	*N=i;
	u = malloc((*N)*sizeof(double));
	v = malloc((*N)*sizeof(double));
	tempcoeff = malloc(degree*sizeof(double));

	i=0;
	rewind(fin);
	*meanu = 0.;
	while(fscanf(fin,"%lf\t%lf\n",&temp[0],&temp[1]) != EOF){
//		if(i==3) {printf("temp[0] = %lf\n",temp[0]);}		
		if(temp[1] != 0){
			u[i] = temp[0];
			if(v[i] > aux) {aux=v[i]; *meanu=u[i];}
			//*meanu += temp[0];
			//v[i] = temp[1];			
			v[i] = log10(temp[1]);
			if (fabs(u[i]) > (*xmax)) {*xmax = fabs(u[i]);}
			if (pow(10,v[i]) > (*ymax)) {*ymax = pow(10,v[i]);}			
			fprintf(datatest,"%lf\t%lf\n",u[i],v[i]);
			i++;
			}
	}

	//*meanu = *meanu/((double)(*N));

	if(polynomialfit(*N,degree,u,v,tempcoeff) != true) {printf("Error in fit \n"); return 1;}
	for(i=0;i<degree;i++) {coeff[i]=tempcoeff[i];}

	fprintf(fit,"f(x) = ");	
	for(i=0;i<degree;i++) {fprintf(fit,"+(%2.20lf*x**%d)",coeff[i],i);}
	fprintf(fit,"\n\nset term x11 0\nset xr [-%lf:%lf]\nset yr [-2*10**(11):0]\nplot f(x), './datatest.dat'",(*xmax)*1.2,(*xmax)*1.2);
	fprintf(fit,"\n\nset term x11 1\nset xr [-%lf:%lf]\nset yr [-2*10**(1):0]\nplot f(x), './datatest.dat'",(*xmax)*1.1,(*xmax)*1.1);
	//fprintf(fit,"f[x_] := ");	
	//for(i=0;i<degree;i++) {fprintf(fit,"+(%2.15lf*x^%d)",coeff[i],i);}	

	fclose(fit);
	fclose(datatest);
	//system("gedit './fit.dat'");
	//system("gnuplot load 'fit.gp'");
	system("gnuplot -p 'fit.gp'");
}

//---------------------------------------------------
void histogram(FILE *fin, FILE *fout, int Nbins){
	int max, aux,N;
	int i, j;
	double temp, mean,sum,sigma=0.;
	double lastbin=0., firstbin=0., lbin;

	float *count = NULL;
	count = malloc(Nbins*sizeof(float));

	for(i=0;i<Nbins;i++) {count[i]=0.;}

	N = 0;
	mean=0.;
	while(fscanf(fin,"%lf\n",&temp) != EOF){
		if(temp < firstbin) {firstbin = temp;}
		if(temp > lastbin) {lastbin = temp;}
		mean+=temp;		
		N++;
		}

	mean = mean/((double) N);

	rewind(fin);
	sigma=0.;
	while(fscanf(fin,"%lf\n",&temp) != EOF){
		sigma+=pow(temp-mean,2);
		}

	sigma=sigma/((double)N);
	sigma=sqrt(sigma);

	printf("mean=%lf\nsigma=%lf\n",mean,sigma);

	//sigma=1.;  //only for test, ERASE IT!!

	firstbin = firstbin/sigma;
	lastbin = lastbin/sigma;
	lbin=(lastbin-firstbin)/Nbins;


	rewind(fin);
	while(fscanf(fin,"%lf\n",&temp) != EOF){
		//temp = (temp-mean)/sigma;
		temp = temp/sigma;		
		aux = (int)((temp-firstbin)/lbin);
		//printf("aux = %d\n",aux);
		if (aux<0) {printf("aux negativo %d\n",aux);}
		if(aux>Nbins) {printf("Alerta de aux > Nbins!!\n");}
		if(aux==Nbins) {count[(int)(Nbins-1)]++;}
		else {count[aux]++;} 		
		}

	for(j=0;j<Nbins;j++){
		count[j]=count[j]/((double)N*lbin);
//		//sum+=count[j]*lbin;
		}		

	for(i=0;i<Nbins;i++) {
//		if(count[i] != 0) {
			fprintf(fout,"%lf\t%2.9lf\n",firstbin+i*lbin,count[i]);
//			}
		}
}

double hflatness_by_pdf(FILE *fin, int q){
	int N=0;
	double min,max,lbin,flat;	
	double temp[2],frac[2];
	frac[0] = 0.; frac[1]=0.;
	
	while(fscanf(fin,"%lf\t%lf\n",&temp[0],&temp[1]) != EOF){
		frac[0] += pow(fabs(temp[0]),q)*temp[1];
		frac[1] += pow(fabs(temp[0]),3)*temp[1];
		if(N==0) {min=temp[0];}
		N++;		
		}

	max = temp[0];
	lbin = (max-min)/((double)N);
	frac[0] = frac[0]*lbin;
	frac[1] = frac[1]*lbin;

	//printf("lbin = %lf,frac[0] = %lf, frac[1] = %lf\n",lbin,frac[0],frac[1]);

	flat = frac[0]/pow(frac[1],(double)q/3.);
	printf("hflatness = %lf, (double)q/3. = %lf\n",flat,(double)q/3.);
	rewind(fin);
	return flat;
}

//---------------------------------------

double cascade_steps_poisson(FILE *f1, FILE *f2, int q){
	double flat1, flat2, deltaN;
	flat1=hflatness_by_pdf(f1,q);
	flat2=hflatness_by_pdf(f2,q);
	printf("flat1 = %lf, flat2 = %lf\n",flat1,flat2);	

	deltaN = log(flat1/flat2)/(log(3./2.)*(2./9.*(double)q-2.*(1.-pow(3./2.,-(double)q/3.))));
	return deltaN;
}

//---------------------------------------

double cascade_steps_gauss(FILE *f1, FILE *f2, int q){
	double flat1, flat2, deltaN;
	flat1=hflatness_by_pdf(f1,q);
	flat2=hflatness_by_pdf(f2,q);
	printf("flat1 = %lf, flat2 = %lf\n",flat1,flat2);	

	deltaN = log(flat1/flat2)/(log(2.)*0.1*(q*q/9. - q/3.));
	rewind(f1);
	rewind(f2);
	return deltaN;
}

