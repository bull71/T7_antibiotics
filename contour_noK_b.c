// file contour_noK_b.c to run recurrence equations of ODEs and output to a file for upload and analysis; accompanies Wilcox et al.:  Phage T7 adaptation to two antibiotics
// Here burst size (b) is varied to change phage growth rate (Fig. 5B)
// 
#include <stdio.h>
#include <math.h>
FILE *fp; //*fopen();

int main()
{

char out[30] ="contour_noK_b.csv" ;
int counter,counter1, I,J,i,j,N,innerloop,C,printcycles,K,klim;

double 
z,lk,Y,Z,a,delta,die,h, r,b,k,L,B,Bx,Pi,Pix,Pf,Pfx,Pf0,B0,tiny,D;

//D = 1.0E8;  // is the cell density used to measure phage growth rate
//printf("input h, N\n");
//scanf("%lf %d",&h,&N);

fp = fopen(out,"a+");

h = .01;
N = 2000;
klim = (int)(1.0/h);
D = 1.0E8;

printf("input log(k),  lysis rate, decay rate\n");
scanf("%lf %lf %lf",&lk,&L,&delta);
printf("data input\n");
 
fprintf(fp,"log(k) = %1.3f, lysis rate = %1.4f decay rate = %1.3f\n",lk,L,delta);
fprintf(fp,"growth rate, cells, Pf0\n");

k = pow(10.0,lk);

for(B0 = 4.0;B0<8.5; B0+= 0.2)				{
for(b = 1.0;b<100;b*=1.2)		{
r = 0.000241;
tiny = 0.0;
die = .001;



for(Pf0 = 1;Pf0<1E12;Pf0 *=2.0)				{ // Pf0 loop

B=Bx= 0.0;
Pi=Pix=Pf=Pfx=0.0;
Pf = Pf0;
B =  pow(10.0,B0);;

for(I=1;I<=N;I++)					{

for(K=1;K<=klim;K++)				{   // K loop reciprocal of h

Bx = B + h*( r*B   - k*Pf*B - die*B +tiny);

Pix = Pi + h*(k*Pf*B - L*Pi - die*Pi + tiny);
Pfx = Pf + h*(b*L*Pi - k*Pf*B - delta*Pf) ;

B = Bx;

Pi = Pix;
Pf = Pfx;

} // end of K loop (Euler)

Z = 0.5*(-die - D*k - delta - L + 
sqrt((die + D*k + delta + L)*(die +D*k + delta +L) -4.0*(-b*L*D*k +D*die*k +die*delta + D*k*L +L*delta) ));
// Z is growth rate for these equations assuming a fixed cell density (D)

z = Z*60.0/(log(2.0));  // convert growth rate to doublings/hr
if(B < .01*pow(10.0,B0)) 	{
		fprintf(fp,"%1.3f, %1.4f, %1.4f\n",z,B0,log10(Pf0));
		Pf0 = 1E13;
		I=N+1;
		}
} // end of N loop

}  //P0
}  //delta loop
}  // k loop

fclose(fp);



}


