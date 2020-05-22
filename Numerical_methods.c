#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float ff(float, float);
float gg(float);
float f (float); 
float g(float, float);

int main(int argc, const char * argv[]) 
{
    //First declaring all variables needed for the Runge-Kutta, Adams moulton and Gaus-Seidel Methods
	float k0, k1, k2, k3, h=0.05;
	float l0=0, l1=0, l2=0, l3=0, n0=0, n1=0, n2=0 ,n3=0, p0=0, p1=0 ,p2=0 ,p3=0;
	float thA=0.3, thB=0, thX=0, thY=0;
	float *a, *b, *x, *y;
	float p, c, q, r, s, v, w, z, t;
	int i, j, n, nsteps=4, msteps=201, npoints;
	float uu[11][11], pp[11][11];
	float cc, dd, uuu, error1, error2=0, maxerror =1;
	float pi=acos(-1);
	npoints = nsteps + msteps;
	
	a = (float *) malloc (npoints * sizeof(float));
	b = (float *) malloc (npoints * sizeof(float));
	x = (float *) malloc (npoints * sizeof(float));
	y = (float *) malloc (npoints * sizeof(float));
	
	//Runge Kutta 4th Order Method to obtain initial values to be used in Adams Moulton Method
	for (i=0; i<nsteps; i++){
	
	t=i*h;
	
	k0 = h*f(thB);
	k1 = h*f(thB + 0.5*l0);
	k2 = h*f(thB + 0.5*l1);
	k3 = h*f(thB + l2);
	
	l0 = h*g(thA, thX);
	l1 = h*g(thA + 0.5*k0, thX + 0.5*n0);
	l2 = h*g(thA + 0.5*k1, thX + 0.5*n1);
	l3 = h*g(thA + k2, thX + n2);
	
	n0 = h*f(thY);
    n1 = h*f(thY + 0.5*p0);
	n2 = h*f(thY + 0.5*p1);
	n3 = h*f(thY + p2);
	
	p0 = h*g(thX, thA);
	p1 = h*g(thX + 0.5*n0, thA + 0.5*k0);
	p2 = h*g(thX + 0.5*n1, thA + 0.5*k1);
	p3 = h*g(thX + n2, thA + k2);
	
	a[i] = thA; b[i] = thB; x[i] = thX; y[i] = thY;
	printf("%f    %f %f %f %f\n", t, a[i], b[i], x[i], y[i]);
	
	thA = thA + (k0+2*k1+2*k2+k3)/6.0;
	thB = thB + (l0 +2*l1+2*l2+l3)/6.0;
	thX = thX + (n0+2*n1+2*n2+n3)/6.0;
	thY = thY + (p0+2*p1+2*p2+p3)/6.0;
	}
		printf("\n");
	//Beginning of Adams Moulton Method.
	for (i=4; i<msteps; i++){
	
	t=i*h;
	//Initial Predicting formula	
	a[i] = a[i-1] + (h*(55*f(b[i-1]) - 59*f(b[i-2]) + 37*f(b[i-3]) - 9*f(b[i-4])))/24.0;
	b[i] = b[i-1] + (h*(55*g(a[i-1],x[i-1]) - 59*g(a[i-2],x[i-2]) + 37*g(a[i-3],x[i-3]) - 9*g(a[i-4],x[i-4])))/24.0;
	x[i] = x[i-1] + (h*(55*f(y[i-1]) - 59*f(y[i-2]) + 37*f(y[i-3]) - 9*f(y[i-4])))/24.0;
	y[i] = y[i-1] + (h*(55*g(x[i-1],a[i-1]) - 59*g(x[i-2],a[i-2]) + 37*g(x[i-3],a[i-3]) - 9*g(x[i-4],a[i-4])))/24.0;
	p=a[i];
	q=b[i];
	s=x[i];
	w=y[i];
	//First Correcting formula
	a[i] = a[i-1] + (h*(9*f(b[i]) + 19*f(b[i-1]) - 5*f(b[i-2]) + f(b[i-3])))/24.0;
	b[i] = b[i-1] + (h*(9*g(a[i],x[i]) + 19*g(a[i-1],x[i-1]) - 5*g(a[i-2],x[i-2]) + g(a[i-3],x[i-3])))/24.0;
	x[i] = x[i-1] + (h*(9*f(y[i]) + 19*f(y[i-1]) - 5*f(y[i-2]) + f(y[i-3])))/24.0;
	y[i] = y[i-1] + (h*(9*g(x[i],a[i]) + 19*g(x[i-1],a[i-1]) - 5*g(x[i-2],a[i-2]) + g(x[i-3],a[i-3])))/24.0;
	c=a[i];
	r=b[i];
	v=x[i];
	z=y[i];
	
	//While loop to keep correcting the initial correction if the accuracy of the solution is 
	//not good enough.
	while(abs(p-c)>0.001||abs(q-r)>0.001||abs(s-v)>0.001||abs(w-z)>0.001){
	p=c;
	q=r;
	s=v;
	w=z;
    a[i] = a[i-1] + (h*(9*f(b[i]) + 19*f(b[i-1]) - 5*f(b[i-2]) + f(b[i-3])))/24.0;
	b[i] = b[i-1] + (h*(9*g(a[i],x[i]) + 19*g(a[i-1],x[i-1]) - 5*g(a[i-2],x[i-2]) + g(a[i-3],x[i-3])))/24.0;
	x[i] = x[i-1] + (h*(9*f(y[i]) + 19*f(y[i-1]) - 5*f(y[i-2]) + f(y[i-3])))/24.0;
	y[i] = y[i-1] + (h*(9*g(x[i],a[i]) + 19*g(x[i-1],a[i-1]) - 5*g(x[i-2],a[i-2]) + g(x[i-3],a[i-3])))/24.0;
	c=a[i];	
	r=b[i];
	v=x[i];
	z=y[i];
		
	}
	printf("%f   %f %f %f %f\n", t, a[i], b[i], x[i], y[i]);

	}
	
	//Beginning of the Gauss-Seidel Method
	// This sets all the values in the newly created matrix to equal 0, these will be used as 
	// initial guesses for our surface values. 
    for(i=0; i<11; i++)
       {
          for(j=0;j<11;j++)
             {
             	uu[i][j] = 0.00000000;
             	pp[i][j] = 0.00000000;
             }
       }
    // Iterations for the Gauss-Seidel method begin. 
    //While loop to obtain required accuracy from the numerical method
    while (maxerror > 0.00001)
    {
    	for(j=0;j<11;j++)
    	{
    		for(i=0;i<11;i++)
    		{
    			pp[i][j]=uu[i][j];
    		}
    	}
    	for(j=0;j<11;j++)
    	{
    		for(i=j;i<11;i++)
    		{
    			if(j==0)
    			{
    				uu[i][j]=0;
    			}
    			else if(i==j && i==10)
    			{
    				uu[i][j]=(2*pp[i][j-1]+2*uu[i][j-1]-gg(h)*ff(i*h,j*h))*0.25;
    			}
    			else if(i==j && i>0)
    			{
    				uu[i][j]=(2*pp[i+1][j]+2*uu[i][j-1]-gg(h)*ff(i*h,j*h))*0.25;
    			}
    			else if(i==10 && j<10)
    			{
    				uu[i][j]=(uu[i-1][j]+pp[i-1][j]+pp[i][j+1]+uu[i][j-1]-gg(h)*ff(i*h,j*h))*0.25;
    			}
    			
    			else
    			{
    				uu[i][j]=(pp[i+1][j]+uu[i-1][j]+pp[i][j+1]+uu[i][j-1]-gg(h)*ff(i*h,j*h))*0.25;
    			}
    			
    		}
    	}
    		//Two for loops to calculate the maximum error on each run through of all the points.
    		for(j=0;j<11;j++)
    	{
    		for(i=j;i<11;i++)
    		{
    		error1 = error2; 
    		cc = uu[i][j];
    		dd = pp[i][j];
    		error2 = cc-dd;
    		if (error2 >= error1 ){
    		     maxerror = error2;
    		   } 
    		   else {
    		   	maxerror = error1;
    		   }
    		}
    	}
    	//printf("%f\n",maxerror);
    	//printf("%f\n",maxerror);
   }
    //Print statement for the numerical solutions of the Gauss-Seidel Method
  for(j=0; j<11; j++)
       {
          for(i=j;i<11;i++)
             {
             	printf("u(%f,%f) = %f\n",i*h,j*h,uu[i][j]);
             }
       } 
       //Print statements for the analytical solutions
       for(j=0;j<11;j++)
    	{
    		for(i=j;i<11;i++)
    		{
    			uuu = (sin(pi*0.05*i))*(sin(pi*0.05*j));
    			//printf("u(%f,%f) = %f\n",i*0.05,j*0.05,uu);
    			//printf("%f\n", uu);
    		}
    	}
}
//Functions defined for use in main

float ff(float x, float y) 
{
	float pi=acos(-1);
	float val;
	val = -2*gg(pi)*sin(pi*x)*sin(pi*y);
	return val;
}
//Function for squaring values
float gg(float x)
{
	float val1;
	val1 = x*x;
	return val1;
}
float f (float x)
{
	float value;
	value = x;
	return value;
}

float g (float x, float y)
{
	float value;
	value = (-9.8/0.17)*x - (6/0.3)*(x-y);
	return value;
}
// Harry Shuttleworth B128662
