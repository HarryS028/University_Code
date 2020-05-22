#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float f(float, float);
float g(float);

int main(int argc, const char * argv[]) 
{

	float u[11][11],p[11][11];
	int i,j,n;
	float h=0.05;
	// This sets all the values in the newly created matrix to equal 0, these will be used as 
	// initial guesses for our surface values. 
    for(i=0; i<11; i++)
       {
          for(j=0;j<11;j++)
             {
             	u[i][j] = 0.00000000;
             	p[i][j] = 0.00000000;
             }
       }
    // Iterations for the Gauss-Seidel method begin. 
    for(n=0; n<350; n++)
    {
    	for(j=0;j<11;j++)
    	{
    		for(i=0;i<11;i++)
    		{
    			p[i][j]=u[i][j];
    		}
    	}
    	for(j=0;j<11;j++)
    	{
    		for(i=j;i<11;i++)
    		{
    			if(j==0)
    			{
    				u[i][j]=0;
    			}
    			else if(i==j && i==10)
    			{
    				u[i][j]=(2*p[i][j-1]+2*u[i][j-1]-g(h)*f(i*h,j*h))*0.25;
    			}
    			else if(i==j && i>0)
    			{
    				u[i][j]=(2*p[i+1][j]+2*u[i][j-1]-g(h)*f(i*h,j*h))*0.25;
    			}
    			else if(i==10 && j<10)
    			{
    				u[i][j]=(u[i-1][j]+p[i-1][j]+p[i][j+1]+u[i][j-1]-g(h)*f(i*h,j*h))*0.25;
    			}
    			
    			else
    			{
    				u[i][j]=(p[i+1][j]+u[i-1][j]+p[i][j+1]+u[i][j-1]-g(h)*f(i*h,j*h))*0.25;
    			}
    			
    		}
    	}
   }
  for(j=0; j<11; j++)
       {
          for(i=j;i<11;i++)
             {
             	printf("u(%f,%f) = %f\n",i*h,j*h,u[i][j]);
             }
       } 
}

float f(float x, float y) 
{
	float pi=acos(-1);
	float val;
	val = -2*g(pi)*sin(pi*x)*sin(pi*y);
	return val;
}
//Function for squaring values
float g(float x)
{
	float val1;
	val1 = x*x;
	return val1;
}
// Harry Shuttleworth B128662
