// this is the cpp function that implements the filter. It looks quite messy, but it
// will work for any value of k. Basically I just broke the matrix operations
// into different steps.

#include "mex.h"
#include "nr3matlab.h"

void doCalc(int k,MatDoub p, const MatDoub n,mwSize nr, MatDoub &f_out,MatDoub &E_out)
{
    
    int i;
    double Dk=n.nrows();    // k in double type
    for (i=0;i<k;i++)       // Inicialize first filtered probabilities
    {
        E_out[0][i]=1/Dk;
    }
    
    double sumN=0,sumP;
    
    int j;
    int t;
    
    MatDoub B1(k,1);    // this is Coeff.p*E(i-1,:)'.*n(i,:)'
    MatDoub B2(k,1);    // this is ones(k,1)'*(Coeff.p*E(i-1,:)'.*n(i,:)')
    
    for (t=1;t<nr;t++)
    {
        for (i=0;i<k;i++)
        {
            sumN=0;
            sumP=0;
            for (j=0;j<k;j++)
            {
                sumP=p[j][i]*E_out[t-1][j];
                sumN=sumN+sumP;
            }
            
            B1[i][0]=sumN;
            
        }
        sumN=0;
        for (i=0;i<k;i++)
        {
            B2[i][0]=(B1[i][0])*n[i][t];
            
            sumN=sumN+B2[i][0];
        }
        
        f_out[t][0]=sumN;
        
        
        
        for (i=0;i<k;i++)
        {
            E_out[t][i]=B2[i][0]/f_out[t][0];
            
        }
        
    }
    
}


void mexFunction( int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[] )
{
    int k;
    
    mwSize nr;
    nr = mxGetM(prhs[1]);
    k  = mxGetN(prhs[1]);
    
    // Building inputs with NR interface
    
    const MatDoub n(prhs[1]);
    const MatDoub p(prhs[0]);
    
    // Building outputs with NR interface
    
    MatDoub f_out(nr,1,plhs[0]);
    MatDoub E_out(nr,k,plhs[1]);
    
    doCalc(k,p,n,nr,f_out,E_out);   // function for Hamilton's Filter
    
}