
#ifdef __unix__
#include <unistd.h>
#endif

#include <math.h>
#include "m_param.h"
#include "gdatastream.h"

TProfil* TProfil::pm;

const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481;

enum volume_code {  /* Codes of volume parameter ??? */
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};
                            
SPP_SETTING pa_ = {
    "GEM-Selektor v2-PSI: Controls and defaults for numeric modules",
    {
        1,  /* PC */  3,     /* PD */   3,   /* PRD */
        1,  /* PSM  */ 144,  /* DP */   15,   /* DW */
        -2, /* DT */  0,     /* PLLG */   1,   /* PE */
        500,   /* IIM */
        1e-30, /* DG */   1e-8,  /* DHB */  1e-12,  /* DS */
        1e-5,  /* DK */  0.01,  /* DF */  1e-12,  /* DFM */
        1e-6,  /* DFYw */  1e-6,  /* DFYaq */    1e-6,  /* DFYid */
        1e-6,  /* DFYr,*/  1e-6,  /* DFYh,*/   1e-6,  /* DFYc,*/
        1e-12, /* DFYs, */  1e-17,  /* DB */   0.7,   /* AG */
        0.07,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
        0.001, /* GAS */  12.05,  /* DNS */  1e-5,  /* XwMin, */
        1e-7,  /* ScMin, */    1e-20, /* DcMin, */   1e-10, /* PhMin, */
        3e-5,  /* ICmin */   1e-7,  /* EPS */   1e-15,  /* IEPS */
        1e-4,  /* DKIN  */ 0,  /* tprn */
    },
}; /* SPP_SETTING */


void
BASE_PARAM::write(ostream& oss)
{
    oss.write( (char*)&PC, 10*sizeof(short) );
    oss.write( (char*)&DG, 28*sizeof(double) );
    oss.write( (char*)&tprn, sizeof(char*) );
}


void
SPP_SETTING::write(ostream& oss)
{
    oss.write( ver, TDBVERSION );
    p.write( oss );
}


TProfil::TProfil( int szN, int szM, int szK ):
 sizeN(szN), sizeM(szM), sizeK(szK)
{
    pa= pa_;
    int nNd = nNodes();
    multi = new TMulti( nNd );
    pmp = multi->GetPM();
}

/*-----------------------------------------------------------------*/
// Interpolation over tabulated values (array y) using the Lagrange method
// for extracting thermodynamic data in gemipm2k or in gem2mt
// parameters:
//  y[N] - discrete values of argument over rows (ascending order)
//  x[M] - discrete values of arguments over columns (ascending order)
//  d[N][M] - discrete values of a function of x and y arguments
//  xoi - column (x) argument of interest ( x[0] <= xi <= x[M-1] )
//  yoi - row (y) argument of interest  ( y[0] <= yi <= y[N-1] )
//  N - number of rows in y array;
//  M - number of columns in y array.
//  Function returns an interpolated value of d(yoi,xoi) or 0 if
//  yoi or xoi are out of range
//
double TProfil::LagranInterp(float *y, float *x, double *d, float yoi,
                    float xoi, int M, int N)
{
    double s=0,z,s1[21];
    int py, px, i=0, j, j1, k, jy, jy1;

//    if (yoi < y[0])
//        Error( GetName(), "E34RErun: yoi < y[0] (minimal row argument value)");
//    if(xoi < x[0])
//        Error( GetName(), "E35RErun: xoi < x[0] (minimal column argument value)");
    py = N-1;
    px = M-1;

    if(yoi < y[0] || xoi < x[0] || yoi > y[py] || xoi > x[px] )
       return s;  // one of arguments outside the range

    for(j1=0;j1<N;j1++)
        if (yoi >= y[j1] && yoi <= y[j1+1])
            goto m1;
    //z=yoi;
    goto m2;
m1:
    for(i=0;i<M;i++)
        if(xoi >= x[i] && xoi <= x[i+1])
            goto m;
    // z=xoi;
    if(xoi <= x[px])
        goto m;
m2:
    if(yoi <= y[py])
        goto m;
m:
    if(i < M-px)
        j=i;
    else j=M-px-1;
    if( j1 >= N-py)
        j1=N-py-1;
    jy1=j1;
    for(jy=0;jy <= py; jy++)
    {
        s=0.;
        for(i=0;i<=px;i++)
        {
            z=1; //z1=1;
            for(k=0;k<=px;k++)
                if(k!=i)
                    z*=(xoi-x[k+j])/(x[i+j]-x[k+j]);
            s+=d[i+j+(j1)*M]*z;
        }
        s1[jy]=s;
        j1++;
    }
    s=0.;
    for(i=0;i<=py;i++)
    {
        z=1;
        for(k=0;k<=py;k++)
            if(k!=i)
                z*=(yoi-y[k+jy1])/(y[i+jy1]-y[k+jy1]);
        s+=s1[i]*z;
    }
    return(s);
}


// Load Thermodynamic Data from MTPARM to MULTI using LagranInterp
void TProfil::CompG0Load()
{
  int j, jj, k, jb, je=0;
  double Gg, Vv;
  float TC, P;

  DATACH  *dCH = multi->data_CH;
  DATABR  *dBR = multi->data_BR;

  if( dCH->nTp <=1 && dCH->nPp <=1 )
    return;

  TC = dBR->T-C_to_K;
  P = dBR->P;
  for( jj=0; jj<dCH->nTp; jj++)
    if( fabs( TC - dCH->Tval[jj] ) < dCH->Ttol )
    {
        TC = dCH->Tval[jj];
        break;
    }
  for( jj=0; jj<dCH->nPp; jj++)
   if( fabs( P - dCH->Pval[jj] ) < dCH->Ptol )
   {
        P = dCH->Pval[jj];
        break;
   }

//Test outpur ***********************************
//  fstream f_log("CompG0Load.txt", ios::out|ios::app );

//  f_log << "TC = " <<  TC << "  P =  " << P << endl;
//Test outpur ***********************************

 if( fabs( pmp->TC - TC ) < 1.e-10 &&
            fabs( pmp->P - P ) < 1.e-10 )
   return;    //T, P not changed


 pmp->T = pmp->Tc = TC + C_to_K;
 pmp->TC = pmp->TCc = TC;
 pmp->P = pmp->Pc = P;
 pmp->denW = LagranInterp( dCH->Pval, dCH->Tval, dCH->roW,
                          P, TC, dCH->nTp, dCH->nPp );
 //       pmp->denWg = tpp->RoV;
 pmp->epsW = LagranInterp( dCH->Pval, dCH->Tval, dCH->epsW,
                          P, TC, dCH->nTp, dCH->nPp );
 //       pmp->epsWg = tpp->EpsV;

//Test outpur ***********************************
//  f_log << "roW = " <<  pmp->denW << "  epsW =  " << pmp->epsW << endl;
//Test outpur ***********************************

 pmp->RT = R_CONSTANT * pmp->Tc;
 pmp->FRT = F_CONSTANT/pmp->RT;
 pmp->lnP = 0.;
 if( P != 1. ) /* ??????? */
   pmp->lnP = log( P );

 for( k=0; k<pmp->FI; k++ )
 {
   jb = je;
   je += pmp->L1[k];
   /*load t/d data from DC */
    for( j=jb; j<je; j++ )
    {
      jj =  j * dCH->nPp * dCH->nTp;
      Gg = LagranInterp( dCH->Pval, dCH->Tval, dCH->G0+jj,
                          P, TC, dCH->nTp, dCH->nPp );
      pmp->G0[j] = Cj_init_calc( Gg, j, k );
//Test outpur ***********************************
//  f_log << j  << " Gg = " <<  Gg  << "  GO =  " << pmp->G0[j] << endl;
//Test outpur ***********************************
    }
 }
 for( j=0; j<pmp->L; j++ )
 {
    jj =  j * dCH->nPp * dCH->nTp;
    Vv = LagranInterp( dCH->Pval, dCH->Tval, dCH->V0+jj,
                          P, TC, dCH->nTp, dCH->nPp );
    switch( pmp->PV )
    { /* make mol volumes of components */
       case VOL_CONSTR:
                    pmp->A[j*pmp->N] = Vv;
       case VOL_CALC:
       case VOL_UNDEF:
                    pmp->Vol[j] = Vv  * 10.;
                    break;
    }
//Test outpur ***********************************
//  f_log << j  << " Vv = " <<  Vv  << "  VO =  " << pmp->Vol[j] << endl;
//Test outpur ***********************************
 }
}


// GEM IPM calculation of equilibrium state in MULTI
void
TProfil::calcMulti()
{
    // MultiCalcInit( keyp.c_str() );
    //    // realloc memory for  R and R1
    pmp->R = new double[pmp->N*(pmp->N+1)];
    pmp->R1 = new double[pmp->N*(pmp->N+1)];
    memset( pmp->R, 0, pmp->N*(pmp->N+1)*sizeof(double));
    memset( pmp->R1, 0, pmp->N*(pmp->N+1)*sizeof(double));

    CompG0Load();
    if( AutoInitialApprox() == false )
        MultiCalcIterations();


}

void TProfil::outMulti( GemDataStream& ff, gstring& path  )
{
    ff.writeArray( &pa.p.PC, 10 );
    ff.writeArray( &pa.p.DG, 28 );
    multi->to_file( ff, path );
}

// Reading structure MULTI (GEM IPM work structure)
void TProfil::readMulti( GemDataStream& ff )
{
      ff.readArray( &pa.p.PC, 10 );
      ff.readArray( &pa.p.DG, 28 );
      multi->from_file( ff );
}


// ------------------ End of ms_param.cpp -----------------------




