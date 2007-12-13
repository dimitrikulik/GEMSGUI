//-------------------------------------------------------------------
// $Id$
//
// Declaration of TMulti class, configuration, and related functions
// based on the IPM work data structure MULTI that represents chemical
// thermodynamic multisystem
//
// Rewritten from C to C++ by S.Dmytriyeva
// Copyright (C) 1995-2007 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and of the
// standalone GEMIPM2K code (define IPMGEMPLUGIN).
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//
#ifndef _ms_multi_h_
#define _ms_multi_h_

#ifndef IPMGEMPLUGIN

#include "m_param.h"
#include "v_ipnc.h"
// Internal subroutine for ET_translate
typedef int (tget_ndx)( int nI, int nO, int Xplace );

#else
#include <time.h>
#include "m_const.h"

#endif

typedef struct
{  // MULTI is base structure to Project (local values)
  char
    stkey[EQ_RKLEN+1],   // Record key identifying IPM minimization problem
    // NV_[MAXNV], nulch, nulch1, // Variant Nr for fixed b,P,T,V; index in a megasystem
    PunE,         // Units of energy  { j;  J c C N reserved }
    PunV,         // Units of volume  { j;  c L a reserved }
    PunP,         // Units of pressure  { b;  B p P A reserved }
    PunT;         // Units of temperature  { C; K F reserved }
  short
    N,        	// N - number of IC in IPM problem
    NR,       	// NR - dimensions of R matrix
    L,        	// L -   number of DC in IPM problem
    Ls,       	// Ls -   total number of DC in multi-component phases
    LO,       	// LO -   index of water-solvent in IPM DC list
    PG,       	// PG -   number of DC in gas phase
    PSOL,     	// PSOL - number of DC in liquid hydrocarbon phase
    Lads,     	// Lads - number of DC in sorption phases
    FI,       	// FI -   number of phases in IPM problem
    FIs,      	// FIs -   number of multicomponent phases
    FIa,      	// FIa -   number of sorption phases
    FI1,     // FI1 -   number of phases present in eqstate
    FI1s,    // FI1s -   number of multicomponent phases present in eqstate
    FI1a,    // FI1a -   number of sorption phases present in eqstate
    IT,      // It - number of completed IPM iterations
    E,       // PE - flag of electroneutrality constraint { 0 1 }
    PD,      // PD - mode of calling GammaCalc() { 0 1 2 3 4 }
    PV,      // PV - flag of system volume constraint { 0 1 }
    PLIM,    // PU - flag of activation of DC/phase restrictions { 0 1 }
    Ec,    // GammaCalc() return code: 0 (OK) or 1 (error)
    K2,    // Number of Selekt2() loops
    PZ,    // Indicator of IPM-2 precision algorithm activation
//     funT, sysT,
    pNP, //Mode of FIA selection: 0-auto-SIMPLEX,1-old eqstate,-1-user's choice
    pESU,  // Unpack old eqstate from EQSTAT record?  0-no 1-yes
    pIPN,  // State of IPN-arrays:  0-create; 1-available; -1 remake
    pBAL,  // State of reloading CSD:  1- BAL only; 0-whole CSD
    pFAG_,  // reserved SD
    pTPD,  // State of reloading thermod data: 0- all  1 - G0 only  2 - no
    pULR,  // Start recalc kinetic constraints (0-do not, 1-do )internal
    pRR1,  // internal
    FIat,   // max. number of surface site types
    MK,     // PM return code: 0 - continue;  1 - converged
    W1,     // internal IPM-2 indicator
    is,     // is - index of IC for IPN equations ( GammaCalc() )
    js,     // js - index of DC for IPN equations ( GammaCalc() )
    next,
    sitNcat,    // SIT: number of cations
    sitNan      // SIT: number of anions
    ;
  double
    TC,TCc, 	// Temperature T, min.-max. (0,2000 C)
    T,Tc,   	// T, min.-max. K
    P,Pc,   	// Pressure P, min.-max.(0,10000 bar)
    VX_,VXc,    // V(X) - volume of the system, min.-max., cm3
    GX_,GXc,    // Gibbs potential of the system G(X), min.-max. (J)
    AX_,AXc,    // Helmholtz potential of the system F(X), reserved
    UX_,UXc,  	// Internal energy of the system U(X), reserved
    HX_,HXc, 	// Total enthalpy of the system H(X), reserved
    SX_,SXc, 	// Total entropy of the system S(X), reserved
    CpX_,CpXc,  // reserved
    CvX_,CvXc,  // reserved
    T0,         // reserved
    VE,         // reserved
    MBX,        // Total mass of the system, kg
    FX,    	// Current Gibbs potential of the system in IPM, moles
    IC,         // Effective molal ionic strength of aqueous electrolyte
    pH,         // pH of aqueous solution
    pe,         // pe of aqueous solution
    Eh,         // Eh of aqueous solution, V
    DHBM,       // Adjusted balance precision criterion (IPM-2 )
    DSM,        // min value phase DS (IPM-2)
    GWAT,       // used in ipm_gamma()
    YMET,       // reserved
    denW,denWg,  // Density of H2O(l) and steam for Tc,Pc
    epsW,epsWg,  // Diel. constant of H2O(l) and steam for Tc,Pc
    PCI,        // Current value of Dikin criterion of IPM convergence DK>=DX
    DX,         // IPM convergence criterion threshold DX (1e-5)
    lnP,        // log Ptotal
    RT,         // RT: 8.31451*T (J/mole/K)
    FRT,        // F/RT, F - Faraday constant = 96485.309 C/mol
    Yw,         // Current number of moles of solvent in aqueous phase
    ln5551,     // ln(55.508373) = 4.0165339
    aqsTail,    // v_j asymmetry correction factor for aqueous species
    lowPosNum,  // Minimum DC amount defining the Ls set (1e-19)
    logXw,      // work variable
    logYFk,     // work variable
    YFk,        // Current number of moles in a multicomponent phase
    FitVar[5];  // internal; FitVar[0] is T,P-dependent b_gamma parameter
  short
    *L1,    // l_a vector - number of DCs included into each phase [Fi]
    *LsMod, // Number of interaction parameters, max. parameter order (cols in IPx),
        // and number of coefficients per parameter in PMc table [3*FIs]
    *LsMdc, // Number of non-ideality coeffs per one DC in multicomponent phase[FIs]
    *IPx,   // List of indexes of interaction parameters for non-ideal solutions
            // ->LsMod[k,0] x LsMod[k,1]   added 07.12.2006   KD
    *mui,   // IC indices in RMULTS IC list [N]
    *muk,   // Phase indices in RMULTS phase list [FI]
    *muj;   // DC indices in RMULTS DC list [L]
  short  (*SATX)[4]; // New: work table [Lads]: link indexes to surface type [XL_ST];
            // sorbent em [XL_EM]; surf.site [XL-SI] and EDL plane [XL_SP]
  float
    *PMc,    // Non-ideality coefficients f(TP) -> LsMod[k,0] x LsMod[k,2]
    *DMc,    // Non-ideality coefficients f(TPX) for DC -> LsMdc[k]
    *A,      // DC stoichiometry matrix A composed of a_ji [0:N-1][0:L-1]
    *Awt,    // IC atomic (molar) mass, g/mole [0:N-1]
    *Wb,     //Relative Born factors (HKF, reserved) [0:Ls-1]
    *Wabs,   // Absolute Born factors (HKF, reserved) [0:Ls-1]
    *Rion,   // Ionic or solvation radii, A (reserved) [0:Ls-1]
    *HYM,    // reserved
    *ENT,    // reserved no object
    *H0,     // DC p-molar enthalpies, reserved [L]
    *A0,     // DC p-molar Helmholtz energies, reserved [L]
    *U0,     // DC p-molar internal energies, reserved [L]
    *S0,     // DC p-molar entropies, reserved [L]
    *Cp0,    // DC p-molar entropies, reserved [L]
    *Cv0,    // DC p-molar Cv, reserved [L]
    *VL,     // ln mole fraction of end members in phases-solutions
    *Xcond, 	// conductivity of phase carrier, sm/m2   [0:FI-1], reserved
    *Xeps,  	// diel.permeability of phase carrier (solvent) [0:FI-1], reserved
    *Aalp,  	// phase specific surface area m2/g       [0:FI-1]
    *Sigw,  	// st.surface free energy in water,J/m2   [0:FI-1]
    *Sigg;  	// st.surface free energy in air,J/m2     [0:FI-1], reserved

//  Data for surface comlexation and sorption models (new variant [Kulik,2006])
  float  (*Xr0h0)[2];   // mean r & h of particles (- pores), nm  [0:FI-1][2], reserved
  float  (*Nfsp)[MST];  // area fraction of surface types At/A  [FIs][FIat]
  float  (*MASDT)[MST]; // Max. surface species density Gamma_C, mkmol/g [FIs][FIat]
  float  (*XcapF)[MST]; // Capacitance density of Ba EDL layer F/m2 [FIs][FIat]
  float  (*XcapA)[MST]; // Capacitance density of 0 EDL layer, F/m2 [FIs][FIat]
  float  (*XcapB)[MST]; // Capacitance density of B EDL layer, F/m2 [FIs][FIat]
  float  (*XcapD)[MST]; // Eff. cap. density of diffuse layer, F/m2 [FIs][FIat]
  float  (*XdlA)[MST];  // Effective thickness of A EDL layer, nm [FIs][FIat], reserved
  float  (*XdlB)[MST];  // Effective thickness of B EDL layer, nm [FIs][FIat], reserved
  float  (*XdlD)[MST];  // Effective thickness of diffuse layer, nm [FIs][FIat], reserved
  float  (*XlamA)[MST]; // Factor of EDL discretness  A < 1 [FIs][FIat], reserved
  float  (*Xetaf)[MST]; // Permanent charge density at surface type, C/m2 [FIs][FIat]
  float  (*MASDJ)[DFCN];  // Max. density, CD-music and isotherm params [Lads][DFCN]
                          // Contents defined in the enum below this structure
// Other data
  float  *XFs,  // Current quantities of phases X_a at IPM iterations [0:FI-1]
    *Falps,     // Current Karpov criteria of phase stability  F_a [0:FI-1]
    *Fug,       // Demo partial fugacities of gases [0:PG-1]
    *Fug_l,     // Demo log partial fugacities of gases [0:PG-1]
    *Ppg_l;     // Demo log partial pressures of gases [0:PG-1]
  double
    *DUL,     // VG Vector of upper kinetic restrictions to x_j, moles [L]
    *DLL,     // NG Vector of lower kinetic restrictions to x_j, moles [L]
    *GEX,     // Increments to molar G0 values of DCs, normalized [L]
    *PUL,  // Vector of upper restrictions to phases amounts X_a (reserved)[FIs]
    *PLL,  // Vector of lower restrictions to phases amounts X_a (reserved)[FIs]
    *YOF,     // Phase metastability parameter (spec.surf.energy), in J/g [FI !!!!]
    *Vol,     // DC molar volumes, cm3/mol [L]
    *MM,      // DC molar masses, g/mol [L]
    *Pparc,   // DC partial pressures/ pure fugacities, bar (Pc by default) [0:L-1]
    *Y_m,     // Molalities of aqueous species and sorbates [0:Ls-1]
    *Y_la,    // log activity of DC in multi-component phases[0:Ls-1]
    *Y_w,     // Mass concentrations of DC in multi-component phases,%(ppm)[Ls]
    *Gamma,   // DC activity coefficients [0:L-1]
    *lnGmf,   // ln of initial DC activity coefficients [0:L-1]
    *lnGmM,   // ln of DC pure gas fugacity (or metastability) coefficients [0:L-1]
    *EZ,      // Formula charge of DC in multi-component phases [0:Ls-1]
    *FVOL,    // phase volumes, cm3/mol                   [0:FI-1]
    *FWGT,    // phase (carrier) masses, g                [0:FI-1]
//
    *G,    // Normalized DC energy function c(j), mole/mole [0:L-1]
    *G0,   // Input normalized g0_j(T,P) for DC at unified standard scale[L]
    *lnGam, // ln of DC activity coefficients [0:L-1]
    *lnGmo; // Copy of lnGam from previous IPM iteration (reserved)
  double  (*lnSAC)[4]; // former lnSAT ln surface activity coeff and Coulomb's term  [Lads][4]
  double  *B,  // Input bulk chem. compos. of the system - b vector, moles of IC[N]
    *U,  // IC chemical potentials u_i (mole/mole) - dual IPM solution [N]
    *U_r,  // IC chemical potentials u_i (J/mole) [0:N-1]
    *C,    // Calculated IC mass-balance deviations (moles) [0:N-1]
    *IC_m, // Total IC molalities in aqueous phase (excl.solvent) [0:N-1]
    *IC_lm,	// log total IC molalities in aqueous phase [0:N-1]
    *IC_wm,	// Total dissolved IC concentrations in g/kg_soln [0:N-1]
    *BF,    //Output bulk compositions of multicomponent phases bf_ai[FIs][N]
    *BFC,   //Total output bulk compositions of all solid phases[1][N]
    *XF,    // Output total number of moles of phases Xa[0:FI-1]
    *YF,    // Copy of X_a from previous IPM iteration [0:FI-1]
    *XFA,   // Quantity of carrier in asymmetric phases Xwa, moles [FIs]
    *YFA,   // Copy of Xwa from previous IPM iteration [0:FIs-1]
    *Falp;  // Karpov phase stability criteria F_a [0:FI-1]
// EDL models (data for electrostatic activity coefficients)
   double (*XetaA)[MST]; // Total EDL charge on A (0) EDL plane, moles [FIs][FIat]
   double (*XetaB)[MST]; // Total charge of surface species on B (1) EDL plane, moles[FIs][FIat]
   double (*XetaD)[MST]; // Total charge of surface species on D (2) EDL plane, moles[FIs][FIat]
   double (*XpsiA)[MST]; // Relative potential at A (0) EDL plane,V [FIs][FIat]
   double (*XpsiB)[MST]; // Relative potential at B (1) EDL plane,V [FIs][FIat]
   double (*XpsiD)[MST]; // Relative potential at D (2) plane,V [FIs][FIat]
   double (*XFTS)[MST];  // Total number of moles of surface DC at surface type [FIs][FIat]
//
   double *X,  // DC quantities at eqstate x_j, moles - primal IPM solution [L]
    *Y,   // Copy of x_j from previous IPM iteration [0:L-1]
    *XY,  // Copy of x_j from previous loop of Selekt2() [0:L-1]
    *Qp,  // Work variables related to non-ideal phases FIs*(QPSIZE=60)
    *Qd,  // Work variables related to DC in non-ideal phases FIs*(QDSIZE=60)
    *MU,  // mu_j values of differences between dual DC chem.potentials [L]
    *EMU, // Exponents of DC increment to F_a criterion for phase [L]
    *NMU, // DC increments to F_a criterion for phase [L]
    *W,   // Weight multipliers for DC (incl restrictions) in IPM [L]
    *Fx,  // Dual DC chemical potentials defined via u_i and a_ji [L]
    *Wx,  // Mole fractions Wx of DC in multi-component phases [L]
    *F,   //Primal DC chemical potentials defined via g0_j, Wx_j and lnGam_j[L]
    *F0;  // Excess Gibbs energies for (metastable) DC, mole/mole [L]
   double (*D)[MST];  // Reserved; new work array for calc. surface act.coeff.
// Name lists
  char  (*sMod)[6];   // Codes of models of multicomponent phases [FIs]
  char  (*SB)[MAXICNAME+MAXSYMB]; // List of IC names in the system [N]
  char  (*SB1)[MAXICNAME]; // List of IC names in the system [N]
  char  (*SM)[MAXDCNAME];  // List of DC names in the system [L]
  char  (*SF)[MAXPHNAME+MAXSYMB];  // List of phase names in the system [FI]
  char  (*SM2)[MAXDCNAME];  // List of multicomp. phase DC names in the system [Ls]
  char  (*SM3)[MAXDCNAME];  // List of adsorption DC names in the system [Lads]
  char  *DCC3;   // Classifier of DC in sorption phases [Lads]
  char  (*SF2)[MAXPHNAME+MAXSYMB]; // List of multicomp. phase names in the syst [FIs]
  char  (*SFs)[MAXPHNAME+MAXSYMB];
    // List of phases currently present in non-zero quantities [FI]
  char  *pbuf, 	// Text buffer for table printouts
// Class codes
    *RLC,   // Classifier of restriction types for x_j [L]
    *RSC,   // Classifier of restriction scales for x_j [L]
    *RFLC,  // Classifier of restriction types for XF_a [FIs]
    *RFSC,  // Classifier of restriction scales for XF_a [FIs]
    *ICC,   // Classifier of IC { e o h a z v i <int> } [N]
    *DCC,   // Classifier of DC { TESWGVCHNIJMDRAB0123XYZPQO } [L]
    *PHC;   // Classifier of phases { a g f p m l x d h } [FI]
  char  (*SCM)[MST]; // Classifier of adsorption models for sur types [FIs][FIat]
  char  *SATT,  // classifier of methods of SAT calculation [Lads]
    *DCCW;  // internal DC class codes [L]
//  short
//     *sitXcat, // SIT: indices of cations (may be changed soon)
//     *sitXan;  // SIT: indices of anions
//  float
//     *sitE;    // pointer to SIT coeff. table (may be changed soon)
  short ITF,        // Number of completed EFD iterations
        ITG;        // Number of completed GEM IPM iterations
  clock_t t_start, t_end;
  double t_elap_sec;  // work variables for determining IPM calculation time
#ifdef IPMGEMPLUGIN
  float *Guns;  //  mu.L work vector of uncertainty space increments to tp->G + sy->GEX
  float *Vuns;  //  mu.L work vector of uncertainty space increments to tp->Vm
  double *tpp_G; // Partial molar(molal) Gibbs energy g(TP) (always), J/mole 
  float *tpp_S;    // Partial molar(molal) entropy s(TP), J/mole/K
  double *tpp_Vm;   // Partial molar(molal) volume Vm(TP) (always), J/bar
#endif  
}
MULTI;

enum { // link indexes to surface type [XL_ST] sorbent em [XL_EM]
//  surf.site [XL-SI] and EDL plane [XL_SP]
   XL_ST = 0, XL_EM, XL_SI, XL_SP
};

// Data of MULTI
class TMulti
#ifndef IPMGEMPLUGIN
  : public TSubModule
#endif
{
    MULTI pm;
    MULTI *pmp;

// Internal arrays for the optimization part (since version 2.0.0)
   int sizeN;
   double *AA;
   double *BB;
   int *arrL;
   int *arrAN;

   void Alloc_A_B( int newN );
   void Free_A_B();
   void Build_compressed_xAN();
   void Free_compressed_xAN();
   void Free_internal();

#ifndef IPMGEMPLUGIN
// This stuff is used only in GEMS-PSI
    SYSTEM *syp;
    MTPARM *tpp;
    RMULTS *mup;

    void multi_sys_dc();
    void multi_sys_ph();
    void ph_sur_param( int k, int kk );
    void ph_surtype_assign( int k, int kk, int jb, int je,
                            short car_l[], int car_c, short Cjs );
    void sm_text_analyze( int nph, int Type, int JB, int JE, int jb, int je );
    void SolModLoad();
//    float *PackSITcoeffs( int k, int JB, int JE, int jb, int je, int nCxA );
    gstring PressSolMod( int nP );
    char *ExtractEG( char *Etext, int jp, int *EGlen, int Nes );
    int find_icnum( char *name, int LNmode );
    int find_dcnum( char *name, int jb, int je, int LNmode );
    int find_phnum( char *name, int LNmode );
    int find_acnum( char *name, int LNmode );
#else

   char PAalp_;
   char PSigm_;
   float EpsW_;
   float RoW_;

#endif

// Internal functions for SCMs
   void getLsModsum( int& LsModSum, int& LsIPxSum );
   void getLsMdcsum( int& LsMdcSum );

// ipm_chemical.cpp
    void XmaxSAT_IPM2();
    void XmaxSAT_IPM2_reset();
    double DualChemPot( double U[], float AL[], int N, int j );
    void Set_DC_limits( int Mode );
    void TotalPhases( double X[], double XF[], double XFA[] );
    double Ej_init_calc( double, int j, int k);
    double  PrimalDC_ChemPot(  double G,  double logY,  double logYF,
                           double asTail,  double logYw,  char DCCW );
    void PrimalChemicalPotentials( double F[], double Y[],
                                  double YF[], double YFA[] );
    double KarpovCriterionDC( double *dNuG, double logYF, double asTail,
                 double logYw, double Wx,  char DCCW );
    void f_alpha();
    double FreeEnergyIncr(   double G,  double x,  double logXF,
                             double logXw,  char DCCW );
    double GX( double LM  );
    double Cj_init_calc( double g0, int j, int k );
    void Mol_u( double Y[], double X[], double XF[], double XFA[] );
    void ConvertDCC();
    int  getXvolume();

// ipm_chemical2.cpp
    void GasParcP();
    void phase_bcs( int N, int M, int jb, float *A, double X[], double BF[] );
    void phase_bfc( int k, int jj );
    double pH_via_hydroxyl( double x[], double Factor, int j);
    void ConCalcDC( double X[], double XF[], double XFA[],
                    double Factor, double MMC, double Dsur, int jb, int je, int k );
    void ConCalc( double X[], double XF[], double XFA[]);
    void GouyChapman(  int jb, int je, int k );
//  Surface activity coefficient terms
    void SurfaceActivityCoeff( int jb, int je, int jpb, int jdb, int k );
    void SurfaceActivityTerm( int jb, int je, int k );  // Obsolete / deleted

// ipm_chemical2.cpp
    void IS_EtaCalc();
    void pm_GC_ods_link( int k, int jb, int jpb, int jdb, int ipb );
    double TinkleSupressFactor( double ag, int ir);
// Main call for calculation of activity coefficients on IPM iterations
    void GammaCalc( int LinkMode );
// Built-in activity coefficient models
//  aqueous electrolyte
    void DebyeHueckel3Hel( int jb, int je, int jpb, int jdb, int k );
    void DebyeHueckel3Karp( int jb, int je, int jpb, int jdb, int k );
    void DebyeHueckel2Kjel( int jb, int je, int jpb, int jdb, int k );
    void DebyeHueckel1LL( int jb, int je, int k );
    void Davies03temp( int jb, int je, int k );
    void SIT_aqac_PSI( int jb, int je, int jpb, int jdb, int k, int ipb );
// fluid mixtures
    void ChurakovFluid( int jb, int je, int jpb, int jdb, int k );
    void CGofPureGases( int jb, int je, int jpb, int jdb, int k );
    void PRSVFluid( int jb, int je, int jpb, int jdb, int k, int ipb );
    void PRSVofPureGases( int jb, int je, int jpb, int jdb, int k, int ipb );
// condensed mixtures
    void RedlichKister( int jb, int je, int jpb, int jdb, int k );
    void MargulesBinary( int jb, int je, int jpb, int jdb, int k );
    void MargulesTernary( int jb, int je, int jpb, int jdb, int k );
// Generic solution model calls
void SolModParPT ( int jb, int je, int jpb, int jdb, int k, int ipb, char ModCode );
void SolModActCoeff( int jb, int je, int jpb, int jdb, int k, int ipb, char ModCode );

// ipm_main.cpp - numerical part of GEM-IPM2
    void MultiCalcMain();
    int EnterFeasibleDomain( );
    int InteriorPointsMethod( );
    void SimplexInitialApproximation( );

// ipm_main.cpp - miscellaneous fuctions of GEM-IPM2
   void MassBalanceResiduals( int N, int L, float *A, double *Y,
                               double *B, double *C );
   double LMD( double LM );
   void ZeroDCsOff( int jStart, int jEnd, int k=-1 );
   void RaiseZeroedOffDCs( int jStart, int jEnd, double sfactor, int k=-1 );
   void LagrangeMultiplier();
   void WeightMultipliers( bool square );
   int SolverLinearEquations( int N, bool initAppr );
   double calcDikin(  int N, bool initAppr );
   double calcLM(  bool initAppr );
   void Restoring_Y_YF();
   double calcSfactor();
   int PhaseSelect( int &k_miss, int &k_unst );
   // IPM_SIMPLEX.CPP Simplex method with two side constraints (Karpov ea 1997)
    void Simplex(int M, int N, int T, double GZ, double EPS,
                 double *UND, double *UP, double *B, double *U,
                 double *AA, int *STR, int *NMB );
    void SPOS( double *P, int STR[],int NMB[],int J,int M,double AA[]);
    void START( int T,int *ITER,int M,int N,int NMB[],
                double GZ,double EPS,int STR[],int *BASE,
                double B[],double UND[],double UP[],double AA[],double *A,
                double *Q );
    void NEW(int *OPT,int N,int M,double EPS,double *LEVEL,int *J0,
             int *Z,int STR[], int NMB[], double UP[],
             double AA[], double *A);
    void WORK(double GZ,double EPS,int *I0, int *J0,int *Z,int *ITER,
              int M, int STR[],int NMB[],double AA[],
              int BASE[],int *UNO,double UP[],double *A,double Q[]);
    void FIN(double EPS,int M,int N,int STR[],int NMB[],
             int BASE[],double UND[],double UP[],double U[],
             double AA[],double *A,double Q[],int *ITER);

public:

#ifndef IPMGEMPLUGIN
// This is used only in GEMS-PSI
    TIArray<IPNCalc> qEp;
    TIArray<IPNCalc> qEd;

    TMulti( int nrt, SYSTEM* sy_, MTPARM *tp_, RMULTS *mu_ );
    ~TMulti()
    {  Free_internal(); };


    void ods_link( int i=0);
    void dyn_set( int i=0);
    void dyn_kill( int i=0);
    void dyn_new( int i=0);
    void set_def( int i=0);
//    void sit_dyn_new();

    // ms_muleq.cpp
    void packData();
    void packData( TCIntArray PHon, TCIntArray DCon );
    void setSizes();
    void loadData( bool newRec );
    void unpackData();

    void EqstatExpand( const char *key );
    void MultiRemake( const char *key );
    void ET_translate( int nOet, int nOpex, int JB, int JE, int jb, int je,
     tget_ndx *get_ndx = 0 );

   class UserCancelException {};
#else
// this variant is used only in GEMIPM2K
   TMulti()
   { pmp = &pm;
     sizeN = 0;
     AA = 0;
     BB = 0;
     arrL = 0;
     arrAN = 0;
     pmp->Guns = 0;
     pmp->Vuns = 0;
     pmp->tpp_G = 0;  
     pmp->tpp_S = 0; 
     pmp->tpp_Vm = 0;
     
   }

    void multi_realloc( char PAalp, char PSigm );
    void multi_free();

#endif

    MULTI* GetPM()
    { return &pm; }

    const char* GetName() const
    {  return "Multi";  }

   //mass transport
    void to_file( GemDataStream& ff, gstring& path  );
    void to_text_file( const char *path );
    void from_file( GemDataStream& ff );
    void to_text_file_gemipm( const char *path, bool addMui, bool with_comments = true );
    void from_text_file_gemipm( const char *path );

    // EXTERNAL FUNCTIONS
    // MultiCalc
    void Alloc_internal();
    void MultiCalcInit( const char *key );
    bool AutoInitialApprox();
    void MultiCalcIterations();
    void CompG0Load();

    // connection to Unspace
    double pb_GX( double *Gxx  );
};

void inArray( fstream& ff, char *name, float* arr, int size );
void inArray( fstream& ff, char *name, short* arr, int size );
void inArray( fstream& ff, char *name, double* arr, int size );
void inArray( fstream& ff, char *name, char* arr,
                              int size, int arr_siz );

/*
void outArray( fstream& ff, char *name, short* arr,
                            int size, int l_size=-1 );
void outArray( fstream& ff, char *name,  float* arr,
                             int size, int l_size=-1 );
void outArray( fstream& ff, char *name,  double* arr,
                             int size, int l_size=-1 );
void outArray( fstream& ff, char *name, char* arr,
                              int size, int arr_siz );
*/
#endif   //_ms_multi_h

