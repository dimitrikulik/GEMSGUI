#ifndef _m_const_h
#define _m_const_h

// #include <iostream>
#include <ctype.h>
#include <fstream>

//#ifndef IPMGEMPLUGIN
// #define IPMGEMPLUGIN
//#endif

using namespace std;
#include "verror.h"
#include "v_user.h"

const int MST =  6,
          DFCN = 6; // number of columns in MASDJ table

const unsigned int
    MAXICNAME =      6,
    MAXSYMB =        4,
    TDBVERSION =     64;


const int
    MAXDCNAME =      16,
    MAXPHNAME =      16,
    EQ_RKLEN = 58;


enum solmod_switches { /* indexes of keys of model solution*/
    SPHAS_TYP, DCOMP_DEP, SPHAS_DEP, SGM_MODE, DCE_LINK, SCM_TYPE,
    /* link state */
    LINK_UX_MODE, LINK_TP_MODE, LINK_FIA_MODE,
    /* Posible values of �of keys of model solution - DCOMP_DEP, SPHAS_DEP */
    SM_UNDEF = 'N', SM_TPDEP = 'T', SM_UXDEP = 'X', SM_PRIVATE_ = 'P',
    SM_PUBLIC = 'U',
    // Posible modes calculating of activity coefficients SGM_MODE
    SM_STNGAM = 'S', SM_NOSTGAM = 'N',
    /* This code (one upper-case letter or digit) defines type of mixing
    and default method of calculation of mixing properties before and at
    IPM iterations.   Possible values: (SPHAS_TYP) */
    SM_IDEAL =  'I', // ideal solution or single-component phase;
    SM_REDKIS = 'G', // built-in Guggenheim (Redlich-Kister) binary SS model
                     // (with 3 coeffs)
    SM_MARGB = 'M',  // built-in binary Margules SS (subreg. w. 3 coeff U,P,T )
    SM_MARGT = 'T',  // built-in ternary Margules SS (reg. w. 3 coeff U,P,T )
    SM_FLUID = 'F',  // built-in multicomp. EOS fluid Churakov&Gottschalk 2003
    SM_RECIP = 'R',  // reciprocal solution (to be done), formerly regular solution
    SM_AQDAV = 'D',  // built-in Davies equation (with 0.3) added KD 25.01.02
    SM_AQDH1 = '1',  // built-in limiting Debye-Hueckel law for aqueous species
    SM_AQDH2 = '2',  // built-in 2-term Debye-Hueckel (Kielland)
    SM_AQDH3 = '3',  // built-in 3-d approximation of Debye-Hueckel
    SM_AQDHH = 'H',  // built-in 3-d approximation of Debye-Hueckel (Helgeson)
    SM_AQSIT = 'S',  // built-in SIT model for aq activity coeffs (reserved)
    SM_AQPITZ = 'P', // Pitzer's model for aqueous brines (reserved)
    SM_IONEX = 'E',  // ion exchange (Donnan, Nikolskii) (reserved)
    SM_SURCOM = 'A', // models of surface complexation at solid-aqueous interface
    SM_USERDEF = 'U', // user-defined mixing model (in Phase record)
    SM_OTHER = 'O'   //  other models of non-ideal solutions (reserved)
};

#ifndef _chbr_classes_h_

typedef enum {  /* Classifications of DC */
    /* Type of input data for */
    SRC_DCOMP =  'd',  // the key points to existing PDB record in DCOMP chain
    SRC_REACDC = 'r',  // the key points to existing PDB record in REACDC chain
    SRC_NEWDC =  'n',  // the key new reaction-defined component
    SRC_NEWISO = 'i',  // the same as n, but this component is an isotopic form
    SRC_FICT =   'f',  // fictive species
    /*Aqueous electrolyte phase:*/
    DC_AQ_PROTON   = 'T',      // hydrogen ion H+
    DC_AQ_ELECTRON = 'E',      // electron (as a DC)
    DC_AQ_SPECIES  = 'S',      // other aqueous species (ions, complexes and ion pairs)
    DC_AQ_SOLVENT  = 'W',      // water H2O (major solvent)
    DC_AQ_SOLVCOM  = 'L',      // other components of a solvent (eg. alcohol)
    /*Gas phase ( G code can be used for all gases; V,C,H,N codes are reserved
    for future use of the built-in equations of state in FGL module): */
    DC_GAS_COMP    = 'G',   // other gases
    DC_GAS_H2O     = 'V',   // H2O steam
    DC_GAS_CO2     = 'C',   // CO2 (carbon dioxide)
    DC_GAS_H2      = 'H',   // H2 hydrogen
    DC_GAS_N2      = 'N',   // N2 nitrogen
    /* Solid/liquid non-electrolyte multicomponent phases:*/
    DC_SOL_IDEAL   = 'I',   // end-member component with ideal behaviour
    DC_SOL_MINOR   = 'J',   // junior component (Henry's Law)
    DC_SOL_MAJOR   = 'M',   // major component (Raoult's Law)
    /* Sorption phases and poly(oligo)electrolytes */
    DC_SUR_CARRIER = 'Q',   /* Principal end-member of solid carrier */
    DC_SUR_MINAL   = 'P',   /* Minor end-member of solid carrier */
    DC_PEL_CARRIER = 'R',   /* Carrier of poly(oligo)electrolyte */

    DC_SSC_A0 = '0', DC_SSC_A1 = '2', DC_SSC_A2 = '4', DC_SSC_A3 = '6',
    DC_SSC_A4 = '8', /* Strong surface complex on site type 0,1,2,3,4 - A plane */
    DC_WSC_A0 = '1', DC_WSC_A1 = '3', DC_WSC_A2 = '5', DC_WSC_A3 = '7',
    DC_WSC_A4 = '9', /* Weak surface complex on site type 0,1,2,3,4 - B plane */
    DC_IESC_A  = 'A', /* Strong exchange ion const-charge plane */
    DC_IEWC_B  = 'B', /* Weak exchange ion const-charge plane */

    /* Aliaces for 1-site model */
    DC_SUR_GROUP    = 'X',  /* Surface group on A plane -> '0' */
    DC_SUR_COMPLEX = 'Y',  /* Strong sur. complex A plane -> '0' */
    DC_SUR_IPAIR   = 'Z',  /* Weak sur complex B plane -> '1' */

    /* Single-component phases:*/
    DC_SCP_CONDEN  = 'O',   // DC forming a single-component phase

    /* �������⥫�� ���� ��� ���� १���⮢ ���� ࠢ����ᨩ */
    DC_AQA_CATION   = 'c',  /* ��⨮� */
    DC_AQA_ANION    = 'n',  /* ����� */
    DC_AQA_LIGAND   = 'l',  /* �࣠���᪨� ������ */
    DC_AQA_COMPLEX  = 'x',  /* �������� � ������ ��⢮� */
    DC_AQA_IONPAIR  = 'p',  /* ������ ��� */
    DC_AQA_GAS      = 'g',  /* ���⢮७�� ��� (����ࠫ�.) */
    DC_AQA_ACID     = 'a',  /* ��᫮� */
    DC_AQA_BASE     = 'b',  /* ������ */
    DC_AQA_SALT     = 's',  /* ���� (�������. ���஫��) */
    DC_AQA_HYDROX   = 'o'   /* ��� ���பᨫ� */

} DC_CLASSES;


//    This code defines standard state and reference scale of concentra-
// tions for components of this phase. It is used by many subroutines
// during calculations of equilibrium states
enum PH_CLASSES{  /* Possible values */
    PH_AQUEL    = 'a',  // aqueous electrolyte
    PH_GASMIX   = 'g',  // mixture of gases
    PH_FLUID    = 'f',  // fluid phase
    PH_PLASMA   = 'p',  // plasma
    PH_SOLUTION = 'l',  // non-electrolyte liquid (melt)
    PH_SIMELT   = 'm',  // silicate (magmatic) melt or non-aqueous electrolyte
    PH_SORPTION = 'x',  // dilspersed solid with adsorption (ion exchange) in aqueous
    PH_POLYEL = 'y',    // colloidal poly- (oligo)electrolyte
    PH_SINCOND  = 's',  // condenced solid phase, also multicomponent
    PH_SINDIS   = 'd',  // dispersed solid phase, also multicomponent
    PH_HCARBL   = 'h'   // mixture of condensed hydrocarbons
};

#else

//    This code defines standard state and reference scale of concentra-
// tions for components of this phase. It is used by many subroutines
// during calculations of equilibrium states
enum PH_CLASSES2{  /* Possible values */
    PH_PLASMA   = 'p',  // plasma
    PH_SIMELT   = 'm',  // silicate (magmatic) melt or non-aqueous electrolyte
    PH_HCARBL   = 'h'   // mixture of condensed hydrocarbons
};

#endif

typedef enum {  /* Limits on DC and phases */
    /* type of lmits */
    NO_LIM = 'O', LOWER_LIM ='L', UPPER_LIM = 'U', BOTH_LIM ='B',
    /* mode recalc of limits Set_DC_Limits() */
    DC_LIM_INIT = 0, DC_LIM_CURRENT = 1
} DC_LIMITS;


enum sorption_control {
    /* EDL interface models - separate for site types in v. 3.1 */
    SC_DDLM = 'D',  SC_CCM = 'C',     SC_TLM = 'T',   SC_MTL = 'M',
    SC_MXC = 'E',   SC_NNE = 'X',     SC_IEV  = 'I',  SC_BSM = 'S',
    SC_NOT_USED = 'N',
    /* Methods of Surface Activity Terms calculation */
    SAT_COMP = 'C', SAT_NCOMP = 'N', SAT_SITE = 'S', SAT_INDEF = 'I',
/* New methods for surface activity coefficient terms (2004) */
 SAT_L_COMP = 'L', SAT_QCA_NCOMP = 'Q', SAT_QCA1_NCOMP = '1',
 SAT_QCA2_NCOMP = '2', SAT_QCA3_NCOMP = '3', SAT_FREU_NCOMP = 'f',
 SAT_QCA4_NCOMP = '4', SAT_BET_NCOMP = 'B', SAT_VIR_NCOMP = 'W',
 SAT_FRUM_NCOMP = 'F', SAT_FRUM_COMP = 'R', SAT_PIVO_NCOMP = 'P',
    /* Assignment of surtype to carrier (end-member) */
    CCA_VOL = 'V', CCA_0 = '0', CCA_1, CCA_2, CCA_3, CCA_4, CCA_5,
    CCA_6, CCA_7, CCA_8, CCA_9, SPL_0='0', SPL_B, SPL_C, SPL_D,
    SDU_N = 'n', SDU_m = 'm', SDU_M = 'M', SDU_g = 'g',
    CST_0 = '0', CST_1, CST_2, CST_3, CST_4, CST_5, // surface type index
    CSI_0 = '0', CSI_1, CSI_2, CSI_3, CSI_4, CSI_5, // surface site index
// Number of parameters per surface species in the MaSdj array
// MCAS = 6 = DFCN ; position index    added by KD 25.10.2004
// Column indices of surface species allocation table MCAS
   SA_MCA=0, SA_EMX, SA_STX, SA_PLAX, SA_SITX, SA_UNIT,
// Column indices of MaSdj table of adsorption parameters
   PI_DEN=0, PI_CD0, PI_CDB, PI_P1, PI_P2, PI_P3
};

typedef enum { // Units of measurement of quantities and concentrations
    /* number of components and phases */
    QUAN_MKMOL = 'Y',  QUAN_MMOL = 'h',  QUAN_MOL = 'M', // NUMBER OF MOLES
    QUAN_MGRAM = 'y',  QUAN_GRAM = 'g',  QUAN_KILO = 'G',// MASS
    /* concentrations of components and phases*/
    CON_MOLFR = 'n', CON_MOLPROC = 'N', CON_pMOLFR = 'f', // MOLE FRACTION
    CON_VOLFR = 'v', CON_VOLPROC = 'V', CON_pVOLFR = 'u', // VOLUME FRACTION
    CON_MOLAL = 'm', CON_MMOLAL =  'i', CON_pMOLAL = 'p', // MOLALITY
    CON_MOLAR = 'L', CON_MMOLAR =  'j', CON_pMOLAR = 'q', // MOLARITY
    CON_WTFR  = 'w', CON_WTPROC =  '%', CON_PPM =    'P', // MASS FRACTION
    CON_AQWFR = 'C', CON_AQWPROC = 'c', CON_AQPPM =  'a', // CONCENTRATION
    // aqueous species
    CON_AQGPL = 'd', CON_AQMGPL = 'e', CON_AQMKGPL = 'b',//VOLUME CONCENTRATION

    /* �����孮���� ����稭 */
    SUR_AREAM = 's',  SUR_AREAN = 'S', /* ��.�����孮��� �2/�, �2/��� */
    SUR_CONA  = 'A', /* ������. ����⮪ ���/�2 */
    SUR_CONM  = 'D', /* ������. ����⮪ ���/� */
    SUR_CONN  = 'E', /* ������. ����⮪ ���/��� */

    //Units of measurement of pressure Pr, P  { b B p P A }'
    PVT_BAR =  'b', // bar - default, 1 bar = 1e5 Pa
    PVT_KBAR = 'B', // kbar, 1 kbar = 1000 bar
    PVT_PASC = 'p', // Pascal (Pa)
    PVT_KPASC = 'P',// MPa, 1 MPa = 10 bar = 1e6 Pa
    PVT_ATM =  'A', // atm, 1 atm = 1.013 bar
    //Attention: Only b code can be used in this version!

    //Units of measurement of molar volume  { c j a L m }'
    PVT_CM3 =  'c', /*cm3, cm3/mole*/
    PVT_LITR =  'L', // liters (L) - volume of the system only, 1 L = 1000 cm3
    PVT_JBAR =  'j', // J/bar, 10 cm3/mole = 1 J/bar
    PVT_CBAR = 'a',  // (cal/bar), 41.84 cm3/mole = 1 cal/bar
    // m  - reserved.
    //Attention: only j code can be used in this version!

    //Units of measurement of reference temperature Tr { C K F }'
    PVT_CELS = 'C', /*degrees Celsius (C)*/
    PVT_KELVIN = 'K', /*Kelvins (K), 0 C = 273.15 K*/
    PVT_FAREN = 'F',  /*degrees Fahrenheit (F)*/
    //Attention: Only C code can be used in this version.

    // Units of measurement of energy values { j c J C n N }
    TDAT_JOUL = 'j',/* Joules (J/mole)*/ TDAT_KJOUL = 'J', /*kilojoules (kJ/mole)*/
    TDAT_CAL = 'c', /* calories (cal/mole); 1 cal = 4.184 J; */
    TDAT_KCAL = 'C', /*kilocalories (kcal/mole)*/
    TDAT_NORM = 'N' /*normalized (mole/mole, J/mole/RT, cal/mole/RT)*/
                // Attention: Only j code can be used in this version!
} SPPA_UNITS;

const char S_OFF = '-',
                   S_ON = '+',
                          S_REM = '*',
                                  A_NUL ='?';

#endif
