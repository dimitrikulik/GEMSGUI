//-------------------------------------------------------------------
// TNode class - implements a
// simple C/C++ interface between GEM IPM and FMT codes
// Works with DATACH and work DATABR structures
// without using the nodearray class
//
// Written by S.Dmytriyeva, D.Kulik
// Copyright (C) 2006 S.Dmytriyeva, D.Kulik
//
// This file is part of GEMIPM2K and GEMS-PSI codes for
// thermodynamic modelling by Gibbs energy minimization

// This file may be distributed under the licence terms
// defined in GEMIPM2K.QAL
//
// See also http://les.web.psi.ch/Software/GEMS-PSI
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#ifndef _node_h_
#define _node_h_

#include "m_param.h"
#include "datach.h"
#include "databr.h"

#ifndef IPMGEMPLUGIN
class QWidget;
#endif

class TNode
{
    MULTI* pmm;   // Pointer to GEMIPM work data structure (TMULTI in ms_multi.h)

#ifdef IPMGEMPLUGIN
                  // This is used in isolated GEMIPM2K module for coupled codes
    TMulti* multi;
    TProfil* profil;

#endif

protected:

    DATACH* CSD;     // Pointer to chemical system data structure CSD

    DATABR* CNode;   // Pointer to a work node data bridge structure   CNode
                     // used for sending input data to and receiving results from GEM IPM

    void allocMemory();
    void freeMemory();

    // Functions that maintain DATACH and DATABR memory structures
    void datach_realloc();
    void datach_free();
    void databr_realloc();

    DATABR* databr_free( DATABR* data_BR_ =0);  // deletes fields of DATABR structure indicated by data_BR_
                                                // and sets the pointer data_BR_ to NULL
    // Binary i/o functions
    // including file i/o using GemDataStream class (for binary data of different endianness)
    void datach_to_file( GemDataStream& ff );   // writes CSD (DATACH structure) to a binary file
    void datach_from_file( GemDataStream& ff ); // reads CSD (DATACH structure) from a binary file
    void databr_to_file( GemDataStream& ff );   // writes node (work DATABR structure) to a binary file
    void databr_from_file( GemDataStream& ff ); // reads node (work DATABR structure) from a binary file

    // Text i/o functions
    void datach_to_text_file( fstream& ff );    // writes CSD (DATACH structure) to a text file
    void datach_from_text_file( fstream& ff);   // reads CSD (DATACH structure) from a text file
    void databr_to_text_file(fstream& ff );     // writes work node (DATABR structure) to a text file
    void databr_from_text_file(fstream& ff );   // reads work node (DATABR structure) from a text file

    // virtual functions for interaction with nodearray class
    virtual void  setNodeArray( int , int*  ) { }
    virtual void  checkNodeArray( int, int*, const char* ) { }

    virtual int nNodes()  const   // virtual call for interaction with nodearray class
    { return 1; }

#ifndef IPMGEMPLUGIN
    // Integration in GEMS - prepares DATACH and DATABR files for reading into the coupled code
    void makeStartDataChBR(
         TCIntArray& selIC, TCIntArray& selDC, TCIntArray& selPH,
         short nTp_, short nPp_, float Ttol_, float Ptol_,
         float *Tai, float *Pai );

    void getG0_V0_H0_Cp0_matrix();  // creates arrays of thermodynamic data for interpolation
                                    // which are written into DATACH file

    // virtual function for interaction with nodearray class
    virtual void  setNodeArray( gstring& , int , bool ) { }
#endif

public:

#ifndef IPMGEMPLUGIN
   TNode( MULTI *apm );   // constructor for integration in GEMS
#else

  static TNode* na;   // static pointer to this class
                      // for the isolated GEMIPM2K module
  TNode();      // constructor for GEMIPM2K
#endif

    virtual ~TNode();      // destructor

// call sequence ----------------------------------------------


    // For separate coupled FMT-GEM programs that use GEMIPM2K module
    // Reads in the MULTI, DATACH and DATABR files prepared from GEMS
    // and fills out nodes in node arrays according to distribution vector
    // nodeTypes ( only for TNodeArray )
    int  GEM_init( const char *ipmfiles_lst_name,
                   int *nodeTypes = 0, bool getNodT1 = false);

#ifdef IPMGEMPLUGIN
//  Wrapper calls for direct coupling of an FMT code with GEMIPM2K

   // Loads GEM input data from FMT part provided in parameters
   // into the DATABR work structure for the subsequent GEM calculation
   void GEM_from_MT(
       short  p_NodeHandle,   // Node identification handle
       short  p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                        //                                     GEM input output  FMT control
       double p_T,      // Temperature T, K                        +       -      -
       double p_P,      // Pressure P, bar                         +       -      -
       double p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
       double p_Ms,     // Mass of reactive subsystem, kg          -       -      +
//       double p_dt,     // actual time step	         			  +
//       double p_dt1,    // previous time step                                   +
       double *p_dul,   // upper kinetic restrictions [nDCb]       +       -      -
       double *p_dll,   // lower kinetic restrictions [nDCb]       +       -      -
       double *p_bIC    // bulk mole amounts of IC [nICb]          +       -      -
   );

   // Copies GEM input data from already loaded DATABR work structure
   // into parameters provided by the FMT part
   void GEM_restore_MT(
   short  &p_NodeHandle,   // Node identification handle
   short  &p_NodeStatusCH, // Node status code;  see typedef NODECODECH
                    //                                     GEM input output  FMT control
   double &p_T,      // Temperature T, K                        +       -      -
   double &p_P,      // Pressure P, bar                         +       -      -
   double &p_Vs,     // Volume V of reactive subsystem, cm3     -       -      +
   double &p_Ms,     // Mass of reactive subsystem, kg          -       -      +
//       double p_dt,     // actual time step	         			  +
//       double p_dt1,    // previous time step                                   +
   double *p_dul,   // upper kinetic restrictions [nDCb]       +       -      -
   double *p_dll,   // lower kinetic restrictions [nDCb]       +       -      -
   double *p_bIC    // bulk mole amounts of IC [nICb]          +       -      -
   );

#endif

   // Main call for GEM IPM calculation
   int  GEM_run();   // calls GEM for a work node

   // For examining GEM calculation results:
   // Prints current multi and/or work node structures to files with
   // names given in the parameter list (if any parameter is NULL
   // then writing the respective file is skipped)
   void  GEM_printf( const char* multi_file, const char* databr_text,
                         const char* databr_bin );

#ifdef IPMGEMPLUGIN
//  Wrapper calls for direct coupling of an FMT code with GEMIPM2K

   // Copies GEM calculation results into parameters provided by the
   // FMT part (dimensions and order of elements in arrays must correspond
   // to the DATACH structure )
   void GEM_to_MT(
       short &p_NodeHandle,    // Node identification handle
       short &p_NodeStatusCH,  // Node status code (changed after GEM calculation); see typedef NODECODECH
       short &p_IterDone,      // Number of iterations performed by GEM IPM
                         //                                     GEM input output  FMT control
       // Chemical scalar variables
       double &p_Vs,    // Volume V of reactive subsystem, cm3     -      -      +     +
       double &p_Ms,    // Mass of reactive subsystem, kg          -      -      +     +
       double &p_Gs,    // Gibbs energy of reactive subsystem (J)  -      -      +     +
       double &p_Hs,    // Enthalpy of reactive subsystem (J)      -      -      +     +
       double &p_IC,    // Effective molal aq ionic strength       -      -      +     +
       double &p_pH,    // pH of aqueous solution                  -      -      +     +
       double &p_pe,    // pe of aqueous solution                  -      -      +     +
       double &p_Eh,    // Eh of aqueous solution, V               -      -      +     +
       double &p_denW,
       double &p_denWg, // Density of H2O(l) and steam at T,P      -      -      +     +
       double &p_epsW,
       double &p_epsWg, // Diel.const. of H2O(l) and steam at T,P  -      -      +     +
       // Dynamic data - dimensions see in DATACH.H and DATAMT.H structures
       // exchange of values occurs through lists of indices, e.g. xDC, xPH
       double  *p_xDC,    // DC mole amounts at equilibrium [nDCb]      -      -      +     +
       double  *p_gam,    // activity coeffs of DC [nDCb]               -      -      +     +
       double  *p_xPH,  // total mole amounts of phases [nPHb]          -      -      +     +
       double  *p_vPS,  // phase volume, cm3/mol        [nPSb]          -      -      +     +
       double  *p_mPS,  // phase (carrier) mass, g      [nPSb]          -      -      +     +
       double  *p_bPS,  // bulk compositions of phases  [nPSb][nICb]    -      -      +     +
       double  *p_xPA,  // amount of carrier in phases  [nPSb] ??       -      -      +     +
       double  *p_rMB,  // MB Residuals from GEM IPM [nICb]             -      -      +     +
       double  *p_uIC  // IC chemical potentials (mol/mol)[nICb]       -      -      +     +
   );

#endif

    DATACH* pCSD() const  // get pointer to chemical system data structure
    {     return CSD;   }

    DATABR* pCNode() const  // get pointer to work node data structure
    {        return CNode;
    }  // usage on the level of nodearray is not recommended !

    // These methods get contents of fields in the work node structure
    double cT() const     // get current Temperature T, K
    {        return CNode->T;   }

    double cP() const     // get current Pressure P, bar
    {        return CNode->P;   }

    void setNodeHandle( int jj )   // setup Node identification handle
    {      CNode->NodeHandle = (short)jj;  }

   // Return DCH index of IC by Name or -1 if illegal name
   int IC_name_to_x( const char *Name );
   // Return DCH index of DC by Name or -1 if illegal name
   int DC_name_to_x( const char *Name );
   // Return DCH index of Ph by Name or -1 if illegal name
   int Ph_name_to_x( const char *Name );

   // Return DBR index of IC by Name or -1 if illegal name
   int IC_name_to_xDB( const char *Name )
   { return IC_xCH_to_xDB( IC_name_to_x( Name ) ); }
   // Return DBR index of DC by Name or -1 if illegal name
   int DC_name_to_xDB( const char *Name )
   { return DC_xCH_to_xDB( DC_name_to_x( Name ) ); }
   // Return DBR index of Ph by Name or -1 if illegal name
   int Ph_name_to_xDB( const char *Name )
   { return Ph_xCH_to_xDB( Ph_name_to_x( Name ) ); }

   // Return for IComp DBR index from DCH index
   // or -1 if not used in the data bridge
   int IC_xCH_to_xDB( const int xCH );
   // Return for DComp DBR index from DCH index
   // or -1 if not used in the data bridge
   int DC_xCH_to_xDB( const int xCH );
   // Return for Phase DBR index from DCH index
   // or -1 if not used in the data bridge
   int Ph_xCH_to_xDB( const int xCH );

   // Return for IComp DCH index from DBR index
   int IC_xBR_to_xCH( const int xBR )
   { return CSD->xIC[xBR]; }
   // Return for DComp DCH index from DBR index
   int DC_xBR_to_xCH( const int xBR )
   { return CSD->xDC[xBR]; }
   // Return for Phase DCH index from DBR index
   int Ph_xBR_to_xCH( const int xBR )
   { return CSD->xPH[xBR]; }

    // Data exchange methods between GEMIPM and work node DATABR structure
    void packDataBr();      //  packs GEMIPM calculation results into work node structure
    void unpackDataBr();    //  unpacks work node structure into GEMIPM data structure


#ifndef IPMGEMPLUGIN
// These calls are used only inside GEMS-PSI GEM2MT module

    // Makes start DATACH and DATABR data using GEMS internal data (MULTI and other)
    // interaction variant (the user must select ICs, DCs and phases to be included
    // in DATABR lists)
    void MakeNodeStructures( QWidget* par, bool select_all,
             float *Tai, float *Pai, short nTp_ = 1 ,
             short nPp_ = 1 , float Ttol_ = 1., float Ptol_ =1. );

    // Overloaded variant - takes lists of ICs, DCs and phases according to
    // already existing index vectors axIC, axDC, axPH (with anICb, anDCb,
    // anPHb, respectively)
    void MakeNodeStructures(  short anICb, short anDCb,  short anPHb,
             short* axIC, short* axDC,  short* axPH,
             float* Tai, float* Pai,  short nTp_,
             short nPp_, float Ttol_, float Ptol_  );

#endif

};

#endif   // _node_h_
