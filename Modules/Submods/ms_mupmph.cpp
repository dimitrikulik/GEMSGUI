//-------------------------------------------------------------------
// $Id: ms_mupmph.cpp 1379 2009-07-29 14:53:22Z gems $
//
// (c) 1992-2007 S.Dmitrieva, D.Kulik, K.Chudnenko, I.Karpov
//
// Implementation of parts of the Interior Points Method (IPM) module
// for convex programming Gibbs energy minimization, described in:
// (Karpov, Chudnenko, Kulik (1997): American Journal of Science
//  v.297 p. 767-806)
//
// This file is part of a GEM-Selektor (GEMS) v.3.1.x program
// environment for thermodynamic modeling in geochemistry
// Uses: GEM-Selektor GUI GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the GPL v.3 license

//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include "m_param.h"
#include "s_formula.h"
#include "v_mod.h"

// Aseembling base GEM-IPM structure to calculate equilibrium states
//
// 1) Makes full  stoichiometry matrix A using the sy.Dcl list
// 2) Looks through Phase records and loads data into mixing model arrays
//    (optionally translating non-ideality Phase scripts)
// 3) Continues looking through one-component Phases, to seek for
//    specific surfaces and related surface parameters
// 4) Performs necessary re-scaling and conversions of data
// 5) Tests whether it is possible to take old primal (x, gamma) and dual
//   (u) solutions as an initial approximation
//

// Load System data
// Aseembling base GEM-IPM structure to calculate equilibrium states
//
// 1) Makes full  stoichiometry matrix A using the sy.Dcl list
// 2) Looks through Phase records and loads data into mixing model arrays
//    (optionally translating non-ideality Phase scripts)
// 3) Continues looking through one-component Phases, to seek for
//    specific surfaces and related surface parameters
void TMulti::MultiSystemInit() {
  int i, N;
  short ii;
  RMULTS *mup = TRMults::sm->GetMU();
  SYSTEM *syp = TSyst::sm->GetSY();

  pm.N = syp->N;
  pm.NR = syp->N;
  pm.L = syp->L;
  pm.Ls = syp->Ls;
  pm.LO = syp->Lw /*-1*/;
  pm.PG = syp->Lg;
  pm.PSOL = syp->Lhc;
  pm.Lads = syp->Lsor;
  pm.FI = syp->Fi;
  pm.FIs = syp->Fis;
  pm.FIa = syp->Fia; // new

  if (syp->DLLim == S_ON || syp->DULim == S_ON)
    pm.PLIM = 1;
  else
    pm.PLIM = 0;

  if (syp->PE != S_OFF)
    pm.E = 1;
  else
    pm.E = 0;

  if (syp->NsTm > 0 && syp->PNfsp != S_OFF)
    pm.FIat = MNST; // SCM Modernization: MNST=6
  else
    pm.FIat = 0;

  pm.MBX = syp->MBX;
  pm.GWAT = syp->Mwat;
  if (syp->Mbel > 1e-6) // workaround 21.06.13 DK - needs more investigation!
    pm.IC = syp->Mbel;
  pm.GWAT = syp->Mwat * H2O_mol_to_kg; // constant corrected 30.08.2008

  // Reallocating memory, if necessary
  if (!pm.pBAL) {
    pmp->pKMM = -1;
    pmp->ITau = -1;
    pmp->kTau = 0.;
    pmp->kdT = 0.;
    dyn_new();
  }
  if (pm.pBAL == 2)
    goto NEXT2;

  // loading parameters for ICs (independent components)
  for (N = 0, i = -1, ii = 0; ii < mup->N; ii++) {
    if (syp->Icl[ii] == S_OFF)
      continue;
    i++;
    N++;
    pm.mui[i] = ii;
    pm.ICC[i] = mup->ICC[ii];
    pm.Awt[i] = mup->BC[ii];
    pm.B[i] = syp->B[ii];
    pm.U[i] = pm.U_r[i] = pm.C[i] = 0.0;
    memcpy(pm.SB[i], mup->SB[ii], MAXICNAME + MAXSYMB);
    memcpy(pm.SB1[i], mup->SB[ii], MAXICNAME);
  }

  if (N != pm.N)
    vfMessage(window(), "Multi make error: N != pm.N",
              "Please, press the BCC button first!");
NEXT2:
  // loading data for DCs (dependent components)
  multi_sys_dc();

  // loading data for phases, incl. sorption and solution models
  multi_sys_ph();

  if (pm.pBAL < 2)
    ConvertDCC(); // Loading generic species codes

  if (!pm.pBAL && pm.pTPD > 1)
    pm.pTPD = 1; // DC_ LoadThermodynamicData();

  // Tests on integrity of CSD can be added here
  pm.pBAL = 2;
}

// Load System data to define lookup arrays
void TMulti::SystemToLookup() {
  SYSTEM *syp = TSyst::sm->GetSY();

  // copy intervals for minimizatiom
  pm.Pai[0] = syp->Pmin;
  pm.Pai[1] = syp->Pmax;
  pm.Pai[2] = syp->Pinc;
  pm.Pai[3] = 0.1;
  pm.Tai[0] = syp->Tmin;
  pm.Tai[1] = syp->Tmax;
  pm.Tai[2] = syp->Tinc;
  pm.Tai[3] = 0.1;
  pm.Fdev1[0] = 0.;
  pm.Fdev1[1] = syp->Tdev1;
  pm.Fdev2[0] = 0.;
  pm.Fdev2[1] = syp->Tdev2;

  if (pm.tMin == G_TP_) // set default intervals
  {
    if (pm.Pai[0] > pm.Pc || pm.Pai[1] < pm.Pc) {
      pm.Pai[0] = pm.Pai[1] = pm.Pc;
      pm.Pai[2] = 0.;
    }
    if (pm.Tai[0] > pm.TC || pm.Tai[1] < pm.TC) {
      pm.Tai[0] = pm.Tai[1] = pm.TC;
      pm.Tai[2] = 0.;
    }
  }
}

// Load data for DC from Modelling Project definition to MULTI structure
//
void TMulti::multi_sys_dc() {
  int j, ii, L, iZ = 0;
  short jj, jja, ja, kk;
  double a, *A; // , Vv =0.;
  double mm;
  std::vector<TFormula> aFo;
  std::string form;
  RMULTS *mup = TRMults::sm->GetMU();
  SYSTEM *syp = TSyst::sm->GetSY();

  ErrorIf(!TMTparm::sm->GetTP()->G, "Multi", "Multi make error: !tpp->G");
  ///    int xVol = 0;   // SD 09/02/2007
  ///    if( pm.PV == VOL_CONSTR )
  ///      xVol = getXvolume();

  if (!pm.pBAL) { // making full stoichiometry matrix from DC formula list
    // load formulae
    for (ii = 0; ii < mup->L; ii++) {
      aFo.push_back(TFormula());
      form = aFo[ii].form_extr(ii, mup->L, mup->DCF);
      aFo[ii].SetFormula(form.c_str()); // and ce_fscan
    }

    A = new double[mup->N * mup->L];
    fillValue(A, 0., (mup->N * mup->L));
    for (ii = 0; ii < mup->L; ii++)
      aFo[ii].Stm_line(mup->N, A + ii * mup->N, (char *)mup->SB, mup->Val);
    aFo.clear();
    decomp_descr.clear();
    // loading data for dependent components
    for (L = 0, j = -1, jj = 0; jj < mup->L; jj++) {
      if (syp->Dcl[jj] == S_OFF)
        continue;
      j++;
      L++;
      pm.DCC[j] = mup->DCC[jj];
      pm.muj[j] = jj;
      pm.lnGam[j] = 0.0;
      pm.lnGmo[j] = 0.0;
      pm.X[j] = 0.0;
      pm.Y[j] = 0.0;
      pm.XY[j] = 0.0;
      memcpy(pm.SM[j], mup->SM[jj] + MAXSYMB + MAXDRGROUP, MAXDCNAME);
      decomp_descr.push_back(TRMults::sm->component_description(jj));
      if (j < syp->Ls)
        memcpy(pm.SM2[j], mup->SM[jj] + MAXSYMB + MAXDRGROUP, MAXDCNAME);
      if (j < syp->Ls &&
          j >= syp->Ls - syp->Lsor) { // assembling DC name list for sorption
                                      // surface species
        ja = (j - (syp->Ls - syp->Lsor));
        if (!(pm.SM3 &&
              pm.DCC3)) { // Here a workaround of crash detected on 7.02.2011
          char buf[32];
          memcpy(buf, mup->SM[jj], MAXSYMB + MAXDRGROUP + MAXDCNAME);
          buf[MAXSYMB + MAXDRGROUP + MAXDCNAME] = 0;
          vfMessage(window(), "Multi make warning: Lsor != pm.Lads, DC:", buf);
        } else {
          memcpy(pm.SM3[ja], mup->SM[jj] + MAXSYMB + MAXDRGROUP, MAXDCNAME);
          // surface species DC class code (usage will be obsolete)
          pm.DCC3[ja] = mup->DCC[jj];
        }
      }
      copyValues(pm.dcMod[j], TMTparm::sm->GetTP()->DCMod[jj], 6);
      mm = 0.0;
      for (ii = 0; ii < pm.N; ii++) // compressing the stoichiometry matrix
      {
        a = A[jj * mup->N + pm.mui[ii]];
        pm.A[j * pm.N + ii] = a;
        mm += a * pm.Awt[ii]; // calculating molar mass of DC
      }
      pm.MM[j] = mm;
    }
    delete[] A;
    // A = 0;
    //     ErrorIf( L != pm.L, "Multi", "Multi make error: L != pm.L" );
    if (L != pm.L)
      vfMessage(window(), "Multi make error: L != pm.L", "Please, press BCC!");
  }

  if (pm.E) { // index of charge
    for (ii = 0; ii < pm.N; ii++)
      if (pm.ICC[ii] == IC_CHARGE)
        goto CH_FOUND;
    // error no charge
    if (vfQuestion(window(), "Multi",
                   "No DC with formula charge are included:\n"
                   "Proceed (Y), Cancel (N)?")) {
      pm.E = 0;
      syp->PE = S_OFF; // Change Sveta 14/10/2002
    } else
      Error("Multi", "No DC with formula charge are included");
  CH_FOUND:
    iZ = ii;
  }

  for (j = -1, jj = 0; jj < mup->L; jj++) // Main loop for data loading
  {
    if (syp->Dcl[jj] == S_OFF)
      continue;
    j++;
    pm.DCC[j] = mup->DCC[jj];
    if (pm.E)
      pm.EZ[j] = pm.A[pm.N * j + iZ];
    pm.Pparc[j] = pm.P; // SD  15/05/22 pm.Pc
    pm.F0[j] = 0.0;
    pm.XY[j] = 0.0;

    ///        if( tpp->PtvVm == S_ON )
    ///            switch( pm.PV )
    ///            { // calculating actual molar volume of DC
    ///            case VOL_CONSTR:
    ///              if( syp->Vuns )
    ///                       Vv = syp->Vuns[jj];
    ///                pm.A[j*pm.N+xVol] = tpp->Vm[jj] + Vv;
    ///            case VOL_CALC:
    ///            case VOL_UNDEF:
    ///                if( syp->Vuns )
    ///                       Vv = syp->Vuns[jj];
    ///                pm.Vol[j] = ( tpp->Vm[jj] + Vv )* 10.; // Check!
    ///               break;
    ///            }
    ///        else pm.Vol[j] = 0.0;

    // testing kinetic/metastability constraints
    if (pm.PLIM && jj < mup->L) { // KM constraints necessary!
      pm.RLC[j] = syp->RLC[jj];
      pm.RSC[j] = syp->RSC[jj];
      if (syp->DUL)
        pm.DUL[j] = syp->DUL[jj];
      else
        pm.DUL[j] = 1e6;
      ja = (j - (pm.Ls - pm.Lads));
      if (pm.lnSAC && ja < pm.Lads && ja >= 0 &&
          pm.DUL[j] < 1e6 - TProfil::pm->pa.p.DKIN)
        pm.lnSAC[ja][3] = pm.DUL[j]; // Copy of DUL for SACT refining
      if (syp->DLL)
        pm.DLL[j] = syp->DLL[jj];
      else
        pm.DLL[j] = 0.0;
    } else {
      pm.RLC[j] = NO_LIM;
      pm.RSC[j] = QUAN_MOL;
      pm.DUL[j] = 1e6; // tc_DCmaxm trivial restrictions!
      pm.DLL[j] = 0.0;
    }
    // additional data
    pm.fDQF[j] = 0.0; // DQF parameters or pure-gas fugacities zeroed off

    if (syp->PParc != S_OFF)
      pm.Pparc[j] = syp->Pparc[jj];
    pm.lnGmM[j] = 0.0;
    if (jj < mup->L && syp->PlnGf != S_OFF) // was jj < mup->Ls
    {
      pm.lnGmf[j] = syp->lnGmf[jj];
    } else
      pm.lnGmf[j] = 0.0;

    if (!(j < pm.Ls && j >= pm.Ls - pm.Lads))
      continue; // the following is not done for non-surface species

    if (!pm.SATC || !pm.MASDJ || pm.Lads <= 0)
      continue;

    ja = (j - (pm.Ls - pm.Lads));
    jja = (jj - (mup->Ls - mup->Lads));
    // SCM Modernization: Loading MASDJ and mCDc for surface species
    if (syp->PMaSdj != S_OFF) {
      for (kk = 0; kk < DNFCN; kk++)
        if (!IsFloatEmpty(syp->MASDJ[jja][kk]))
          pm.MASDJ[ja][kk] = syp->MASDJ[jja][kk];
        else
          pm.MASDJ[ja][kk] = 0.0;

      if (syp->XSMD && pm.XSMD) {
        for (kk = 0; kk < NSURFSITE; kk++)
          pm.XSMD[ja][kk] = syp->XSMD[jja][kk] - '0';
      }

      if (syp->mCDc && pm.mCDc) {
        for (kk = 0; kk < DMCDC; kk++)
          if (!IsFloatEmpty(syp->mCDc[jja][kk]))
            pm.mCDc[ja][kk] = syp->mCDc[jja][kk];
          else
            pm.mCDc[ja][kk] = 0.0;
      }
    }
    if (syp->PSATT !=
        S_OFF) { // Loading SATC - codes of methods for SAT calculation
      pm.SATT[ja] = syp->SATC[jja][SURSP_INDIV_SACT_CODE];

      // Loading allocation to surface types
      switch (syp->SATC[jja][SURSP_ALLOC_SURFTYPE]) {
      default:
      case CST_0:
        pm.SATC[ja][XL_ST] = 0;
        break;
      case CST_1:
        pm.SATC[ja][XL_ST] = 1;
        break;
      case CST_2:
        pm.SATC[ja][XL_ST] = 2;
        break;
      case CST_3:
        pm.SATC[ja][XL_ST] = 3;
        break;
      case CST_4:
        pm.SATC[ja][XL_ST] = 4;
        break;
      case CST_5:
        pm.SATC[ja][XL_ST] = 5;
        break;
      }

      // Loading allocation to carrier end-members
      switch (syp->SATC[jja][SURSP_ALLOC_SORBENT_EM]) {
      case CCA_0:
        pm.SATC[ja][XL_EM] = 0;
        break;
      case CCA_1:
        pm.SATC[ja][XL_EM] = 1;
        break;
      case CCA_2:
        pm.SATC[ja][XL_EM] = 2;
        break;
      case CCA_3:
        pm.SATC[ja][XL_EM] = 3;
        break;
      case CCA_4:
        pm.SATC[ja][XL_EM] = 4;
        break;
      case CCA_5:
        pm.SATC[ja][XL_EM] = 5;
        break;
      case CCA_VOL:
        pm.SATC[ja][XL_EM] = -1; // whole sorbent
        break;
      default:
        pm.SATC[ja][XL_EM] = -2; // main DC of the sorbent
        break;
      }

      // SCM Modernization: Loading allocation to surface sites from XSMD
      pm.SATC[ja][XL_SI] = 0;
      if (syp->XSMD) {
        for (int s = 0; s < NSURFSITE; s++) {
          if (syp->XSMD[jja][s] > 0) {
            pm.SATC[ja][XL_SI] = s;
            break;
          }
        }
      }
      pm.SATC[ja][XL_SP] = 0; // Obsolete: plane charges now handled via mCDc
    }
  }
}

//==========================================
bool TMulti::CompressPhaseIpxt(int kPH) {
  int jj, jb, cnt = 0;
  TCIntArray aDCused;
  TPhase *aPH = dynamic_cast<TPhase *>(aMod[RT_PHASE].get());
  RMULTS *mup = TRMults::sm->GetMU();
  TCStringArray form_array;
  TFormula aFo;

  for (jj = 0, jb = 0; jj < kPH; jj++)
    jb += mup->Ll[jj];

  for (jj = 0, cnt = 0; jj < mup->Ll[kPH]; jj++) {
    if (TSyst::sm->GetSY()->Dcl[jj + jb] == S_OFF)
      aDCused.push_back(-1);
    else {
      aDCused.push_back(cnt);
      cnt++;
      form_array.push_back(aFo.form_extr(jj + jb, mup->L, mup->DCF));
    }
  }

  if (cnt < mup->Ll[kPH])
    return aPH->CompressRecord(cnt, aDCused, form_array, true);
  else
    return true;
}

//============================================

// Loading data for phases and mixing models into the structure MULTI
//
void TMulti::multi_sys_ph() {
  int k, i;
  bool non_sorption_phase, is_ss;
  short kk, j, je, jb;
  time_t crt;
  //    double G;
  double PMM;                    // Phase mean mol. mass
  int Cjs, car_l[32], car_c = 0; // current index carrier sorbent
  TPhase *aPH = dynamic_cast<TPhase *>(aMod[RT_PHASE].get());
  RMULTS *mup = TRMults::sm->GetMU();
  SYSTEM *syp = TSyst::sm->GetSY();

  aPH->ods_link(0);

  pm.Ls = 0;
  pm.PG = 0;
  pm.FIs = syp->Fis;
  pm.FIa = syp->Fia;
  pm.Lads = 0; // Is this really needed?
               //    ja=0;

  int JB, JE = 0;
  phase_descr.clear();
  for (kk = 0, /*FI=0,*/ k = -1, /*jb=0,*/ je = 0; kk < mup->Fi; kk++) {
    JB = JE;
    JE += mup->Ll[kk];
    if (TSyst::sm->GetSY()->Pcl[kk] == S_OFF)
      continue;
    k++; // FI++;  //PMM = 0.0;
    pm.L1[k] = syp->Ll[kk];
    jb = je;
    je += pm.L1[k];
    pm.muk[k] = kk;
    pm.PHC[k] = mup->PHC[kk];
    memcpy(pm.SF[k], mup->SF[kk], MAXSYMB);
    memcpy(pm.SF[k] + MAXSYMB, mup->SF[kk] + MAXSYMB + MAXPHSYMB, MAXPHNAME);
    if (kk < mup->Fis) // bugfix 4.04.14 DK was k < syp.Fis
    {
      memcpy(pm.SF2[k], mup->SF[kk], MAXSYMB);
      memcpy(pm.SF2[k] + MAXSYMB, mup->SF[kk] + MAXSYMB + MAXPHSYMB, MAXPHNAME);
    }
    if (pm.FIa > 0 && pm.SF3) {
      int ka = k - (pm.FIs - pm.FIa);
      if (ka >= 0 && ka < pm.FIa) {
        memcpy(pm.SF3[ka], mup->SF[kk], MAXSYMB);
        memcpy(pm.SF3[ka] + MAXSYMB, mup->SF[kk] + MAXSYMB + MAXPHSYMB,
               MAXPHNAME);
      }
    }
    if (pm.pNP && pm.pBAL && pm.pKMM > 0)
      goto PARLOAD;

    aPH->TryRecInp(mup->SF[kk], crt, 0); // Now reading all phase records!
    phase_descr.push_back(char_array_to_string(aPH->php->name, MAXFORMULA) +
                          " (Phase)");

    // New stuff for TKinMet
    // Resetting vector of TSolMod pointers to avoid problems with activity
    // coeffs
    //     in TSolMod calculations after switching phases on/off
    if (pm.ITau < 0) {
      if (phKinMet[k])
        delete phKinMet[k];
      phKinMet[k] = nullptr;
      pm.kMod[k][0] = aPH->php->kin_t[2];
      pm.kMod[k][1] = aPH->php->kin_t[3];
      pm.kMod[k][2] = aPH->php->kin_t[4];
      pm.kMod[k][3] = aPH->php->kin_t[5];
      pm.kMod[k][4] = aPH->php->kin_t[6];
      pm.kMod[k][5] = aPH->php->kin_t[7];
    }

    if (aPH->php->PrpCon == S_ON) {
      pm.LsKin[k * 6] = aPH->php->nPRk;
      pm.LsKin[k * 6 + 1] = aPH->php->nSkr;
      pm.LsKin[k * 6 + 2] = aPH->php->nrpC;
      pm.LsKin[k * 6 + 3] = aPH->php->naptC;
      pm.LsKin[k * 6 + 4] = aPH->php->nAscC;
      pm.LsKin[k * 6 + 5] = aPH->php->nFaces;
    } else {
      pm.LsKin[k * 6] = pm.LsKin[k * 6 + 1] = pm.LsKin[k * 6 + 2] = 0;
      pm.LsKin[k * 6 + 3] = pm.LsKin[k * 6 + 4] = pm.LsKin[k * 6 + 5] = 0;
    }
    /*
    long int
    *xSKrC,  ///< new: Collected array of aq/gas/sorption species indexes used
    in activity products (-> += LsKin[k][1]) *ocPRkC; ///< new: Collected array
    of operation codes for kinetic parallel reaction terms (-> += LsKin[k][0])
    */

    if (kk < mup->Fis &&
        mup->Ll[kk] >
            1) { // getting data and parameters for a multicomponent phase
                 //           if( pm.pNP && pm.pBAL )
                 //                goto PARLOAD;
      // Resetting vector of TSolMod pointers to avoid problems with activity
      // coeffs
      //     in TSolMod calculations after switching phases on/off
      //     (DK 25.05.2009)
      if (phSolMod[k])
        delete phSolMod[k];
      phSolMod[k] = nullptr;
      //          aPH->TryRecInp( mup->SF[kk], crt, 0 );
      // read informations from phase-solution
      memcpy(pm.sMod[k], aPH->php->sol_t, 6);
      pm.sMod[k][6] = aPH->php->kin_t[0];
      pm.sMod[k][7] = aPH->php->kin_t[1];
      // 10/12/2019 added for multi-site mixed moodels
      if (aPH->php->nMoi > 0) {
        TFormula aFo;
        TCStringArray form_array;

        // build formula list
        for (j = JB; j < JE; j++) {
          form_array.push_back(aFo.form_extr(j, mup->L, mup->DCF));
        }

        // get moiety full structure from phase
        aPH->MakeSublatticeLists(form_array);
      }
      // Added SD 20/01/2010
      if (aPH->php->Ppnc == S_ON && aPH->php->npxM > 0)
        CompressPhaseIpxt(kk);
      if (aPH->php->Ppnc == S_ON) {
        pm.LsMod[k * 3] = aPH->php->ncpN;
        pm.LsMod[k * 3 + 2] = aPH->php->ncpM;
        pm.LsMod[k * 3 + 1] = aPH->php->npxM;
      } else
        pm.LsMod[k * 3] = pm.LsMod[k * 3 + 1] = pm.LsMod[k * 3 + 2] = 0;
      if (aPH->php->Psco == S_ON)
        pm.LsMdc[k * 3] = aPH->php->nscM;
      else
        pm.LsMdc[k * 3] = 0;

      pm.LsPhl[k * 2] = aPH->php->nlPh;
      pm.LsPhl[k * 2 + 1] = aPH->php->nlPc;
      pm.LsMdc2[k * 3] = aPH->php->ndqf;
      pm.LsMdc2[k * 3 + 1] = aPH->php->nrcp;

      pm.LsUpt[k * 2] = aPH->php->numpC;
      if (aPH->php->PumpCon && aPH->php->lICu &&
          (aPH->php->kin_t[4] == KM_UPT_ENTRAP ||
           aPH->php->kin_t[4] == KM_IEX_FAST ||
           aPH->php->kin_t[4] == KM_IEX_SLOW))
        pm.LsUpt[k * 2 + 1] =
            pm.L1[k]; // bugfix 04.03.2015 DK
                      //  aPH->php->nFaces;     // provisional 8.04.2013 DK

      // Postponed TSorpMod (branch no_tsorpmod_intermed)
      // if(phSorpMod[k])
      //    delete phSorpMod[k];
      // phSorpMod[k] = nullptr;
      /*
      long int
      *LsESmo, ///< new: number of EIL model layers; EIL params per layer; CD
      coefs per DC; reserved  [Fis][4] *LsISmo, ///< new: number of surface
      sites; isotherm coeffs per site; isotherm coeffs per DC; max.denticity of
      DC [Fis][4] *xSMd;   ///< new: denticity of surface species per surface
      site (site allocation) (-> L1[k]*LsISmo[k][3]] )
      */
    } else // nothing to read in the Phase record
      if (k < pm.FIs) {
        pm.LsMod[k * 3] = 0;
        pm.LsMod[k * 3 + 1] = 0;
        pm.LsMod[k * 3 + 2] = 0;
        pm.LsMdc[k * 3] = 0;
        memset(pm.sMod[k], ' ', 8);
        pm.LsPhl[k * 2] = 0;
        pm.LsPhl[k * 2 + 1] = 0;
        pm.LsMdc2[k * 3] = 0;
        pm.LsMdc2[k * 3 + 1] = 0;
        pm.LsUpt[k * 2] = 0;
        pm.LsUpt[k * 2 + 1] = 0;
        // New stuff for TSorpMod
      }

  PARLOAD:
    if (k < syp->Fis)
      pm.Ls += pm.L1[k];
    non_sorption_phase = true;
    switch (pm.PHC[k]) { // loading necessary parameters
    case PH_AQUEL:
      break;
    case PH_GASMIX:
    case PH_PLASMA:
    case PH_FLUID:
      //            pm.PG = pm.L1[k];
      pm.PG += pm.L1[k];
      break;
    case PH_HCARBL:
      pm.PSOL = pm.L1[k];
      break;
    case PH_LIQUID:
    case PH_SINCOND:
    case PH_SIMELT:
    case PH_IONEX:   // 'i' ClaySor planar IX NEM
    case PH_ADSORPT: // 'q' ClaySor edge NEM
      break;
    case PH_ADSSBEIL:
      //   case PH_ADSORPT:    'q'  ClaySor edge NEM
    case PH_POLYEL:
    case PH_SORPTION:
      non_sorption_phase = false;
      ph_sur_param(k);
      for (i = 0; i < 32; i++)
        car_l[i] = -1;
      car_c = 0; // cQ=0;
      break;
    case PH_SINDIS:
      break;
    default:; // Error message may be needed here!
    }
    // Moved by DK on 25.07.2006
    if (syp->PAalp != S_OFF)
      pm.Aalp[k] = syp->Aalp[kk];
    if (syp->PSigm != S_OFF) {
      pm.Sigw[k] = syp->Sigm[kk][0];
      pm.Sigg[k] = syp->Sigm[kk][1];
    }
    if (syp->PXepsC != S_OFF) {
      pm.Xcond[k] = syp->XEpsC[kk][1];
      pm.Xeps[k] = syp->XEpsC[kk][0];
    }
    if (syp->PXr0h0 != S_OFF) {
      pm.Xr0h0[k][0] = syp->Xr0h0[kk][0];
      pm.Xr0h0[k][1] = syp->Xr0h0[kk][1];
    }
    PMM = 0.0;
    Cjs = -1;
    // Loop over DCs in Phase, calculation of phase mean mol. mass
    for (j = jb; j < je; j++) {
      is_ss = true;
      switch (pm.DCC[j]) {
      case DC_AQ_SOLVENT:
        pm.LO = j;
        is_ss = false;
        PMM = pm.MM[j];
        break;
      case DC_PEL_CARRIER:
      case DC_SUR_CARRIER:
        Cjs = j;
        PMM = pm.MM[j];
        is_ss = true;
        car_l[car_c++] = j;
        pm.Lads++;
        //    ja++;
        break;
      case DC_SUR_MINAL:
        is_ss = true;
        car_l[car_c++] = j;
        pm.Lads++;
        //   ja++;
        break;
        // cleaned up here on 17.10.2026 by DK
      // case DC_IESC_A:
      // case DC_IEWC_B:
      case DC_SUR_SPECIES:
      case DC_SUR_GROUP:
      case DC_SUR_COMPLEX:
      case DC_SUR_IPAIR:
        pm.Lads++;
        //    ja++;
        break; // [[fallthrough]]; ?
      case DC_IESC_A:
      case DC_IEWC_B:
      case DC_SCM_SPECIES:
        is_ss = false; // true;
        PMM += pm.MM[j];
        break;
        /* End extension */
      case DC_SOL_MAJOR:
      case DC_SOL_MAJDEP:
        Cjs = j;
        [[fallthrough]];
      case DC_SOL_MINOR:
      case DC_SOL_MINDEP:
      case DC_SOL_IDEAL:
        PMM += pm.MM[j];
        is_ss = false;
        break;
      default: /* if( isdigit( pm.DCC[j] ))
                              pm.Lads++; else */
        is_ss = false;
        PMM += pm.MM[j]; // /pm.L1[k]; bugfix DK 28.02.2012
        break;
      }
      // Inserted 08.02.2011 by DK to prevent crashes on inconsistent phases
      if (non_sorption_phase && is_ss) {
        char buf[128];
        sprintf(buf, ": %d %2.16s %c  ( DC: %d  %2.16s %c ) ", k, pm.SF[k],
                pm.PHC[k], j, pm.SM[j], pm.DCC[j]);
        vfMessage(window(), "Inconsistent DC in non-sorption phase", buf);
        pm.Lads--;
      }
      if (!non_sorption_phase && !is_ss) {
        char buf[128];
        sprintf(buf, ": %d %2.16s %c  ( DC: %d  %2.16s %c ) ", k, pm.SF[k],
                pm.PHC[k], j, pm.SM[j], pm.DCC[j]);
        vfMessage(window(), "Inconsistent DC in sorption phase", buf);
      }
    } /* j */
    if (PMM > 0.) // cannot calculate better before the composition of the phase
                  // is known
      PMM /= pm.L1[k]; // DK 28.02.2012

    if (pm.PHC[k] == PH_POLYEL || pm.PHC[k] == PH_SORPTION ||
        pm.PHC[k] == PH_ADSSBEIL) { // || pm.PHC[k] == PH_ADSORPT)
      /* Analyzing assignment of surface site types */
      ph_sur_site_type_assign(k, jb, je, car_l, car_c, Cjs);
      // ka++; kka++;
    }
    if (syp->PYOF == S_OFF)
      pm.YOF[k] = 0.0;
    else
      switch (pm.PHC[k]) { // convert free energy of metastable phase from J/g
                           // to J/mol */
      case PH_POLYEL:
      case PH_SORPTION:
      case PH_ADSSBEIL:
      case PH_ADSORPT:
        if (Cjs >= 0) {
          PMM = pm.MM[Cjs];
          pm.FWGT[k] = PMM;
        }
        [[fallthrough]];
        /*  leave normalized to 1 gram pm.ln1Na[k] *= PMM; */
      case PH_AQUEL:
        if (pm.FWGT[k] < 1e-33) // bugfix 26.04.13 DK
        {
          pm.YOF[k] = syp->YOF[kk] * PMM / pm.RT;
          pm.FWGT[k] = PMM;
        }
        break;
      default:
        if (Cjs >= 0) // if there is a MAJOR component in solution phase
          PMM = pm.MM[Cjs];
        if (pm.FWGT[k] < 1e-33) // bugfix 26.04.13 DK
        {
          pm.YOF[k] = syp->YOF[kk] * PMM / pm.RT;
          pm.FWGT[k] = PMM;
        }
        break;
      }

    // test of kinetic constraints (units of measurement!)
    if (pm.PLIM && k >= pm.FIs) { /* restrictions set! */
      /* one-component phases */
      if (syp->PLLim != S_OFF || syp->PULim != S_OFF) {
        pm.RLC[pm.Ls + k - pm.FIs] = syp->RFLC[kk];
        pm.RSC[pm.Ls + k - pm.FIs] = syp->RFSC[kk];
      }
      if (syp->PUL)
        pm.DUL[pm.Ls + k - pm.FIs] = syp->PUL[kk];
      if (syp->PLL)
        pm.DLL[pm.Ls + k - pm.FIs] = syp->PLL[kk];
    }
    if (pm.PLIM && k < pm.FIs)
    /*    ( syp->PLLim != S_OFF || syp->PULim != S_OFF ) )  */
    { /* solutions */
      if (syp->RFLC)
        pm.RFLC[k] = syp->RFLC[kk];
      else
        pm.RFLC[k] = NO_LIM;
      if (syp->RFSC)
        pm.RFSC[k] = syp->RFSC[kk];
      else
        pm.RFSC[k] = QUAN_MOL;
      if (syp->PUL)
        pm.PUL[k] = syp->PUL[kk];
      else
        pm.PUL[k] = 1e6;
      if (syp->PLL)
        pm.PLL[k] = syp->PLL[kk];
      else
        pm.PLL[k] = 0.0;
    }
    if (!pm.PLIM && k < pm.FIs) { /* no set restrictions */
      pm.RFLC[k] = NO_LIM;
      pm.RFSC[k] = QUAN_MOL;
      pm.PUL[k] = 1e6;
      pm.PLL[k] = 0.0;
    }

  } /* k */
}

// SCM Modernization: Loading surface type and site parameters
// Corrected for keeping sizes of sorption input arrays pm.FIa
// instead of pm.FIs
void TMulti::ph_sur_param(int k /*, int ka, int kk, int kka */) {
  int i;
  RMULTS *mup = TRMults::sm->GetMU();
  SYSTEM *syp = TSyst::sm->GetSY();

  int ka = k - (pm.FIs - pm.FIa);
  int kk = pm.muk[k];
  int kka = kk - (mup->Fis - mup->Fia);

  if (!pm.SCMC || pm.FIat <= 0 || k < 0 || k >= pm.FIs || ka < 0 ||
      ka >= pm.FIa || kka < 0 || kka >= mup->Fia)
    return; // some problem with dimensions

  // Tier 1: Surface types / patches [MNST=6]
  for (i = 0; i < MNST; i++) {
    if (syp->PNfsp != S_OFF && syp->SCMC) {
      pm.SCMC[ka][i][0] = syp->SCMC[kka][i][0];
      pm.SCMC[ka][i][1] = syp->SCMC[kka][i][1];
    } else {
      pm.SCMC[ka][i][0] = SC_NOT_USED;
      pm.SCMC[ka][i][1] = '\0';
    }

    if (pm.SCMC[ka][i][0] == SC_NOT_USED ||
        pm.SCMC[ka][i][0] == '\0') { /* this type surface is not used */
      pm.SCMC[ka][i][0] = SC_NOT_USED;
      pm.SCMC[ka][i][1] = '\0';
      for (int p = 0; p < NMSDT; p++)
        pm.MASDT[ka][i][p] = 0.0;
      pm.XcapA[ka][i] = 0.0;
      pm.XcapB[ka][i] = 0.0;
      pm.XcapC[ka][i] = 0.0;
      pm.XcapZ[ka][i] = 0.0;
      pm.XcapD[ka][i] = 0.0;
      pm.XdlA[ka][i] = 0.0;
      pm.XdlB[ka][i] = 0.0;
      pm.XdlC[ka][i] = 0.0;
      pm.XdlD[ka][i] = 0.0;
      pm.XlamA[ka][i] = 0.0;
      for (int p = 0; p < NCIEC; p++)
        pm.XfIEC[ka][i][p] = 0.0;
      continue;
    }

    // Tier 1 surface type/patch parameters: MASDT [NMSDT=4]
    for (int p = 0; p < NMSDT; p++) {
      if (syp->PMaSdt != S_OFF && syp->MASDT &&
          !IsFloatEmpty(syp->MASDT[kka][i][p]))
        pm.MASDT[ka][i][p] = syp->MASDT[kka][i][p];
      else
        pm.MASDT[ka][i][p] = 0.0;
    }

    if (syp->PXcapA != S_OFF && syp->XcapA && !IsFloatEmpty(syp->XcapA[kka][i]))
      pm.XcapA[ka][i] = syp->XcapA[kka][i];
    else
      pm.XcapA[ka][i] = 0.0;

    if (syp->PXcapB != S_OFF && syp->XcapB && !IsFloatEmpty(syp->XcapB[kka][i]))
      pm.XcapB[ka][i] = syp->XcapB[kka][i];
    else
      pm.XcapB[ka][i] = 0.0;

    // EIL C3 (4LM) capacitance density and slip plane parameter f_zeta
    if (syp->PXcapC != S_OFF && syp->XcapC && !IsFloatEmpty(syp->XcapC[kka][i]))
      pm.XcapC[ka][i] = syp->XcapC[kka][i];
    else
      pm.XcapC[ka][i] = 0.0;

    if (syp->PXcapC != S_OFF && syp->XcapZ && !IsFloatEmpty(syp->XcapZ[kka][i]))
      pm.XcapZ[ka][i] = syp->XcapZ[kka][i];
    else
      pm.XcapZ[ka][i] = 0.0;

    if (syp->PXfIEC != S_OFF &&
        syp->XfIEC) { /* Permanent surface charge density / IEC */
      for (int p = 0; p < NCIEC; p++) {
        if (!IsFloatEmpty(syp->XfIEC[kka][i][p]))
          pm.XfIEC[ka][i][p] = syp->XfIEC[kka][i][p];
        else
          pm.XfIEC[ka][i][p] = 0.0;
      }
    } else {
      for (int p = 0; p < NCIEC; p++)
        pm.XfIEC[ka][i][p] = 0.0;
    }

    if (syp->PXlam != S_OFF && syp->Xlam && !IsFloatEmpty(syp->Xlam[kka][i]))
      pm.XlamA[ka][i] = syp->Xlam[kka][i];
    else
      pm.XlamA[ka][i] = 0.0;

    pm.XdlA[ka][i] = 0.0;
    pm.XdlB[ka][i] = 0.0;
    pm.XdlC[ka][i] = 0.0;
    pm.XdlD[ka][i] = 0.0;
    pm.XcapD[ka][i] = 0.0;
  } /* end for i */

  // SCM Modernization: Tier 2: Surface Site Types [NSURFSITES=8][NFSIT=4]
  for (int s = 0; s < NSURFSITES; s++) {
    for (int p = 0; p < NFSIT; p++) {
      if (syp->PNfsp != S_OFF && syp->Nfsp &&
          !IsFloatEmpty(syp->Nfsp[kka][s][p]))
        pm.Nfsp[ka][s][p] = syp->Nfsp[kka][s][p];
      else
        pm.Nfsp[ka][s][p] = 0.0;
    }
  }
}

// SCM Modernization: Assigning surface site types to carrier and loading Tier 3
// species data
void TMulti::ph_sur_site_type_assign(int k, int jb, int je, int car_l[],
                                     int car_c, int Cjs) {
  int j, jcl;
  long int ist, isi;
  RMULTS *mup = TRMults::sm->GetMU();
  SYSTEM *syp = TSyst::sm->GetSY();

  if (!pm.SATC || !pm.SATT || !pm.MASDJ || pm.Lads <= 0)
    return;
  int ka = k - (pm.FIs - pm.FIa);
  int kk = pm.muk[k];

  for (j = jb; j < je; j++) {
    if (pm.DCC[j] == DC_PEL_CARRIER || pm.DCC[j] == DC_SUR_CARRIER ||
        pm.DCC[j] == DC_SUR_MINAL)
      continue;

    int ja = j - (pm.Ls - pm.Lads);
    int jj = pm.muj[j];
    int jja = jj - (mup->Ls - mup->Lads);

    if (ja < 0 || ja >= pm.Lads || jja < 0 || jja >= mup->Lads)
      continue;

    // Tier 3: SACT codes and allocations from SATC
    if (syp->PSATT != S_OFF && syp->SATC) {
      pm.SATT[ja] = syp->SATC[jja][SURSP_INDIV_SACT_CODE];

      // Allocation of surface species to surface type (patch)
      switch (syp->SATC[jja][SURSP_ALLOC_SURFTYPE]) {
      default:
      case CST_0:
        pm.SATC[ja][XL_ST] = 0;
        break;
      case CST_1:
        pm.SATC[ja][XL_ST] = 1;
        break;
      case CST_2:
        pm.SATC[ja][XL_ST] = 2;
        break;
      case CST_3:
        pm.SATC[ja][XL_ST] = 3;
        break;
      case CST_4:
        pm.SATC[ja][XL_ST] = 4;
        break;
      case CST_5:
        pm.SATC[ja][XL_ST] = 5;
        break;
      }

      // Allocation of surface species to sorbent (or its end-members)
      switch (syp->SATC[jja][SURSP_ALLOC_SORBENT_EM]) {
      case CCA_0:
        jcl = 0;
        break;
      case CCA_1:
        jcl = 1;
        break;
      case CCA_2:
        jcl = 2;
        break;
      case CCA_3:
        jcl = 3;
        break;
      case CCA_4:
        jcl = 4;
        break;
      case CCA_5:
        jcl = 5;
        break;
      case CCA_VOL:
      case 'v':
        jcl = -1; // whole sorbent
        break;
      default:
        jcl = -2; // main DC of the sorbent
        break;
      }

      pm.SATC[ja][XL_SP] = 0; // Obsolete: plane charges now handled via mCDc
    } else {
      pm.SATT[ja] = SAT_INDEF;
      pm.SATC[ja][XL_ST] = 0;
      jcl = -1;
      pm.SATC[ja][XL_SP] = 0;
    }

    // Tier 3: Loading allocation to surface sites from XSMD (species denticity
    // map)
    pm.SATC[ja][XL_SI] = 0;
    if (syp->XSMD) {
      for (int s = 0; s < NSURFSITE; s++) {
        if (syp->XSMD[jja][s] > 0) {
          pm.SATC[ja][XL_SI] = s;
          break;
        }
      }
      if (pm.XSMD) {
        for (int s = 0; s < NSURFSITE; s++) {
          pm.XSMD[ja][s] = syp->XSMD[jja][s] - '0';
        }
      }
    }

    // Tier 3: Loading MASDJ (parameters for SACT and Frumkin/isotherms)
    if (syp->PMaSdj != S_OFF && syp->MASDJ) {
      for (int p = 0; p < DNFCN; p++) {
        if (!IsFloatEmpty(syp->MASDJ[jja][p]))
          pm.MASDJ[ja][p] = syp->MASDJ[jja][p];
        else
          pm.MASDJ[ja][p] = 0.0;
      }
    } else {
      for (int p = 0; p < DNFCN; p++)
        pm.MASDJ[ja][p] = 0.0;
    }

    // Tier 3: Loading mCDc (multi-plane charge allocation)
    if (syp->mCDc && pm.mCDc) {
      for (int p = 0; p < DMCDC; p++) {
        if (!IsFloatEmpty(syp->mCDc[jja][p]))
          pm.mCDc[ja][p] = syp->mCDc[jja][p];
        else
          pm.mCDc[ja][p] = 0.0;
      }
    }

    // Assigning carrier end-member DC
    if (jcl < car_c && jcl >= 0)
      pm.SATC[ja][XL_EM] = car_l[jcl];
    else if (jcl >= car_c || jcl == -1)
      pm.SATC[ja][XL_EM] = -1; /* whole carrier */
    else
      pm.SATC[ja][XL_EM] = Cjs; /* main carrier DC */

    ist = pm.SATC[ja][XL_ST];
    isi = pm.SATC[ja][XL_SI];
    if (pm.D && ist >= 0 && ist < MNST && isi >= 0 && isi < MNST)
      pm.D[ist][isi] = 0.0; // cleaning the totals for sites

    // Convert surface density from 1/nm2 to mkmol/g using specific surface area
    // Aalp
    if (pm.SATT[ja] != SAT_COMP && syp->PAalp != S_OFF && syp->Aalp &&
        !IsFloatEmpty(syp->Aalp[kk]))
      pm.MASDJ[ja][PI_DEN] *= syp->Aalp[kk] * 1.66054; /* 1/nm2 to mkmol/g */

    if (approximatelyZero(pm.MASDJ[ja][PI_DEN])) {
      if (pm.SATT[ja] == SAT_NCOMP && ist >= 0 && ist < MNST && ka >= 0 &&
          ka < pm.FIa && pm.MASDT)
        pm.MASDJ[ja][PI_DEN] = pm.MASDT[ka][ist][0];
      // SCM Modernization: patch reactive area fraction
      else if (pm.SATT[ja] == SAT_COMP)
        pm.MASDJ[ja][PI_DEN] = 1.0;
    }
  } /* for j */
}

// End of file ms_mupmph.cpp
// --------------------------------------------------------------------
