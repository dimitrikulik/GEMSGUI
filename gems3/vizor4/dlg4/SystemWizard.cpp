//-------------------------------------------------------------------
// $Id: SystemWizard.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of SystemWizard class
//
// Copyright (C) 2010  S.Dmytriyeva, D.Kulik
// Uses  gstring class (C) A.Rysin 1999
//
// This file is part of the GEM-Vizor library which uses the
// Qt v.2.x GUI Toolkit (Troll Tech AS, http://www.trolltech.com)
// according to the Qt Duo Commercial license
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------

const char *SYST_HTML = "gems_sys";

#include <qcheckbox.h>

#include "SystemWizard.h"
#include "visor_w.h"

void SystemWizard::languageChange()
{
    retranslateUi(this);
}

void SystemWizard::CmBack()
{
    stackedWidget->setCurrentIndex ( stackedWidget->currentIndex()-1 );
    resetNextButton();
    resetBackButton();
}

void SystemWizard::CmNext()
{
	stackedWidget->setCurrentIndex ( stackedWidget->currentIndex()+1 );
    resetNextButton();
    resetBackButton();
}

void 	SystemWizard::resetNextButton()
{
	if( stackedWidget->currentIndex() == stackedWidget->count() - 1 )
	{	
		pNext->disconnect();
		connect( pNext, SIGNAL(clicked()), this, SLOT(accept()) );
		pNext->setText("&Finish");
	}
	else
	{	
		pNext->disconnect();
		connect( pNext, SIGNAL(clicked()), this, SLOT(CmNext()) );
		pNext->setText("&Next>");
	}
}

void 	SystemWizard::resetBackButton()
{
	pBack->setEnabled( stackedWidget->currentIndex() > 0 );
}

SystemWizard::SystemWizard( const char* pkey, char flgs[40],
    gstring name, gstring comment, gstring EQkey,  QWidget* parent):
    QDialog( parent )
{

    //setFinishEnabled( WizardPage3, true);
    setupUi(this);
    gstring str1= "GEM-Selektor System Setup:  ";
            str1 += pkey;
            setWindowTitle( str1.c_str() );

    QObject::connect( pHelp, SIGNAL(clicked()), this, SLOT(help()));
    QObject::connect( pBack, SIGNAL(clicked()), this, SLOT(CmBack()));
    QObject::connect( pNext, SIGNAL(clicked()), this, SLOT(CmNext()));
    stackedWidget->setCurrentIndex (0);
    resetNextButton();
    resetBackButton();

// page1
    pName->setText(name.c_str());
    pComment->setText(comment.c_str());

    if( flgs[13] != '-' )
            pPbPh->setChecked( true );
    else pPbPh->setChecked( false );
    if( flgs[15] != '-' || flgs[17] != '-')
            pPLim->setChecked( true );
    else pPLim->setChecked( false );

    if( flgs[0] != '-' )
            pPE->setChecked( true );
    else pPE->setChecked( false );

    if( flgs[1] != '-' )
            pPV->setChecked( true );
    else pPV->setChecked( false );

    if( flgs[9] != '-' )
            pPBeq->setChecked( true );
    else pPBeq->setChecked( false );

    if( flgs[24] != '-' )
            pPAalp->setChecked( true );
    else pPAalp->setChecked( false );

    if( flgs[25] != '-' )
            pPSigm->setChecked( true );
    else pPSigm->setChecked( false );

    if( flgs[26] != '-' )
            pXr0h0->setChecked( true );
    else pXr0h0->setChecked( false );

    if( flgs[27] != '-' )
            pPXepsC->setChecked( true );
    else pPXepsC->setChecked( false );

    if( flgs[28] != '-' )
            pPNfsp->setChecked( true );
    else pPNfsp->setChecked( false );

    if( flgs[34] != '-' )
            pPParc->setChecked( true );
    else pPParc->setChecked( false );

    if( flgs[35] != '-' )
            pPDelB->setChecked( true );
    else pPDelB->setChecked( false );

    if( flgs[36] != '-' )
            pPXlam->setChecked( true );
    else pPXlam->setChecked( false );


// Page 2
    pKey->setText(EQkey.c_str());

}

SystemWizard::~SystemWizard()
{}

gstring SystemWizard::getName(  )
{
  return gstring(pName->text().toLatin1().data());
}

gstring SystemWizard::getComment(  )
{
  return gstring(pComment->text().toLatin1().data());
}

gstring SystemWizard::getEQkey(  )
{
  return gstring(pKey->text().toLatin1().data());
}


void SystemWizard::getFlags( char flgs[40] )
{
// Page 1 
    if( pPbPh->isChecked() )
          flgs[13] = '+';
    else flgs[13] = '-';

    if( pPLim->isChecked() )
          flgs[15] = flgs[17] = '+';
    else flgs[15] = flgs[17] = '-';

    if( pPE->isChecked() )
          flgs[0] = '+';
    else flgs[0] = '-';

    if( pPV->isChecked() )
          flgs[1] = '+';
    else flgs[1] = '-';

    if( pPBeq->isChecked() )
          flgs[9] = '+';
    else flgs[9] = '-';

    if( pPAalp->isChecked() )
          flgs[24] = '+';
    else flgs[24] = '-';

    if( pPSigm->isChecked() )
          flgs[25] = '+';
    else flgs[25] = '-';

    if( pXr0h0->isChecked() )
          flgs[26] = '+';
    else flgs[26] = '-';

    if( pPXepsC->isChecked() )
          flgs[27] = '+';
    else flgs[27] = '-';

    if( pPNfsp->isChecked() )
          flgs[28] = '+';
    else flgs[28] = '-';

    if( pPParc->isChecked() )
          flgs[34] = '+';
    else flgs[34] = '-';

    if( pPDelB->isChecked() )
          flgs[35] = '+';
    else flgs[35] = '-';

    if( pPXlam->isChecked() )
          flgs[36] = '+';
    else flgs[36] = '-';

}

void
SystemWizard::help()
{
  pVisorImp->OpenHelp( SYST_HTML, 0, this/*, true*/ );
}

//--------------------- End of SystemWizard.cpp ---------------------------