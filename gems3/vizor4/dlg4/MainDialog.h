//-------------------------------------------------------------------
// $Id: MainDialog.h 968 2007-12-13 13:23:32Z gems $
//
// Declaration of MainDialog class
//
// Copyright (C) 1996-2007  A.Rysin, S.Dmytriyeva
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

#ifndef MainDialog_included
#define MainDialog_included

#include <QWidget>

#include "ui_MainDialog4.h"

class MainDialog : public QWidget, public Ui::MainDialogData
{
    Q_OBJECT

public:

    MainDialog( QWidget* parent = NULL );

    virtual ~MainDialog();

protected slots:
    virtual void languageChange();

    virtual void CmSettings();
    virtual void CmEQDialog();
    virtual void CmDBDialog();
    virtual void CmAbout();
    virtual void CmHelp();
    virtual void CmExit();

};

#endif // MainDialog_included