//-------------------------------------------------------------------
// $Id: SelectDialog.h 968 2007-12-13 13:23:32Z gems $
//
// Declaration of SelectDialog class
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

#ifndef SelectDialog_included
#define SelectDialog_included

#include <QDialog>
#include <QListWidget>

#include "ui_SelectDialog4.h"
#include "gstring.h"
#include "v_user.h"

class SelectDialog : public QDialog, public Ui::SelectDialogData
{
    Q_OBJECT

    bool multi;
    bool sel_to_all;

protected slots:
    virtual void languageChange();

    virtual void CmSelectAll();
    virtual void CmClearAll();

public:

    /*! single select constructor
    */

    SelectDialog(QWidget* parent, const char* title,
                 TCStringArray& list, int sel);

    SelectDialog(QWidget* parent, const char* title,
                 TCStringArray& list, int sel, bool all_ );

    /*! multiple select constructor
    */
    SelectDialog(QWidget* parent, const char* title,
                 TCStringArray& list, TCIntArray& sel);

    virtual ~SelectDialog();

    /*! returns single selection
        returns '-1' if nothing selected
    */
    int selected()
    {
        return (( result() ) ? pList->currentRow() : -1);
    }

    int selected( bool& all_ )
    {
        all_ = sel_to_all;
        return ( ( result() ) ? pList->currentRow() : -1);
    }

    TCIntArray allSelected();

};

#endif // SelectDialog_included