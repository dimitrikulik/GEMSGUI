//-------------------------------------------------------------------
// $Id: RTparmWizard.h 968 2007-12-13 13:23:32Z gems $
//
// Declaration of RTparmWizard class
//
// Copyright (C) 2005-2008  S.Dmytriyeva, D.Kulik
//
// This file is part of the GEM-Selektor GUI library which uses the
// Qt v.4 cross-platform App & UI framework (https://qt.io/download-open-source)
// under LGPL v.2.1 (http://www.gnu.org/licenses/lgpl-2.1.html)
//
// This file may be distributed under the GPL v.3 license

//
// See http://gems.web.psi.ch/ for more information
// E-mail gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef RTparmWizard_included
#define RTparmWizard_included

#include <QDialog>
#include "EquatSetupWidget.h"

namespace Ui {
class RTparmWizardData;
}

class RTparmWizard : public QDialog
{
    Q_OBJECT

    Ui::RTparmWizardData *ui;
    bool isDC;
    string calcScript;
    EquatSetup *pageScript;

    void 	resetNextButton();
    void 	resetBackButton();
    void        definePTArray();
    void        initPTTable();

public:

    RTparmWizard( const char* pkey, char flgs[10], int sizes[7], double val[6],
    const char *acalcScript, const char* aXname, const char* aYname,
    QWidget* parent = nullptr);
    virtual ~RTparmWizard();

    void   getSizes( int size[7] );
    void   getFlags( char flgs[10], string& xName );
    void   getFloat( double val[6] );

    string getScript() const
    { return pageScript->getScript(); }
    TCStringArray getNames( string& xName, string& yName ) const
    { return pageScript->getNames(xName, yName); }

protected slots:

    void help();
    void CmNext();
    void CmBack();
    void resetPageList(const char* aXname, const char* aYname);
    void PChange();
    void TChange();
    void PairwiseChecked();
    void objectChanged(){}
};


#endif // RTparmWizard_included
