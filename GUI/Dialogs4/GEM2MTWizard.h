//-------------------------------------------------------------------
// $Id: GEM2MTWizard.h 487 2005-03-04 15:44:28Z sveta $
//
// Declaration of GEM2MTWizard class
//
// Copyright (C) 2005-2007  S.Dmytriyeva
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

#ifndef GEM2MTWizard_included
#define GEM2MTWizard_included

#include "EquatSetupWidget.h"
#include "model_w.h"

namespace Ui {
class GEM2MTWizardData;
}

enum gem2mt_remake_pages {
    mode_RMT = 0,
    size_RMT,
    fluxes_transport,
    auto_script,
    graphic_script,
    vtk_format,
    gems3k_exchange,
    T_P_lookup,
    other_options
};

class GEM2MTWizard : public QDialog
{
    Q_OBJECT

    Ui::GEM2MTWizardData *ui;
    string calcScript;
    string outScript;
    EquatSetup *pageScript;

    QButtonGroup *GroupBox1;

    // for VTK window
    int cPage;
    std::vector<pagesSetupData> stData;
    std::vector<pagesSetupData> pgData;
    std::vector<scriptSetupData>  scriptData;
    QList<QListWidget *> pLists;
    TCIntArray pNdx;


    void resetNextButton();
    void resetBackButton();
    void resetVTKList();
    inline int findVTKarr( int vtk1  );

    TObjectTable* PTable;
    TObjectTable* TTable;

    void setupPTArrays();
    void defineTArray();
    void definePArray();
    void initPTable();
    void initTTable();
    void showPTTable();

    // internal functions
    int  tableFindRow( int nO, int ndx );
    void tableInsertRow( int nO, int ndx, const char * andName );
    void tableDeleteRow( int row );
    void listUpdate( );
    //void emptyScriptTable();

public:

    GEM2MTWizard( const char* pkey, char flgs[32], int sizes[20],
    double Tai[4], double Pai[4], double Tau[3],
    const char *acalcScript, const char *aoutScript,
    const char* aXname, const char* aYname,
    TCIntArray vtk1, TCIntArray vtk2, QWidget* parent = nullptr);

    virtual ~GEM2MTWizard();

    void   getSizes( int size[20] );
    void   getFlags( char flgs[32] );
    void   getTdata( double Tai[4] );
    void   getPdata( double Pai[4] );
    void   getTaudata( double Tau[3] );

    string getCalcScript() const;
    string getOutScript() const
    { return pageScript->getScript(); }

    TCStringArray getNames( string& xName, string& yName ) const
    { return pageScript->getNames(xName, yName); }

    void getVTK( TCIntArray& vtk1, TCIntArray& vtk2  );
    void setVTK( TCIntArray vtk1, TCIntArray vtk2  );

protected slots:

    void arrayChecked( bool check );
    void UnCheckSolid( bool check );
    void UnCheckAqGas( bool check );

    void help();
    void CmNext();
    void CmBack();
    void resetPageList(const char* aXname, const char* aYname);
    void setTScript( bool c_PvMSt_checked );
    void ScriptChange(QAbstractButton*);
    void changePage( int nPage );
    void changeTable(const QItemSelection & selected, const QItemSelection & deselected);
    void disableVTK();
    void objectChanged();

};

#endif // GEM2MTWizard_included
