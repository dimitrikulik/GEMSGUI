//-------------------------------------------------------------------
// $Id: GraphDialog.h 968 2007-12-13 13:23:32Z gems $
//
// Declaration of GraphDialog, PlotTypeBtn classes (Plotting system)
//
// Copyright (C) 1996-2008  A.Rysin, S.Dmytriyeva
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
#ifndef graph_dialog_h
#define graph_dialog_h

#include <qcolor.h>
#include <QResizeEvent>
#include <QDialog>

#include "ui_GraphDialog4.h"

#include  "array.h"
#include  "v_module.h"
#include  "pshape.h"
#include  "graph.h"

// Widgets for legend and plotting dialogs

class PlotTypeBtn:  public QPushButton
{
    TPlotLine plLine;

protected:
    void drawButtonLabel(QPainter* paint);
    void paintEvent( QPaintEvent * event );
    
public:
    PlotTypeBtn(TPlotLine& ln, QWidget* parent=0, const char* name=0):
            QPushButton( parent ),
	    plLine(ln)
    { }


    void setName( const char* );
    void setData( TPlotLine& ln);
    TPlotLine& getData()
    {
        return plLine;
    }

    void setColor( QColor& );
    QColor getColor() const
    {
        return QColor( plLine.red, plLine.green, plLine.blue);
    }

    static void drawSymbol(QPainter* painter, const QPoint& center,
         int type, int size, const QColor& color, int width=1);
};

class SymbolLabel;
class DragLabel;

class GraphDialog: public QDialog, public Ui::GraphDialogData
{
    Q_OBJECT

    TPlotWin* plot;
    TCModule *pModule;

    // work part
    bool isFragment;
    float minX, minY;
    float maxX, maxY;
//    TOArray<SymbolLabel*> aSymbolLabels;
//    TOArray<DragLabel*> aLegendLabels;
    QGridLayout *pGrp;
 //   TIArray<SymbolLabel> aSymbolLabels;
 //   TIArray<DragLabel> aLegendLabels;

    void Show();
    void ShowPlots();

    bool isoline_put;
    bool lines_put;

    void ShowLegend();
    void ShowIsolineLegend();


protected slots:
    virtual void languageChange();

    virtual void CmFragment();
    virtual void CmLegend();
    virtual void CmPrint();
    virtual void CmSave();

protected:
    void resizeEvent(QResizeEvent* qpev);

public:

    GraphData gr_data;

    GraphDialog(TCModule *pmodule, GraphData& data);
    ~GraphDialog();


    void ShowNew();
    void Apply();    // Update changes
    void AddPoint( int nPlot, int nPoint, bool no_mt ); // Add new point to graph

    QColor getColor(int ii) const
    {
        return QColor( gr_data.lines[ii].red,
        	    gr_data.lines[ii].green, gr_data.lines[ii].blue);
    }

    QColor getColorIsoline(int ii) const
    {
        return QColor( gr_data.scale[ii].red,
        	    gr_data.scale[ii].green, gr_data.scale[ii].blue);
    }

    gstring getTextIsoline(int ii);

    QColor getBackgrColor();
    void setBackgrColor( QColor color );
};


#endif   // graph_dialog_h