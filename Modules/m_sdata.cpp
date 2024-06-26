//-------------------------------------------------------------------
// $Id: m_sdata.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of TSData class, config and calculation functions
//
// Rewritten from C to C++ by S.Dmytriyeva
// Copyright (C) 1995-2001 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
// Uses: GEM-Selektor GUI GUI DBMS library, gems/lib/gemvizor.lib
//
// This file may be distributed under the GPL v.3 license

//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include "m_sdata.h"
#include "service.h"
#include "visor.h"
//#include "filters_data.h"

TSData* TSData::pm;

TSData::TSData( uint nrt ):
        TCModule( nrt )
{
    refs = 0;
    title = abstr = 0;
    aFldKeysHelp.push_back("Script name or author's name for data source");
    aFldKeysHelp.push_back("Script number or year of publication { 1990 }");
    aFldKeysHelp.push_back("Script module name { DComp ReacDC } or data source type");
    //startKeyEdit = 0;
    set_def();
    start_title = " Scripts or Bibliographic References ";
    aObj[o_ssname]->set_empty_string("Please, enter here a title explaining what this chemical system is");
    aObj[o_ssnotes]->set_empty_string("Please, enter here a comment about the purpose of this system definition");
}

// link values to objects
void TSData::ods_link( int )
{
    aObj[o_sdnref]->SetPtr( &nREf );
    aObj[o_sdauth]->SetPtr( auth );
    aObj[o_sdyear]->SetPtr( year );
    aObj[o_sdclass]->SetPtr( type );
    aObj[o_sdvoly]->SetPtr( volyr );
    aObj[o_sdpage]->SetPtr( pages );
    aObj[o_sdauthr]->SetPtr( authors );
    aObj[o_sdedit]->SetPtr( editn );
    aObj[o_sdnote]->SetPtr( notes );
}

// set default data
void TSData::set_def(int)
{
    nREf = 0;
    dyn_new();
    strcpy( auth, "Einstein_ea" );
    strcpy( year, "2003" );
    strcpy( type, "wow" );
    strcpy( authors, "Einstein, A. and others" );
    strcpy( title, "Forgotten discovery of greatest importance" );
    strcpy( volyr, "2003, v.33" );
    strcpy( pages, "222-333" );
    strcpy( editn, "Journal of Systematic Research" );
    strcpy( notes, "Something wonderful" );
    strcpy( abstr, "Precise description of Universe" );
}

// reallocate dynamic memory
void TSData::dyn_new(int)
{
    if( title==0 )
        title = (char *)aObj[o_sdtitle]->Alloc( 1, 128, S_ );
    if( abstr==0 )
        abstr = (char *)aObj[o_sdabstr]->Alloc( 1, 128, S_ );
    if( nREf != 0 )
        refs = (char *)aObj[o_sdrefs]->Alloc( nREf, 1, 32 );
    else  if( refs )
        refs = (char *)aObj[o_sdrefs]->Free();
}

// set dynamic objects ptr to values
void TSData::dyn_set(int)
{
    title = (char *)aObj[o_sdtitle]->GetPtr();
    abstr = (char *)aObj[o_sdabstr]->GetPtr();
    refs = (char *)aObj[o_sdrefs]->GetPtr();
    nREf = (short)aObj[o_sdrefs]->GetN();
}

// free dynamic memory in objects and values
void TSData::dyn_kill(int)
{
    title = (char *)aObj[o_sdtitle]->Free();
    abstr = (char *)aObj[o_sdabstr]->Free();
    refs = (char *)aObj[o_sdrefs]->Free();
}

//Rebuild record structure before calc
int TSData::RecBuild( const char *key, int mode  )
{
AGAIN:
    int ret = TCModule::RecBuild( key, mode );
    if( ret == VF_CANCEL )
        return ret;
    if( nREf < 0 || nREf > 7 )
        goto AGAIN;
    dyn_new();
    pVisor->Update();
    return ret;
}

//Save record structure - added 09.04.03 KD
void TSData::RecSave( const char *key, bool onOld )
{
    TCModule::RecSave( key, onOld );
}

const char* TSData::GetHtml()
{
   return GM_SDREF_HTML;
}

// virtual callback for printing the record
// should be redefined in subclasses (modules)
// to make some plotting

void TSData::RecordPrint( const char* /*key*/ )
{
   std::string text_fmt = "line %s \"SDref record: \", %s rkey, %11s date, %6s time\n"
                      "line %s \"\"\n"
                      "line %s \"Authors: \", %s #SDauth\n"
                      "line %s \"Title: \", %s #SDtitl\n"
                      "line %s \"Edition: \", %s #SDedit\n"
                      "line %s \"Year: \", %s #SDvoly, %s \"Pages: \", %s #SDpage \n"
                      "line %s \"Comment\"\n"
                      "line %s #SDnote\n"
                      "line %s \"Abstract\"\n"
                      "line %s #SDabst\n"
                      "line %s \"------\"\n";
    PrintSDref( "pscript-script:0000:sdref:", text_fmt.c_str() );
}


void TSData::CopyRecords( const char *prfName, std::set<std::string>& SDlist )
{
    int Rnum;
    uint ii;

    // open selected kernel files
    int fnum_ = db->GetOpenFileNum( prfName );

    // add scripts for list to copy
    TCIntArray anR;
    TCStringArray aScripts;
    db->GetKeyList( "?script*:*:*:", aScripts, anR );
    for(ii=0; ii<aScripts.size(); ii++ )
       SDlist.insert(aScripts[ii]);

    // rename records from list
    for( const auto& sdkey: SDlist )
    {
      Rnum = db->Find( sdkey.c_str() );
      if( Rnum < 0 )
        continue;
      RecInput( sdkey.c_str() );

      // changing record key
      std::string str= std::string(db->FldKey( 2 ), 0, db->FldLen( 2 ));

      // we can change/add only one last (7) symbol
      // first (6) critical for scripts
      ChangeforTempl( str, "*", "*_", db->FldLen( 2 ));
      str += ":";
      std::string str1 = std::string(db->FldKey( 1 ), 0, db->FldLen( 1 ));
      strip( str1 );
      str = str1 + ":" + str;
      str1 = std::string(db->FldKey( 0 ), 0, db->FldLen( 0 ));
      strip( str1 );
      str = str1 + ":" + str;

      //Point SaveRecord
      Rnum = db->Find( str.c_str() );
      if(Rnum < 0 )
        AddRecordTest( str.c_str(), fnum_ );

    }

    // close all no project files
    TCStringArray names1;
    names1.push_back(prfName);
    db->OpenOnlyFromList(names1);
}

// Test record with key
void TSData::TryRecInp( const char *key_, time_t& time_s, int q )
{
    TCStringArray aSDlist;
    TCIntArray anRDc;

    if( ! MessageToSave() )
        return;

    TDBKey dbKey(db->GetDBKey());
    dbKey.SetKey(key_);
    std::string fld2 = std::string(dbKey.FldKey(2), 0, dbKey.FldLen(2));
    strip( fld2 );
    fld2 += "*";
    gui_logger->debug("TSData {} ", fld2);
    dbKey.SetFldKey(2,fld2.c_str());
    std::string str_key( dbKey.UnpackKey(), 0, db->KeyLen());
    RecStatus iRet = db->Rtest( str_key.c_str(), 1 );
    std::string msg;

    switch( iRet )
    {
    case MANY_:    // Get Key list
        db->GetKeyList(  str_key.c_str(), aSDlist, anRDc );
        db->Get(anRDc[0]);
        [[fallthrough]];
    case ONEF_:
        dyn_set(q);
        time_s = db->Rtime();
        return;
    case UNDF_:
    case NONE_:
    case EMPC_:
        {
            msg = "E16DCrun: In database chain  ";
            msg +=  GetName();
            msg += ": Data record not found, \n"
                   " key  '";
            msg += char_array_to_string( key_, db->KeyLen() );
            msg += "'.\n Maybe, a database file is not linked.\n";
            Error( GetName(), msg );
        }
        break;
    case FAIL_:
        msg = "E22DCrun: Failure!!! Database chain ";
        msg += GetName();
        msg += " is corrupt,\n"
               "Data record key '";
        msg += char_array_to_string( key_, db->KeyLen() );
        msg += "'\n Try to backup/restore or compress files in this database chain!";
        Error( GetName(),  msg );
    }
}

TSData* pmSData;

// --------------------- End of m_sdata.cpp --------------------------

