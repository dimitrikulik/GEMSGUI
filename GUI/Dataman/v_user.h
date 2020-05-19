//-------------------------------------------------------------------
// $Id$
/// \file v_user.h
/// Declaration of platform-specific utility functions and classes
//
// Copyright (C) 1996,2001,2012 A.Rysin, S.Dmytriyeva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#ifndef _v_user_h_
#define _v_user_h_

#include <algorithm>
#include <iostream>
using namespace std;
#include "GEMS3K/v_detail.h"
#include "array.h"
#include "gstring.h"

#ifdef __APPLE__

#ifndef __unix
#define __unix
#endif
#ifndef __FreeBSD
#define __FreeBSD
#endif

typedef unsigned int uint;
#endif

const int MAXKEYWD = 6+1;

#ifndef  __unix

typedef unsigned int uint;

#endif //  __noborl


void Gcvt(double number, size_t ndigit, char *buf);
double NormDoubleRound(double aVal, int digits);
void NormDoubleRound(double *aArr, int size, int digits);
void NormFloatRound(float *aArr, int size, int digits);

inline
int ROUND(double x )
{
    return int((x)+.5);
}


inline bool IsSpace(char ch)
{
    return ( (ch == ' ') || (ch == '\t') );
}

void StripLine(gstring& line);
void KeyToName(std::string& line);
TCStringArray split(const gstring& str, const gstring& delimiters);

// Added by SD on 22/12/2001
// Change string on templates
void
ChangeforTempl( gstring& data_str,  const gstring& from_templ1,
                const gstring& to_templ1, uint len_ );

// Returns string representation of current date in dd/mm/yyyy format
gstring curDate();

// Returns string representation of current date in dd/mm/yy format
std::string curDateSmol(char ch = '/');

// Returns string representation of current time in HH:MM  format
gstring curTime();

// Returns string representation of current date and time
inline
gstring curDateTime()
{
    return curDate() + curTime();
}

// reads line to gstring class from istream with a delimiter
istream& u_getline(istream& instream, gstring& dst_string, char delimit = '\n');
istream& f_getline(istream& is, std::string& str, char delim);


#ifdef __FreeBSD
// replacement for missing function in FreeBSD
inline char* gcvt(double num, int digit, char* buf)
{
    sprintf(buf, "%*g", digit, num);
    return buf;
}

#endif  // __FreeBSD


// dynamically allocates temporary 'char*'
// for simple string manipulations
// (used instead of stack char[] allocation to avoid stack problems)
struct vstr
{
    char* p;
    vstr(int ln): p(new char[ln+1])
    { }

    vstr(int ln, const char* s): p(new char[ln+1])    {
        strncpy(p, s, ln);
        p[ln]='\0';
    }

    vstr(const char* s): p(new char[strlen(s)+1])    {
       strcpy(p, s);
    }

    ~vstr()    {
        delete[] p;
    }

    operator char* ()    {
        return p;
    }

private:
    vstr (const vstr&);
    const vstr& operator= (const vstr&);

};


#define fileNameLength 64
/// Get Path of file and Reading list of file names from it, return number of files
char  (* f_getfiles(const char *f_name, char *Path, 
		long int& nElem, char delim ))[fileNameLength];



#endif // _v_user_h_

