#ifndef IDENT_H
#define IDENT_H


#include <string.h>
#include <stdio.h>

extern const char *global_id;
/** SVN ID string class : 
  The constructor is designed to prevent gcc from removing the ID string from compiled binary
    due to its not being used anywhere in the program.

The intended use is to have each source code have a line (after #include "ident.h") as follows:

static Ident MyFile_id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/ident.h 1925 2013-09-12 19:45:59Z MRequa $);

The variable name MyFile_id does not matter, but should be unique. To work, SVN's propset command must be used once on any file which uses Ident. For this file the following command was used:

svn propset svn:keywords "Header" ident.h

Subsequent "svn commit" followed by "svn update" will then replace the $Header:...$ string with the latest svn database information.

see example usage right after the end of this class definition.

 */
class Ident {
 public:
  Ident(const char *ident){
    if(global_id){
      register int len = strlen(global_id) + strlen(ident);
      char *newstring = new char[len+2];
      sprintf(newstring,"%s %s",global_id,ident);
      delete [] global_id;
      global_id = newstring;
    } else {
      register int len = strlen(ident);
      char *newstring = new char[len+1];
      sprintf(newstring,"%s",ident);
      global_id = newstring;
    }
  }
  ~Ident() {
    if(global_id){
      delete [] global_id;
      global_id = 0;
    }
  }
};

static Ident ident_h_id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/ident.h 1925 2013-09-12 19:45:59Z MRequa $");

#endif
