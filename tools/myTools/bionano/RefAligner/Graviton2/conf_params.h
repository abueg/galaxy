#ifndef CONFPARAMS_H
#define CONFPARAMS_H

#include "ident.h"

static Ident conf_params_h_Id("$Header: http://svn.bnm.local:81/svn/Informatics/RefAligner/branches/11442/conf_params.h 8243 2018-12-08 21:59:07Z tanantharaman $");

//these are the maximum sizes for each axis (array)
#define SVCONFLEN 50
#define SVSIZELEN 50

extern int sv_conflen;
extern int sv_sizelen;
extern int sv_ppvbin_conf[SVCONFLEN];
extern int sv_ppvbin_sizekb[SVSIZELEN];
extern float sv_ppvbin_del[SVCONFLEN][SVSIZELEN];
extern float sv_ppvbin_ins[SVCONFLEN][SVSIZELEN];

#define INV_PAIRED_LEN 5
extern const double inv_paired_vec[INV_PAIRED_LEN];/* for 0, 1, 2, 3, 4+ label intervals */

#endif
