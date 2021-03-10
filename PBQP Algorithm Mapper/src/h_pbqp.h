/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Header file for heuristical PBQP Solver
 *****************************************************************************
 $Id: h_pbqp.h,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $
 *****************************************************************************
 $Log: h_pbqp.h,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/
#ifndef __H_PBQP_H__
#define __H_PBQP_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "vec_mat.h"

#ifndef PBQP_TYPE
#define PBQP_TYPE
struct pbqp;
typedef struct pbqp pbqp;
#endif

/*****************
 * PBQP routines *
 *****************/

/* add node costs */
void h_add_nodecosts(pbqp *this,int u, vec *costs);

/* add edge mat */
void h_add_edgecosts(pbqp *this,int u,int v,mat *costs);

/* solve PBQP problem */
void h_solve_pbqp(pbqp *this);

/* get minimum of PBQP */
num h_get_min(pbqp *this);

/* get solution of a node */
int h_get_solution(pbqp *this,int u);

/* alloc PBQP */
pbqp *h_alloc_pbqp(int num);

/* free PBQP */
void h_free_pbqp(pbqp *this);

/* get number of nodes */
int h_get_numnodes(pbqp *this);

/* is optimal */
boolean h_is_optimal(pbqp *this);

/* set dump file */
void h_set_dumpfile(pbqp *this,FILE *file);

#endif
