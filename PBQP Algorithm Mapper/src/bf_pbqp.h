/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Header file for brute force PBQP Solver
 *****************************************************************************
 $Id: bf_pbqp.h,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $
 *****************************************************************************
 $Log: bf_pbqp.h,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/
#ifndef __BF_PBQP_H__
#define __BF_PBQP_H__

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
void bf_add_nodecosts(pbqp *this,int u, vec *costs);

/* add edge mat */
void bf_add_edgecosts(pbqp *this,int u,int v,mat *costs);

/* solve PBQP problem */
void bf_solve_pbqp(pbqp *this);

/* get minimum of PBQP */
num bf_get_min(pbqp *this);

/* get solution of a node */
int bf_get_solution(pbqp *this,int u);

/* alloc PBQP */
pbqp *bf_alloc_pbqp(int num);

/* free PBQP */
void bf_free_pbqp(pbqp *this);

/* get number of nodes */
int bf_get_numnodes(pbqp *this);

/* is optimal */
boolean bf_is_optimal(pbqp *this);

/* set solution of a node */
void bf_set_solution(pbqp *this,int u,int sol);

/* set min */
void bf_set_min(pbqp *this,num min);

/* validate solution */
boolean bf_validate_pbqp(pbqp *this);

#endif
