/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Header file for generic PBQP solver interface

 note: following definitions represent an abstract interface for various PBQP 
       approaches. 
 *****************************************************************************
 $Id: interface.h,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $
 *****************************************************************************
 $Log: interface.h,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/
#ifndef __INTERFACE_H__
#define __INTERFACE_H__

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

/**************
 * call table *
 **************/

struct pbqp_interface {
   /* add node costs */
   void (*add_nodecosts)(pbqp *th,int u, vec *costs);
   /* add edge mat */
   void (*add_edgecosts)(pbqp *th,int u,int v,mat *costs);
   /* solve PBQP problem */
   void (*solve_pbqp)(pbqp *th);
   /* get minimum of PBQP */
   num (*get_min)(pbqp *th);
   /* get solution of a node */
   int (*get_solution)(pbqp *th,int u);
   /* alloc PBQP */
   pbqp *(*alloc_pbqp)(int num);
   /* free PBQP */
   void (*free_pbqp)(pbqp *th);
   /* get number of nodes */
   int (*get_numnodes)(pbqp *th);
   /* is solution optimal */
   boolean (*is_optimal)(pbqp *th);
};


extern struct pbqp_interface _pbqp_calltable[];

/* Computation methods for PBQP 
 */

#define PBQP_BRUTEFORCE  (0)
#define PBQP_HEURISTICAL (1)
#define PBQP_BRANCHBOUND (2)

/* Generic Interface
 * =================
 *
 * First parameter of following method calls determines solving method.
 * A solving method is a entry in the calltable and can be either
 * brutefore, heuristical or branch&bound.
 */


/* add node costs */
#define ADD_NODECOSTS(m,th,u,costs) (*_pbqp_calltable[m].add_nodecosts)(th,u,costs)

/* add edge mat */
#define ADD_EDGECOSTS(m,th,u,v,costs) (*_pbqp_calltable[m].add_edgecosts)(th,u,v,costs)

/* solve PBQP problem */
#define SOLVE_PBQP(m,th) (*_pbqp_calltable[m].solve_pbqp)(th)

/* get minimum of PBQP */
#define GET_MIN(m,th) (*_pbqp_calltable[m].get_min)(th)

/* get solution of a node */
#define GET_SOLUTION(m,th,u) (*_pbqp_calltable[m].get_solution)(th,u)

/* alloc PBQP */
#define ALLOC_PBQP(m,num) (*_pbqp_calltable[m].alloc_pbqp)(num)

/* free PBQP */
#define FREE_PBQP(m,th) (*_pbqp_calltable[m].free_pbqp)(th)

/* get number of nodes */
#define GET_NUMNODES(m,th) (*_pbqp_calltable[m].get_numnodes)(th)

/* is solution optimal */
#define IS_OPTIMAL(m,th) (*_pbqp_calltable[m].is_optimal)(th)

#endif
