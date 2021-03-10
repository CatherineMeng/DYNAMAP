/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Generic PBQP solver interface
 *****************************************************************************
 $Id: interface.c,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $
 *****************************************************************************
 $Log: interface.c,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "vec_mat.h"
#include "bf_pbqp.h"
#include "h_pbqp.h"

#include "interface.h"

/**************
 * call table *
 **************/

struct pbqp_interface _pbqp_calltable[] = 
{ { /* Brute Force calls */ 
   bf_add_nodecosts,
   bf_add_edgecosts,
   bf_solve_pbqp,
   bf_get_min,
   bf_get_solution,
   bf_alloc_pbqp,
   bf_free_pbqp,
   bf_get_numnodes,
   bf_is_optimal}
  ,
  { /* Heuristical calls */ 
   h_add_nodecosts,
   h_add_edgecosts,
   h_solve_pbqp,
   h_get_min,
   h_get_solution,
   h_alloc_pbqp,
   h_free_pbqp,
   h_get_numnodes,
   h_is_optimal},
};

/* end of interface.c */
