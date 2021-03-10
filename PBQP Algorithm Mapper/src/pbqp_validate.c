/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Program for validating a PBQP solution. 
 *****************************************************************************
 $Id: pbqp_validate.c,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $ 
 *****************************************************************************
 $Log: pbqp_validate.c,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources


 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "vec_mat.h"
#include "bf_pbqp.h"


/*****************************************************************************
 * I/O operations 
 ****************************************************************************/

/* read PBQP problem */
pbqp *read_pbqp(FILE *f)
{
   int num_nodes;
   int num_edges;
   pbqp *this;
   int u,v,e;
   num min;

   ASSERT(f != NULL);
  
   ASSERT(fscanf(f,"%d %d",&num_nodes, &num_edges)==2);
   this = bf_alloc_pbqp(num_nodes);

   /* read in node costs */
   for(u=0;u<num_nodes;u++) {
      vec *v = v_read(f);
      bf_add_nodecosts(this,u,v); 
      v_free(v);
   }
   /* read in edges and its costs */
   for(e=0;e<num_edges;e++) {
     mat *m;
     ASSERT(fscanf(f,"%d %d",&u,&v) == 2);
     m  = m_read(f);
     ASSERT(u >= 0 && u < num_nodes);
     ASSERT(v >= 0 && v < num_nodes);
     bf_add_edgecosts(this,u,v,m);
     m_free(m);
   }
   /* read solutions of decision vectors */
   for(u=0;u<num_nodes;u++) {
     int sol;
     ASSERT(fscanf(f,"%d",&sol) == 1);
     bf_set_solution(this,u,sol);
   }

   /* set minimum */
   ASSERT(fscanf(f,"%f",&min) == 1);
   bf_set_min(this,min);

   return this;
}

/* test program */
int main(int argc,char **argv)
{
   FILE *f;

   assert (argc != 1);
   assert ((f = fopen(argv[1],"r")) != NULL);
   
   /* read PBQP problem */ 
   pbqp *this = read_pbqp(f);

   return !bf_validate_pbqp(this);
}
