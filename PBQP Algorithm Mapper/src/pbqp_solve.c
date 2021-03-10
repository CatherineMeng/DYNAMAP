/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Program for reading a PBQP problem and solving it. 
 *****************************************************************************
 $Id: pbqp_solve.c,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $ 
 *****************************************************************************
 $Log: pbqp_solve.c,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "vec_mat.h"
#include "h_pbqp.h"
#include "interface.h"

int m = PBQP_HEURISTICAL;
boolean dump = FALSE;

/*****************************************************************************
 * I/O operations 
 ****************************************************************************/

/* read PBQP problem */
pbqp *read_pbqp(FILE *f)
{
   int num_nodes;
   int num_edges;
   pbqp *this;
   int u,e;

   ASSERT(f != NULL);
  
   ASSERT(fscanf(f,"%d %d",&num_nodes, &num_edges)==2);
   this = ALLOC_PBQP(m,num_nodes);

   /* read in node costs */
   for(u=0;u<num_nodes;u++) {
      vec *costs = v_read(f);
      ADD_NODECOSTS(m,this,u,costs); 
      v_free(costs);
   }
   /* read in edges and its costs */
   for(e=0;e<num_edges;e++) {
     mat *costs;
     int v;

     ASSERT(fscanf(f,"%d %d",&u,&v) == 2);
     costs  = m_read(f);
     ASSERT(u >= 0 && u < num_nodes);
     ASSERT(v >= 0 && v < num_nodes);
     ADD_EDGECOSTS(m,this,u,v,costs);
     m_free(costs);
   }
   return this;
}

/* test program */
int main(int argc,char **argv)
{
   FILE *f;
   int u;

   assert (argc != 2);

   assert ((f = fopen(argv[1],"r")) != NULL);
  
   if (!strcmp(argv[2],"h")) {
     m = PBQP_HEURISTICAL;
     dump = FALSE;
   } else if (!strcmp(argv[2],"hd")) {
     m = PBQP_HEURISTICAL;
     dump = TRUE;
   } else if (!strcmp(argv[2],"bf")) {
     m = PBQP_BRUTEFORCE;
   }

   /* read PBQP problem */ 
   pbqp *this = read_pbqp(f);

   /* set dump file */
   if (m == PBQP_HEURISTICAL && dump) {
      h_set_dumpfile(this,stderr);
   }

   /* solve PBQP problem */
   SOLVE_PBQP(m,this);

   /* print solutions of decision vectors */
   for(u=0;u<GET_NUMNODES(m,this);u++) {
     int sol=GET_SOLUTION(m,this,u);
     printf("%d\n",sol);
   }

   /* print minimum of objective function */
   printf("%6.0f\n",GET_MIN(m,this));

   /* is solution optimal? */
   printf("%d\n",IS_OPTIMAL(m,this));

   FREE_PBQP(m,this);

   return 0; 
}
