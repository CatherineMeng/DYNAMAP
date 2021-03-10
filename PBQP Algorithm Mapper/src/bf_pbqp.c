/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Brute Force PBQP Solver
 *****************************************************************************
 $Id: bf_pbqp.c,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $ 
 *****************************************************************************
 $Log: bf_pbqp.c,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "h_pbqp.h"


/**************************************************************************
 * Data Structures 
 **************************************************************************/

/* edge of PBQP graph */
typedef struct adjnode {
  struct adjnode *prev,      /* doubly chained list */ 
                 *succ;
  int adj;                   /* adj. node */
  mat *costs;                /* cost matrix of edge */
} adjnode;

/* data structure of partitioned boolean quadratic problem */
struct pbqp {
  int num_nodes;             /* number of nodes */
  int max_deg;               /* maximal degree of a node */
  boolean solved;            /* flag that indicates whether PBQP has been solved yet */
  num min;
                             /* node fields */
  vec **node_costs;	     /* cost vectors of nodes */
  int *solution;	     /* solution for node */
  int *min_solution;	     /* solution for node */
  adjnode **adj_list;        /* adj. list */
};


/*****************************************************************************
 * allocation/de-allocation of pbqp problem 
 ****************************************************************************/

/* allocate new partitioned booleanean quadratic program problem */
pbqp *bf_alloc_pbqp(int num_nodes)
{
  pbqp *this;
  int u;
  
  ASSERT(num_nodes > 0);
  
  /* allocate memory for pbqp data structure */   
  this = (pbqp *)alloc_mem(sizeof(pbqp));
  
  /* Initialize pbqp fields */
  this->num_nodes = num_nodes;
  this->solved = FALSE;
  this->min = 0;
  this->max_deg = 0;
  
  /* initialize/allocate node fields of pbqp */
  this->adj_list = (adjnode **) alloc_mem(sizeof(adjnode *)*num_nodes);
  this->solution = (int *) alloc_mem(sizeof(int)*num_nodes);
  this->min_solution = (int *) alloc_mem(sizeof(int)*num_nodes);
  this->node_costs = (vec **) alloc_mem(sizeof(vec *)*num_nodes);
  for(u=0;u<num_nodes;u++) {
    this->solution[u]=-1;
    this->min_solution[u]=-1;
    this->adj_list[u]=NULL;
    this->node_costs[u]=NULL;
  }
  
  return this;
}

/* free pbqp problem */
void bf_free_pbqp(pbqp *this)
{
  int u;
  adjnode *adj_ptr,*adj_next;
  
  ASSERT(this != NULL);
  
  /* free node cost fields */
  for(u=0;u < this->num_nodes;u++) {
    v_free(this->node_costs[u]);
  }
  free_mem(this->node_costs);
  
  /* free adj. list */
  ASSERT(this->adj_list != NULL);
  for(u=0;u < this->num_nodes; u++) {
    for(adj_ptr = this->adj_list[u]; adj_ptr != NULL; adj_ptr = adj_next) {
      adj_next = adj_ptr -> succ;
      if (u < adj_ptr->adj) {
	ASSERT(adj_ptr != NULL);
	m_free(adj_ptr->costs);
      }
      free(adj_ptr);
    }
  }
  free_mem(this->adj_list);
  
  /* free other node fields */
  free_mem(this->solution);
  free_mem(this->min_solution);

  /* free pbqp data structure itself */
  free_mem(this);
}


/****************************************************************************
 * adj. node routines 
 ****************************************************************************/

/* find data structure of adj. node of a given node */
static
adjnode *find_adjnode(pbqp *this,int u,int v)
{
  adjnode *adj_ptr;
  
  ASSERT (this != NULL);
  ASSERT (u >= 0 && u < this->num_nodes);
  ASSERT (v >= 0 && v < this->num_nodes);
  ASSERT(this->adj_list != NULL);

  for(adj_ptr = this -> adj_list[u];adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
    if (adj_ptr->adj == v) {
      return adj_ptr;
    }
  }
  return NULL;
}

/* allocate a new data structure for adj. node */
static
adjnode *alloc_adjnode(pbqp *this,int u, mat *costs)
{
  ASSERT(this != NULL);
  ASSERT(costs != NULL);
  ASSERT(u >= 0 && u < this->num_nodes);

  adjnode *p = (adjnode *)alloc_mem(sizeof(adjnode));
  ASSERT(p != NULL);
  
  p->adj = u;
  p->costs = costs;  
  return p;
}

/* insert adjacence node to adj. list */
static
void insert_adjnode(pbqp *this, int u, adjnode *adj_ptr)
{

  ASSERT(this != NULL);
  ASSERT(adj_ptr != NULL);
  ASSERT(u >= 0 && u < this->num_nodes);

  /* if adjacency list of node is not empty -> update
     first node of the list */
  if (this -> adj_list[u] != NULL) {
    ASSERT(this->adj_list[u]->prev == NULL);
    this->adj_list[u] -> prev = adj_ptr;
  }

  /* update doubly chained list pointers of pointers */
  adj_ptr -> succ = this->adj_list[u];
  adj_ptr -> prev = NULL;

  /* update adjacency list pointer of node u */
  this->adj_list[u] = adj_ptr;
}

/* get cost matrix ptr */
static
mat *get_costmatrix_ptr(pbqp *this, int u, int v)
{
  adjnode *adj_ptr;
  mat *m = NULL;

  ASSERT (this != NULL);
  ASSERT (u >= 0 && u < this->num_nodes);
  ASSERT (v >= 0 && v < this->num_nodes); 

  adj_ptr = find_adjnode(this,u,v);

  if (adj_ptr != NULL) {
    m = adj_ptr -> costs;
  } 

  return m;
}

/*****************************************************************************
 * edge functions
 ****************************************************************************/

/* insert edge to graph */
/* (does not check whether edge exists in graph */
static
void insert_edge(pbqp *this, int u, int v, mat *costs)
{
  adjnode *adj_u,
          *adj_v;
  
  /* create adjanceny entry for u */
  adj_u = alloc_adjnode(this,v,costs);
  insert_adjnode(this,u,adj_u);


  /* create adjanceny entry for v */
  adj_v = alloc_adjnode(this,u,costs);
  insert_adjnode(this,v,adj_v);
  
}


/*****************************************************************************
 * cost functions 
 ****************************************************************************/

/* Note: Since cost(u,v) = transpose(cost(v,u)), it would be necessary to store 
   two matrices for both edges (u,v) and (v,u). However, we only store the 
   matrix for the case u < v. For the other case we transpose the stored matrix
   if required. 
*/

/* add costs to cost vector of a node */
void bf_add_nodecosts(pbqp *this,int u, vec *costs)
{
  ASSERT(this != NULL);
  ASSERT(costs != NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  
  if (!this->node_costs[u]) {
    this->node_costs[u] = v_copy(costs);
  } else {
    v_add(this->node_costs[u],costs);
  }
}


/* add costs to cost matrix of an edge */
void bf_add_edgecosts(pbqp *this,int u,int v,mat *costs)
{
  mat *adj_costs;

  ASSERT(this!= NULL);
  ASSERT(costs != NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  ASSERT(v >= 0 && v <= this->num_nodes);
  
  /* does the edge u-v exists ? */
  if (u == v) {
    vec *diag = v_diagonalize(costs);
    bf_add_nodecosts(this,v,diag);
    v_free(diag);
  } else if ((adj_costs = get_costmatrix_ptr(this,u,v))!=NULL) {
    if ( u < v) {
      m_add(adj_costs,costs);
    } else {
      m_addtransposed(adj_costs,costs);
    }
  } else {
    adj_costs = ((u < v) ? m_copy(costs) : m_transpose(costs));
    insert_edge(this,u,v,adj_costs);
  } 
}

/*****************************************************************************
 * check pbqp problem
 ****************************************************************************/
static
void check_pbqp(pbqp *this)
{
  int u,v;
  mat *costs;
  adjnode *adj_ptr;
  
  ASSERT( this != NULL);
  
  for(u=0;u< this->num_nodes; u++) {
    ASSERT (this -> node_costs[u] != NULL);
    for(adj_ptr = this -> adj_list[u];adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
      v = adj_ptr -> adj;
      ASSERT( v>= 0 && v < this->num_nodes);
      if (u < v ) {
	costs = adj_ptr -> costs;
	ASSERT( costs -> rows == this->node_costs[u] -> len &&
		costs -> cols == this->node_costs[v] -> len);
      }           
    }
  }
}

/*****************************************************************************
 * PBQP solve routines 
 ****************************************************************************/

static
void next_solution(pbqp *this)
{
  int u;

  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);
  ASSERT(this->node_costs != NULL);
 
  for(u=0;u<this->num_nodes;u++) {
    ASSERT(this->solution[u] >= 0);
    if (this->node_costs[u]->len != 1 ) {
      this->solution[u]++;
      if (this->solution[u] < this->node_costs[u]->len ) {
         break;
      } else {
         this->solution[u] = 0;
      }
    }
  }
}

static
boolean is_firstsolution(pbqp *this)
{
  int u;
  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);
  ASSERT(this->node_costs != NULL);
 
  for(u=0;u < this->num_nodes;u++) {
    int sol_u = this->solution[u];
    ASSERT(sol_u >= 0 && sol_u < this->node_costs[u]->len);
    if(sol_u > 0) {
      return FALSE;
    }
  }
  return TRUE;
}

static void store_solution(pbqp *this)
{
  int u;
  
  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);
  ASSERT(this->min_solution != NULL);

  for(u=0;u<this->num_nodes;u++) {
     this->min_solution[u] = this->solution[u];
  }
}

static void restore_solution(pbqp *this)
{
  int u;
  
  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);
  ASSERT(this->min_solution != NULL);

  for(u=0;u<this->num_nodes;u++) {
     this->solution[u] = this->min_solution[u];
  }
}

/* get value of pbqp function */
num get_h(pbqp *this)
{
  num h=0.0;
  int u;
  adjnode *adj_ptr;
  
  ASSERT(this->adj_list != NULL);

  for(u=0;u<this->num_nodes;u++) {
    int sol_u = this->solution[u];

    ASSERT(sol_u >= 0 && sol_u < this->node_costs[u]->len);
    h+= V_ELEM(this->node_costs[u],sol_u);

    for(adj_ptr = this -> adj_list[u];adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
      int v = adj_ptr -> adj;
      int sol_v = this->solution[v];
      ASSERT(sol_v >= 0 && sol_v < this->node_costs[v]->len);
      if (u < v) {
         h+=M_ELEM(adj_ptr->costs,sol_u,sol_v); 
      }
    }
  }
  return h;
}

/* solve PBQP problem */
void bf_solve_pbqp(pbqp *this)
{
  num h;
  int u;
  
  ASSERT(this != NULL);
  ASSERT(!this->solved); 
  ASSERT(this->solution != NULL);
  
  /* check vector & matrix dimensions */
  check_pbqp(this);

  /* get minimum of first solution */
  for(u=0;u<this->num_nodes;u++) {
    this->solution[u] = 0;
  }
  this->min = get_h(this);
  store_solution(this);

  /* iterate solution space and find */
  /* cheapest solution */
  for(;;) {
    next_solution(this);

    /* stop loop until first solution appears again */
    if(is_firstsolution(this)) {
       break; 
    }

    /* determine costs of next solution */
    h = get_h(this);

    if (h < this->min)  {
       this->min = h;
       store_solution(this);
    }
  }  
 
  /* restore minimal solution of PBQP */
  restore_solution(this);
  
  this->solved = TRUE;
}

/* validate solution */
boolean bf_validate_pbqp(pbqp *this)
{
   return (get_h(this) == this->min)?TRUE:FALSE;
}

/* set solution of a node */
void bf_set_solution(pbqp *this,int x,int sol)
{
  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);

  ASSERT(sol >= 0 && sol < this->node_costs[x]->len);
  
  this->solution[x]=sol;
}

/* set min */
void bf_set_min(pbqp *this,num min)
{
  ASSERT(this != NULL);
  
  this->min=min;
}


/*****************************************************************************
 * getter functions
 ****************************************************************************/

/* get solution of a node */
int bf_get_solution(pbqp *this,int x)
{
  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);
  ASSERT(this -> solved);
  
  return this->solution[x];
}

/* get number of nodes */
int bf_get_numnodes(pbqp *this)
{
  return this->num_nodes;
}

/* get minimum */
num bf_get_min(pbqp *this)
{
  return this->min;
}

/* is optimal? */
boolean bf_is_optimal(pbqp *this)
{
  return TRUE;
}


/* end of bf_pbqp.h */
