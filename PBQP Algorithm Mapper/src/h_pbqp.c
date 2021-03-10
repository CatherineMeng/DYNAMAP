/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Heuristical PBQP Solver
 *****************************************************************************
 $Id: h_pbqp.c,v 1.2 2003/12/30 13:27:26 scholz Exp $ 
 *****************************************************************************
 $Log: h_pbqp.c,v $
 Revision 1.2  2003/12/30 13:27:26  scholz
 changed text of dump output

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
                 *succ, 
                 *reverse;   /* reverse edge */
  int adj;                   /* adj. node */
  mat *costs;                /* cost matrix of edge */
} adjnode;

/* bucket node */
typedef struct bucketnode {
  struct bucketnode *prev;   /* doubly chained list */
  struct bucketnode *succ;   
  int u;                     /* node */
} bucketnode;

/* data structure of partitioned boolean quadratic problem */
struct pbqp {
  int num_nodes;             /* number of nodes */
  int max_deg;               /* maximal degree of a node */
  boolean solved;            /* flag that indicates whether PBQP has been solved yet */
  boolean optimal;           /* flag that indicates whether PBQP is optimal */
  num min;
  boolean changed;           /* flag whether graph has changed in simplification */

                             /* node fields */
  vec **node_costs;	     /* cost vectors of nodes */
  int *node_deg;	     /* node degree of nodes */
  int *solution;	     /* solution for node */
  adjnode **adj_list;        /* adj. list */
  bucketnode **bucket_ptr;   /* bucket pointer of a node */

                             /* node stack */
  int *stack;	             /* stack of nodes */
  int stack_ptr;             /* stack pointer */

                             /* bucket fields */
  bucketnode **bucket_list;  /* bucket list */

  FILE *dump_file;           /* dump file */
};

/*****************************************************************************
 * prototypes of local function
 *****************************************************************************/

static
void remove_bucket(pbqp *this, bucketnode *bucket);

/*****************************************************************************
 * allocation/de-allocation of pbqp problem 
 ****************************************************************************/

/* allocate new partitioned booleanean quadratic program problem */
pbqp *h_alloc_pbqp(int num_nodes)
{
  pbqp *this;
  int u;
  
  ASSERT(num_nodes > 0);
  
  /* allocate memory for pbqp data structure */   
  this = (pbqp *)alloc_mem(sizeof(pbqp));
  
  /* Initialize pbqp fields */
  this->num_nodes = num_nodes;
  this->solved = FALSE;
  this->optimal = TRUE;
  this->min = 0;
  this->max_deg = 0;
  this->dump_file = NULL;
  this->changed = FALSE;
  
  /* initialize/allocate stack fields of pbqp */ 
  this->stack = (int *) alloc_mem(sizeof(int)*num_nodes);
  this->stack_ptr = 0;
  
  /* initialize/allocate node fields of pbqp */
  this->adj_list = (adjnode **) alloc_mem(sizeof(adjnode *)*num_nodes);
  this->node_deg = (int *) alloc_mem(sizeof(int)*num_nodes);
  this->solution = (int *) alloc_mem(sizeof(int)*num_nodes);
  this->bucket_ptr = (bucketnode **) alloc_mem(sizeof(bucketnode **)*num_nodes);
  this->node_costs = (vec **) alloc_mem(sizeof(vec *)*num_nodes);
  for(u=0;u<num_nodes;u++) {
    this->solution[u]=-1;
    this->adj_list[u]=NULL;
    this->node_deg[u]=0;
    this->bucket_ptr[u]=NULL;
    this->node_costs[u]=NULL;
  }
  
  /* initialize bucket list */
  this->bucket_list = NULL;
  
  return this;
}

/* free pbqp problem */
void h_free_pbqp(pbqp *this)
{
  int u;
  int deg;
  adjnode *adj_ptr,*adj_next;
  bucketnode *bucket,*bucket_next;
  
  ASSERT(this != NULL);
  
  /* free node cost fields */
  for(u=0;u < this->num_nodes;u++) {
    v_free(this->node_costs[u]);
  }
  free_mem(this->node_costs);
  
  /* free bucket list */
  for(deg=0;deg<this->max_deg;deg++) {
    for(bucket=this->bucket_list[deg];bucket!=NULL;bucket=bucket_next) {
      this->bucket_ptr[bucket->u] = NULL;
      bucket_next = bucket-> succ;
      free_mem(bucket);
    }
  }
  free_mem(this->bucket_list);
  
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
  free_mem(this->node_deg);
  free_mem(this->solution);
  free_mem(this->bucket_ptr);

  /* free stack */
  free_mem(this->stack);

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

/* remove entry in an adj. list */
static
void remove_adjnode(pbqp *this, int u, adjnode *adj_ptr)
{
  ASSERT(this!= NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  ASSERT(this->adj_list != NULL);
  ASSERT(adj_ptr != NULL);
  
  if (adj_ptr -> prev == NULL) {
    this->adj_list[u] = adj_ptr -> succ;
  } else {
    adj_ptr -> prev -> succ = adj_ptr -> succ;
  } 

  if (adj_ptr -> succ != NULL) {
    adj_ptr -> succ -> prev = adj_ptr -> prev;
  }

  if(adj_ptr->reverse != NULL) {
    adjnode *rev = adj_ptr->reverse;
    rev->reverse = NULL;
  }

  free_mem(adj_ptr);
}

/*****************************************************************************
 * node functions 
 ****************************************************************************/

/* get degree of a node */
static
int get_deg(pbqp *this,int u)
{
  adjnode *adj_ptr;
  int deg = 0;
  
  ASSERT(this != NULL);
  ASSERT(u >= 0 && u < this->num_nodes);
  ASSERT(this->adj_list != NULL);

  for(adj_ptr = this -> adj_list[u];adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
    deg ++;
  }
  return deg;
}

/* reinsert node */
static
void reinsert_node(pbqp *this,int u)
{
  adjnode *adj_u,
          *adj_v;

  ASSERT(this!= NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  ASSERT(this->adj_list != NULL);

  for(adj_u = this -> adj_list[u]; adj_u != NULL; adj_u = adj_u -> succ) {
    int v = adj_u -> adj;
    adj_v = alloc_adjnode(this,u,adj_u->costs);
    insert_adjnode(this,v,adj_v);
  }
}

/* remove node */
static
void remove_node(pbqp *this,int u)
{
  adjnode *adj_ptr;

  ASSERT(this!= NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  ASSERT(this->adj_list != NULL);

  for(adj_ptr = this -> adj_list[u]; adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
    remove_adjnode(this,adj_ptr->adj,adj_ptr -> reverse);
  }
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
  
  /* create link for reverse edge */
  adj_u -> reverse = adj_v;
  adj_v -> reverse = adj_u;
}

/* delete edge */
static
void delete_edge(pbqp *this,int u,int v)
{
  adjnode *adj_ptr;
  adjnode *rev;
  
  ASSERT(this != NULL);
  ASSERT( u >= 0 && u < this->num_nodes);
  ASSERT( v >= 0 && v < this->num_nodes);

  adj_ptr=find_adjnode(this,u,v);
  ASSERT(adj_ptr != NULL);
  ASSERT(adj_ptr->reverse != NULL);

  m_free(adj_ptr -> costs);
 
  rev = adj_ptr->reverse; 
  remove_adjnode(this,u,adj_ptr);
  remove_adjnode(this,v,rev);
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
void h_add_nodecosts(pbqp *this,int u, vec *costs)
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

/* get cost matrix ptr */
/* Note: only the pointer is returned for 
   cost(u,v), if u < v.
*/ 
static
mat *get_costmatrix(pbqp *this, int u, int v)
{
  adjnode *adj_ptr = find_adjnode(this,u,v);
  
  if (adj_ptr != NULL) {
    if ( u < v) {
      return m_copy(adj_ptr -> costs);
    } else {
      return m_transpose(adj_ptr -> costs);
    }
  } else {
    return NULL;
  }  
}

/* add costs to cost matrix of an edge */
void h_add_edgecosts(pbqp *this,int u,int v,mat *costs)
{
  mat *adj_costs;

  ASSERT(this!= NULL);
  ASSERT(costs != NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  ASSERT(v >= 0 && v <= this->num_nodes);
  
  /* does the edge u-v exists ? */
  if (u == v) {
    vec *diag = v_diagonalize(costs);
    h_add_nodecosts(this,v,diag);
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

/**********************************************************************************
 * pop functions
 **********************************************************************************/

/* pop node of certain degree */
static
int pop_node(pbqp *this,int deg)
{
  bucketnode *bucket;
  int u;

  ASSERT(this != NULL);
  ASSERT(deg >= 0 && deg <= this->max_deg);
  ASSERT(this->bucket_list != NULL);
   
  /* get first bucket of bucket list */
  bucket = this->bucket_list[deg];
  ASSERT(bucket != NULL);

  /* remove bucket */
  remove_bucket(this,bucket);
  u = bucket->u;
  free_mem(bucket);
  return u;
}

/* pop node with maximal degree */
static
int pop_maxnode(pbqp *this)
{
  int deg;
  bucketnode *bucket;
  int u;
  
  for(deg=this->max_deg;deg > 2;deg--) {
    if (this->bucket_list[deg] != NULL) {
      bucket = this->bucket_list[deg];
      remove_bucket(this,bucket);
      u = bucket->u;
      free_mem(bucket);
      return u;
    }
  }
  return -1;
}

/**********************************************************************************
 * reorder functions
 **********************************************************************************/

/* add bucket to bucketlist */
static
void add_to_bucketlist(pbqp *this,bucketnode *bucket, int deg)
{
  bucketnode *old_head;
  
  ASSERT(bucket != NULL);
  ASSERT(this != NULL);
  ASSERT(deg >= 0 && deg <= this->max_deg);
  ASSERT(this->bucket_list != NULL);

  /* store node degree (for re-ordering purposes)*/
  this->node_deg[bucket->u] = deg;
  
  /* put bucket to front of doubly chained list */
  old_head = this->bucket_list[deg];
  bucket -> prev = NULL;
  bucket -> succ = old_head;
  this -> bucket_list[deg] = bucket;
  if (bucket -> succ != NULL ) {
    ASSERT ( old_head -> prev == NULL);
    old_head -> prev = bucket;
  }
}

/* remove bucket from bucket list */
static
void remove_bucket(pbqp *this, bucketnode *bucket)
{
  int u = bucket->u;
  
  ASSERT(this != NULL);
  ASSERT(u >= 0 && u < this->num_nodes);
  ASSERT(this->bucket_list != NULL);
  ASSERT(this->bucket_ptr[u] != NULL);
  
  /* update predecessor node in bucket list 
     (if no preceeding bucket exists, then
     the bucket_list pointer needs to be 
     updated.)
  */    
  if (bucket->prev != NULL) {
    bucket->prev-> succ = bucket->succ; 
  } else {
    this->bucket_list[this->node_deg[u]] = bucket -> succ;
  }
  
  /* update successor node in bucket list */ 
  if (bucket->succ != NULL) { 
    bucket->succ-> prev = bucket->prev;
  }
}

/* reorder node in bucket list according to 
   current node degree */
static
void reorder_node(pbqp *this, int u)
{
  int deg; 
  
  ASSERT(this != NULL);
  ASSERT(u>= 0 && u < this->num_nodes);
  ASSERT(this->bucket_list != NULL);
  ASSERT(this->bucket_ptr[u] != NULL);

  /* get current node degree */
  deg = get_deg(this,u);
  
  /* remove bucket from old bucket list only
     if degree of node has changed. */
  if (deg != this->node_deg[u]) {
    remove_bucket(this,this->bucket_ptr[u]);
    add_to_bucketlist(this,this->bucket_ptr[u],deg);
  } 
}

/* reorder adj. nodes of a node */
static
void reorder_adjnodes(pbqp *this,int u)
{
  adjnode *adj_ptr;
  
  ASSERT(this!= NULL);
  ASSERT(u >= 0 && u <= this->num_nodes);
  ASSERT(this->adj_list != NULL);

  for(adj_ptr = this -> adj_list[u]; adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
    reorder_node(this,adj_ptr->adj);
  }
}

/**********************************************************************************
 * creation functions
 **********************************************************************************/

/* create new bucket entry */
/* consistency of the bucket list is not checked! */
static
void create_bucket(pbqp *this,int u,int deg)
{
  bucketnode *bucket;
  
  ASSERT(this != NULL);
  ASSERT(u >= 0 && u < this->num_nodes);
  ASSERT(this->bucket_list != NULL);
  
  bucket = (bucketnode *)alloc_mem(sizeof(bucketnode));
  ASSERT(bucket != NULL);

  bucket -> u = u;
  this->bucket_ptr[u] = bucket;

  add_to_bucketlist(this,bucket,deg);
}

/* create bucket list */
static
void create_bucketlist(pbqp *this)
{
  int u;
  int max_deg;
  int deg;

  ASSERT(this != NULL);
  ASSERT(this->bucket_list == NULL);

  /* determine max. degree of the nodes */
  max_deg = 2;  /* at least of degree two! */
  for(u=0;u<this->num_nodes;u++) {
    deg = this->node_deg[u] = get_deg(this,u);
    if (deg > max_deg) {
      max_deg = deg;
    }
  }
  this->max_deg = max_deg;
  
  /* allocate bucket list */
  this -> bucket_list = (bucketnode **)alloc_mem(sizeof(bucketnode *)*(max_deg + 1));
  memset(this->bucket_list,0,sizeof(bucketnode *)*(max_deg + 1));
  ASSERT(this->bucket_list != NULL);
  
  /* insert nodes to the list */
  for(u=0;u<this->num_nodes;u++) {
    create_bucket(this,u,this->node_deg[u]);  
  }
}

/*****************************************************************************
 * dump routines 
 ****************************************************************************/

static
void dump_reductiongraph(pbqp *this,int x)
{
  int deg;
  bucketnode *bucket;
  adjnode *adj_ptr;
 
  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  fputs("<p>\n<graph>\n\tgraph input {\n",this->dump_file);
  fprintf(this->dump_file,"\t n%d [color=red];\n",x);

  /* dump active nodes and denote node which is going to be 
     reduced. */
  for(deg=1;deg<this->max_deg;deg++) {
    for(bucket=this->bucket_list[deg];bucket!=NULL;bucket=bucket -> succ) {
      int u = bucket -> u;
      fprintf(this->dump_file,"\t n%d;\n",u);
    }
  }

  /* dump active edges and denote edges which are going to be 
     reduced. */
  for(adj_ptr=this->adj_list[x];adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
     int v = adj_ptr -> adj;
     fprintf(this->dump_file,"\t n%d--n%d;\n",x,v);
  }
  for(deg=1;deg<=this->max_deg;deg++) {
    for(bucket=this->bucket_list[deg];bucket!=NULL;bucket=bucket -> succ) {
      int u = bucket -> u;
      for(adj_ptr=this->adj_list[u];adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
         int v = adj_ptr -> adj;
         fprintf(this->dump_file,"\t n%d--n%d;\n",u,v);
      }
    }
  }
  fputs("\t}\n</graph>\n</p>\n",this->dump_file);
}


static
void dump_simplifyedge(pbqp *this,int u,int v)
{
  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  ASSERT(u >= 0 && u < this->num_nodes);
  ASSERT(v >= 0 && v < this->num_nodes);

  mat *costs = get_costmatrix(this,u,v);
  
  fputs("<p>\n",this->dump_file);

  fputs("<tex>\n",this->dump_file);
  fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",u);
  v_texprint(this->dump_file,this->node_costs[u]);
  fputs("</tex><br>\n",this->dump_file);

  if (costs != NULL) {
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\overline{C}_{%d,%d}=\n",u,v);
    m_texprint(this->dump_file,costs);
    fputs("</tex><br>",this->dump_file);
  }
   
  fputs("<tex>\n",this->dump_file);
  fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",v);
  v_texprint(this->dump_file,this->node_costs[v]);
  fputs("</tex><br>\n",this->dump_file);

  if(costs == NULL) {
    fputs("edge has been eliminated\n",this->dump_file);
  }
  
  fputs("</p>\n",this->dump_file);

  if (costs != NULL) {
    m_free(costs);
  }
}

static
void dump_header(pbqp *this)
{
  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  fputs("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"
        "<html>\n<head>\n<meta http-equiv=\"content-type\" "
        "content=\"text/html; charset=ISO-8859-1\">\n<title>PBQP Dump</title>\n"
        "</head>\n<body>\n",this->dump_file);

}

static 
void dump_footer(pbqp *this)
{
  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  fputs("</body>\n</html>\n",this->dump_file);
}

static 
void dump_section(FILE *f,int level,char *txt)
{
  ASSERT(f!=NULL);

  fprintf(f,"<h%d>%s</h%d>\n",level,txt,level);
}

static 
void dump_graph(pbqp *this)
{
  int u;
  adjnode *adj_ptr;

  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  fputs("<p>\n<graph>\n\tgraph input {\n",this->dump_file);
  ASSERT(this->adj_list != NULL);
  for(u=0;u < this->num_nodes; u++) {
    fprintf(this->dump_file,"\t n%d;\n",u);
  }
  for(u=0;u < this->num_nodes; u++) {
    for(adj_ptr = this->adj_list[u]; adj_ptr != NULL; adj_ptr = adj_ptr->succ) {
      int v = adj_ptr -> adj;
      if (u < v) {
         fprintf(this->dump_file,"\t n%d -- n%d;\n",u,v);
      }
    }
  }
  fputs("\t}\n</graph>\n</p>\n",this->dump_file);
}

static 
void dump_nodecosts(pbqp *this)
{
  int u;

  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  /* dump node costs */ 
  fputs("<p>",this->dump_file);
  for(u=0;u < this->num_nodes; u++) {
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",u);
    v_texprint(this->dump_file,this->node_costs[u]);
    fputs("</tex><br>\n",this->dump_file);
  }
  fputs("</p>",this->dump_file);
}

static
void dump_edgecosts(pbqp *this)
{
  int u;
  adjnode *adj_ptr;

  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  fputs("<p>",this->dump_file);
  for(u=0;u < this->num_nodes; u++) {
    for(adj_ptr = this->adj_list[u]; adj_ptr != NULL; adj_ptr = adj_ptr->succ) {
      int v = adj_ptr -> adj;
      if (u < v) {
        fputs("<tex>\n",this->dump_file);
        fprintf(this->dump_file,"\t\\overline{C}_{%d,%d}=\n",u,v);
        m_texprint(this->dump_file,adj_ptr->costs);
        fputs("</tex><br>",this->dump_file);
      }
    }
  }
  fputs("</p>",this->dump_file);
}

static
void dump_trivialnodes(pbqp *this)
{
  int u;
  int ctr;

  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  dump_section(this->dump_file,1,"2. Remove Trivial Nodes");
  for(ctr=0,u=0;u < this->num_nodes; u++) {
    if(this->node_costs[u]->len == 1) {
      ctr++;
    } 
  }
  if (ctr > 0) {
    fputs("<p>Following nodes have cost vectors of length one:</p>",this->dump_file);
    fputs("<ul type=\"square\">",this->dump_file);
    for(u=0;u < this->num_nodes; u++) {
      if(this->node_costs[u]->len == 1) {
        fprintf(this->dump_file,"<li>n%d</li>\n",u);
      } 
    }
    fputs("</ul>\n",this->dump_file);
    dump_section(this->dump_file,2,"2.1 New Topology");
    dump_graph(this);
    dump_section(this->dump_file,2,"2.2 New Cost Vectors");
    dump_nodecosts(this);
  } else {
    fputs("<p>no trivial nodes<br>",this->dump_file);
  }
}

static
void dump_input(pbqp *this)
{
  ASSERT(this != NULL);
  ASSERT(this->dump_file != NULL);

  dump_section(this->dump_file,1,"1. PBQP Problem");
  dump_section(this->dump_file,2,"1.1 Topology");
  dump_graph(this);
  dump_section(this->dump_file,2,"1.2 Cost Vectors");
  dump_nodecosts(this);
  dump_section(this->dump_file,2,"1.3 Cost Matrices");
  dump_edgecosts(this);
}

/*****************************************************************************
 * PBQP simplification for trivial nodes
 ****************************************************************************/

/* remove trivial node with cost vector length of one */
static
void disconnect_trivialnode(pbqp *this,int u)
{
  int v;
  adjnode *adj_ptr, 
          *next;
  mat *c_uv;
  vec *c_v;
  
  ASSERT(this != NULL);
  ASSERT(this->node_costs != NULL);
  ASSERT(u >= 0 && u < this -> num_nodes);
  ASSERT(this->node_costs[u]->len == 1);
  
  /* add edge costs to node costs of adj. nodes */
  for(adj_ptr = this->adj_list[u]; adj_ptr != NULL; adj_ptr = next){
    next = adj_ptr -> succ;
    v = adj_ptr -> adj;
    ASSERT(v >= 0 && v < this -> num_nodes);
    
    /* convert matrix to cost vector offset for adj. node */
    c_uv = get_costmatrix(this,u,v);
    c_v = v_row(c_uv,0);
    v_add(this->node_costs[v],c_v);
    
    /* delete edge & free vec/mat */
    v_free(c_v);
    m_free(c_uv);
    delete_edge(this,u,v);
  }   
}

/* find all trivial nodes and disconnect them */
static   
void eliminate_trivial_nodes(pbqp *this)
{
   int u;
   
   ASSERT(this != NULL);
   ASSERT(this -> node_costs != NULL);
   
   for(u=0;u < this -> num_nodes; u++) {
     if (this->node_costs[u] -> len == 1) {
       disconnect_trivialnode(this,u); 
     }
   }
}

/*****************************************************************************
 * Normal form for PBQP 
 ****************************************************************************/

/* simplify a cost matrix. If the matrix
   is independent, then simplify_matrix
   returns TRUE - otherwise FALSE. In
   vectors u and v the offset values of
   the decomposition are stored. 
*/

static
boolean normalize_matrix(mat *m,vec *u,vec *v)
{
  int rows, 
      cols;
  int r,
      c;
  
  ASSERT( m != NULL);
  ASSERT( u != NULL);
  ASSERT( v != NULL);
  ASSERT( u->len > 0);
  ASSERT( v->len > 0);
  
  rows = m->rows;
  cols = m->cols;
  
  ASSERT(rows == u->len);
  ASSERT(cols == v->len);

  /* determine u vector */
  for(r=0;r<rows;r++) {
    num min = m_rowmin(m,r);
    u->data[r] += min;
    if (!isinf(min)) {
      m_subrow(m,r,min);
    } else {
      m_setrow(m,r,0);
    }
  }
  
  /* determine v vector */
  for(c=0;c<cols;c++) {
    num min = m_colmin(m,c);
    v -> data[c] += min;
    if (!isinf(min)){
      m_subcol(m,c,min);
    } else {
      m_setcol(m,c,0);
    }
  }
  
  /* determine whether matrix is 
     independent or not. 
    */
  return m_iszero(m);
}

/* simplify single edge */
static
void simplify_edge(pbqp *this,int u,int v)
{
  mat *costs;
  boolean is_zero; 
  
  ASSERT (this != NULL);
  ASSERT (u >= 0 && u <this->num_nodes);
  ASSERT (v >= 0 && v <this->num_nodes);
  ASSERT (u != v);

  /* swap u and v  if u > v in order to avoid un-necessary
     tranpositions of the cost matrix */
  
  if (u > v) {
    int swap = u;
    u = v;
    v = swap;
  }
  
  /* get cost matrix and simplify it */  
  costs = get_costmatrix_ptr(this,u,v);
  is_zero=normalize_matrix(costs,this->node_costs[u],this->node_costs[v]);

  /* delete edge */
  if(is_zero){
    delete_edge(this,u,v);
    this->changed = TRUE;
  }
}

/* normalize cost matrices and remove 
   edges in PBQP if they ary independent, 
   i.e. can be decomposed into two 
   cost vectors.
*/
static
void eliminate_independent_edges(pbqp *this)
{
  int u,v;
  adjnode *adj_ptr;
  int step = 1;
  
  ASSERT(this != NULL);
  ASSERT(this -> adj_list != NULL);

  this->changed = FALSE;
  if (this->dump_file) {
     dump_section(this->dump_file,1,"3. Simplification of Cost Matrices");
     dump_section(this->dump_file,2,"3.1. New Vectors/Matrices");
  }
  for(u=0;u < this->num_nodes;u++) {
    for (adj_ptr = this -> adj_list[u]; adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
      v = adj_ptr -> adj;
      ASSERT(v >= 0 && v < this->num_nodes);
      if (u < v) {
        if (this->dump_file) {
           static char txt[100];
           sprintf(txt,"3.2.%d Simplification of Edge n%d-n%d",step++,u,v);
           dump_section(this->dump_file,3,txt);
           fprintf(this->dump_file,"<p>Input:</p>");
           dump_simplifyedge(this,u,v); 
        }
	simplify_edge(this,u,v);
        if (this->dump_file) {
           fprintf(this->dump_file,"<p>Output:</p>");
           dump_simplifyedge(this,u,v); 
        }
      } 
    }
  }
  if (this->dump_file) {
     if (this->changed) {
        dump_section(this->dump_file,2,"3.2. New Topology");
        dump_graph(this);
     }
  }
}


/*****************************************************************************
 * PBQP reduction rules 
 ****************************************************************************/

/* RI reduction
   This reduction rule is applied for nodes 
   of degree one. */

static
void apply_RI(pbqp *this,int x)
{
  int y,
      i,j; 
  int xlen;
  int ylen;
  mat *c_yx;
  vec *c_x, *delta;
  
  ASSERT(this != NULL);
  ASSERT(x >= 0 && x < this->num_nodes);
  ASSERT(this -> adj_list[x] != NULL);
  ASSERT(this -> adj_list[x] -> succ == NULL);

  /* get adjacence matrix */
  y = this -> adj_list[x] -> adj;
  ASSERT(y >= 0 && y < this->num_nodes);
  
  /* determine length of cost vectors for node x and y */
  xlen = this -> node_costs[x]->len;
  ylen = this -> node_costs[y]->len;

  /* get cost vector c_x and matrix c_yx */
  c_x = this -> node_costs[x];
  c_yx = get_costmatrix(this,y,x); 
  ASSERT (c_yx != NULL);

  if (this->dump_file) {
    fputs("<p>Before reduction:<br>\n<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",y);
    v_texprint(this->dump_file,this->node_costs[y]);
    fputs("</tex></p>\n",this->dump_file);
    fputs("<p>Input:<br>\n<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",x);
    v_texprint(this->dump_file,this->node_costs[x]);
    fputs("</tex><br>\n",this->dump_file);
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\overline{C}_{%d,%d}=\n",y,x);
    m_texprint(this->dump_file,c_yx);
    fputs("</tex>\n</p>\n",this->dump_file);
  }
  
  /* allocate delta vector */
  delta = v_alloc(ylen);

  /* compute delta vector */
  for(i=0;i<ylen;i++) {
    num min =  M_ELEM(c_yx,i,0) + V_ELEM(c_x,0);
    for(j=1;j<xlen;j++) {
      num c =  M_ELEM(c_yx,i,j) + V_ELEM(c_x,j);
      if ( c < min )  
         min = c;
    }
    v_set(delta,i,min); 
  } 

  /* add delta vector */
  v_add(this -> node_costs[y],delta);

  /* dump delta vector for cost vector y */
  if (this->dump_file) {
    fputs("<p>Output:<br>\n<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\vec{\\Delta}_{%d}=\n",y);
    v_texprint(this->dump_file,delta);
    fputs("</tex>\n<br>\n",this->dump_file);
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",y);
    v_texprint(this->dump_file,this->node_costs[y]);
    fputs("</tex>\n</p>\n",this->dump_file);
  }


  /* delete node x */
  remove_node(this,x);

  /* reorder adj. nodes of node x */
  reorder_adjnodes(this,x);

  /* push node x on stack */
  ASSERT(this -> stack_ptr < this -> num_nodes);
  this->stack[this -> stack_ptr++] = x;

  /* free vec/mat */
  m_free(c_yx);
  v_free(delta);
}

/* RII reduction
   This reduction rule is applied for nodes 
   of degree two. */

static
void apply_RII(pbqp *this,int x)
{
  int y,z; 
  int xlen,ylen,zlen;
  int i,j,k;

  mat *c_yx, *c_zx;
  vec *cx;
  mat *delta;
 
  ASSERT(this != NULL);
  ASSERT(x >= 0 && x < this->num_nodes);
  ASSERT(this -> adj_list[x] != NULL);
  ASSERT(this -> adj_list[x] -> succ != NULL);
  ASSERT(this -> adj_list[x] -> succ -> succ == NULL);

  /* get adjacence matrix */
  y = this -> adj_list[x] -> adj;
  z = this -> adj_list[x] -> succ -> adj;
  ASSERT(y >= 0 && y < this->num_nodes);
  ASSERT(z >= 0 && z < this->num_nodes);
  
  /* determine length of cost vectors for node x and y */
  xlen = this -> node_costs[x]->len;
  ylen = this -> node_costs[y]->len;
  zlen = this -> node_costs[z]->len;

  /* get cost vector c_x and matrix c_yx */
  cx = this -> node_costs[x];
  c_yx = get_costmatrix(this,y,x); 
  c_zx = get_costmatrix(this,z,x); 
  ASSERT(c_yx != NULL);
  ASSERT(c_zx != NULL);

  if (this->dump_file) {
    fputs("<p>Before reduction:</p>\n",this->dump_file);
    dump_simplifyedge(this,y,z);
    
    fputs("<p>Input:<br>\n",this->dump_file);
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\overline{C}_{%d,%d}=\n",y,x);
    m_texprint(this->dump_file,c_yx);
    fputs("</tex><br>\n",this->dump_file);
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\vec{c}_{%d}=\n",x);
    v_texprint(this->dump_file,this->node_costs[x]);
    fputs("</tex><br>\n",this->dump_file);
    fputs("<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\overline{C}_{%d,%d}=\n",z,x);
    m_texprint(this->dump_file,c_zx);
    fputs("</tex><br>\n",this->dump_file);
    fputs("</p>\n",this->dump_file);
  }
  /* allocate delta matrix */
  delta = m_alloc(ylen,zlen);

  /* compute delta matrix */
  for(i=0;i<ylen;i++) {
    for(j=0;j<zlen;j++) {
      num min = M_ELEM(c_yx,i,0) + M_ELEM(c_zx,j,0) + V_ELEM(cx,0);
      for(k=1;k<xlen;k++) {
        num c = M_ELEM(c_yx,i,k) + M_ELEM(c_zx,j,k) + V_ELEM(cx,k);
        if ( c < min ) {
          min = c;
        }
      }
      m_set(delta,i,j,min);
    }
  }

  /* add delta matrix */
  h_add_edgecosts(this,y,z,delta);

  /* dump delta matrix and new edge costs */
  if (this->dump_file) {
    mat *costs;

    fputs("<p>Output:<br>\n<tex>\n",this->dump_file);
    fprintf(this->dump_file,"\t\\overline{\\Delta}_{%d,%d}=\n",y,z);
    m_texprint(this->dump_file,delta);
    fputs("</tex>\n<br>\n",this->dump_file);
    fputs("<tex>\n",this->dump_file);
    costs = get_costmatrix(this,y,z); 
    fprintf(this->dump_file,"\t\\overline{C}_{%d,%d}=\n",y,z);
    m_texprint(this->dump_file,costs);
    fputs("</tex>\n<br>\n",this->dump_file);
    
    m_free(costs);
  }
      
  /* delete node x */
  remove_node(this,x);

  /* simplify cost matrix c_yz */
  simplify_edge(this,y,z);

  if (this->dump_file) {
    fputs("<p>Simplified:</p>\n",this->dump_file);
    dump_simplifyedge(this,y,z);
  }

  /* reorder adj. nodes */
  reorder_adjnodes(this,x);

  /* push node x on stack */
  ASSERT(this -> stack_ptr < this -> num_nodes);
  this->stack[this -> stack_ptr++] = x;

  /* free vec/mat */
  m_free(c_yx);
  m_free(c_zx);
  m_free(delta);
}

/* RN reduction */
static
void apply_RN(pbqp *this,int x)
{
  int sol,
      min_sol = 0,
      xlen;
  num min = 0;
  adjnode *adj_ptr;

  ASSERT(this != NULL);
  ASSERT(x >= 0 && x < this->num_nodes);
  ASSERT(this -> node_costs[x] != NULL);

  xlen = this -> node_costs[x] -> len;

  /* after application of RN rule no optimality
     can be guaranteed! */
  this -> optimal = FALSE;

  /* determine local minimum */
  for(sol=0;sol<xlen;sol++) {
    num h = V_ELEM(this -> node_costs[x],sol);
    for(adj_ptr=this->adj_list[x] ;adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
      int y = adj_ptr -> adj;
      mat *c_xy = get_costmatrix(this,x,y);
      vec *v = v_copy(this -> node_costs[y]);
     
      ASSERT(c_xy != NULL);
      v_addrow(v,c_xy,sol);
      h = h + v_min(v);

      v_free(v);
      m_free(c_xy);
    }
    if (h < min || sol == 0) {
      min = h;
      min_sol = sol;
    } 
  }
  ASSERT( min_sol >= 0 && min_sol < xlen);
  this->solution[x] = min_sol;

  if (this->dump_file) {
     fprintf(this->dump_file,"<p>Solution of n%d is %d.</p>",x,min_sol);
  }
  
  /* add solution costs to minimum */
  this->min += V_ELEM(this->node_costs[x],min_sol);
  
  /* add cost vectors to adj. nodes of node x */
  for(adj_ptr=this->adj_list[x] ;adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
    int y = adj_ptr -> adj;
    mat *c_xy = get_costmatrix(this,x,y);
    v_addrow(this->node_costs[y],c_xy,min_sol);
    m_free(c_xy);
  }
  

  /* push node x on stack */
  ASSERT(this -> stack_ptr < this -> num_nodes);
  this->stack[this -> stack_ptr++] = x;

  /* delete node x */ 
  remove_node(this,x);

  /* reorder adj. nodes of node x */
  reorder_adjnodes(this,x);
}

/*****************************************************************************
 * PBQP graph parsing
 ****************************************************************************/
 
/* reduce pbqp problem (first phase) */
static
void reduce_pbqp(pbqp *this)
{
  int u;
  int step=1;
  static char txt[100];
 
  ASSERT(this != NULL);
  ASSERT(this->bucket_list != NULL);

  if (this->dump_file) {
     dump_section(this->dump_file,1,"4. Reductions");
  }
  for(;;){
    if (this->bucket_list[1] != NULL) {
      u = pop_node(this,1);

      if (this->dump_file) {
        sprintf(txt,"4.%d RI-Reduction of Node n%d",step++,u);
        dump_section(this->dump_file,2,txt);
        dump_reductiongraph(this,u);
      }

      apply_RI(this,u); 
    } else if (this->bucket_list[2] != NULL) {
      u = pop_node(this,2);

      if (this->dump_file) {
        sprintf(txt,"4.%d RII-Reduction of Node n%d",step++,u);
        dump_section(this->dump_file,2,txt);
        dump_reductiongraph(this,u);
      }

      apply_RII(this,u);
    } else if ((u = pop_maxnode(this)) != -1) {

      if (this->dump_file) {
        sprintf(txt,"4.%d RN-Reduction of Node n%d",step++,u);
        dump_section(this->dump_file,2,txt);
        dump_reductiongraph(this,u);
      }

      apply_RN(this,u);
    } else {
      break;
    }
  } 
}

/*****************************************************************************
 * PBQP back propagation
 ****************************************************************************/

/* determine solution of a reduced node. Either
   RI or RII was applied for this node. */
static
void determine_solution(pbqp *this,int x)
{
  vec *v = v_copy(this -> node_costs[x]);
  adjnode *adj_ptr;

  ASSERT(this != NULL);
  ASSERT(x >= 0 && x < this->num_nodes);
  ASSERT(this -> adj_list != NULL);
  ASSERT(this -> solution != NULL);

  for(adj_ptr=this->adj_list[x] ;adj_ptr != NULL; adj_ptr = adj_ptr -> succ) {
    int y = adj_ptr -> adj;
    int y_sol = this -> solution[y];

    mat *c_yx = get_costmatrix(this,y,x);
    ASSERT(y_sol >= 0 && y_sol < this->node_costs[y]->len);
    v_addrow(v,c_yx,y_sol);
    m_free(c_yx);
  }
  this -> solution[x] = v_minidx(v);

  v_free(v);
}

/* back popagation phase of PBQP */
static
void back_propagate(pbqp *this)
{
   int i;

   ASSERT(this != NULL);
   ASSERT(this->stack != NULL);
   ASSERT(this->stack_ptr < this->num_nodes);

  if (this->dump_file != NULL && this -> stack_ptr > 0) {
    dump_section(this->dump_file,2,"5.3. Back Propagation"); 
    fputs("<ul type=\"square\">\n",this->dump_file);
  }

  for(i=this -> stack_ptr-1;i>=0;i--) {
    int x = this -> stack[i];
    ASSERT( x >= 0 && x < this -> num_nodes);
    reinsert_node(this,x);
    if (this->solution[x] == -1) {
      determine_solution(this,x);
      if (this->dump_file != NULL) {
       fprintf(this->dump_file,"<li>node n%d is set to %d</li>",x,this->solution[x]+1);
      }
    }


  }

  if (this->dump_file != NULL) {
    fputs("</ul>\n",this->dump_file);
  }
}

/* solve trivial nodes of degree zero */
static
void solve_trivialnodes(pbqp *this)
{
  ASSERT( this != NULL);
  ASSERT( this -> bucket_list != NULL);

  if (this->dump_file != NULL) {
     dump_section(this->dump_file,1,"5. Determine Solution/Minimum"); 
     dump_section(this->dump_file,2,"5.1. Trivial Solution"); 
     fputs("<ul type=\"square\">\n",this->dump_file);
  }

  while (this->bucket_list[0] != NULL) {
    int u = pop_node(this,0);
    ASSERT( u >= 0 && u < this -> num_nodes);
    this->solution[u] = v_minidx(this->node_costs[u]);
    this->min += V_ELEM(this->node_costs[u],this->solution[u]);

    if (this->dump_file != NULL) {
       fprintf(this->dump_file,"<li>node n%d is set to %d<br> <tex> \\vec{c}_{%d}=",u,this->solution[u]+1,u);
       v_texprint(this->dump_file,this->node_costs[u]);
       fputs("</tex></li>",this->dump_file);
    }
  }
  if (this->dump_file != NULL) {
     fputs("</ul>\n",this->dump_file);
     dump_section(this->dump_file,2,"5.2. Minimum"); 
     if (this->optimal) {
       fprintf(this->dump_file,"<p>Minimum is equal to %3.2g.</p>\n",this->min);
     } else {
       fprintf(this->dump_file,"<p>Minimum was not achieved. Function value is equal to %3.2g.</p>\n",this->min);
     }
  }
}

/*****************************************************************************
 * debug facilities
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

/* solve PBQP problem */
void h_solve_pbqp(pbqp *this)
{
  
  ASSERT(this != NULL);
  ASSERT(!this->solved); 
  
  /* check vector & matrix dimensions */
  check_pbqp(this);

  /* dump input */
  if (this->dump_file) {
    dump_header(this);
    dump_input(this);
  }

  /* simplify PBQP problem */  
  
  /* eliminate trivial nodes, i.e.
     nodes with cost vectors of length one.
  */
  eliminate_trivial_nodes(this);
  if (this->dump_file){
     dump_trivialnodes(this);
  }

  /* eliminate edges with independent 
     cost matrices and normalize matrices
  */
  eliminate_independent_edges(this);
  
  /* create bucket list for graph parsing */
  create_bucketlist(this);
  
  /* reduce phase */
  reduce_pbqp(this);
  
  /* solve trivial nodes */
  solve_trivialnodes(this);

  /* back propagation phase */
  back_propagate(this); 
  
  this->solved = TRUE;

  /* dump footer */
  if (this->dump_file) {
     dump_footer(this);
  }
}

/*****************************************************************************
 * getter functions
 ****************************************************************************/

/* get solution of a node */
int h_get_solution(pbqp *this,int x)
{
  ASSERT(this != NULL);
  ASSERT(this->solution != NULL);
  ASSERT(this -> solved);
  
  return this->solution[x];
}

/* get minimum of function */
num h_get_min(pbqp *this)
{
  ASSERT(this != NULL);
  ASSERT(this -> solved);
  
  return this->min;
}

/* get number of nodes */
int h_get_numnodes(pbqp *this)
{
  return this->num_nodes;
}

/* is solution optimal? */
boolean h_is_optimal(pbqp *this)
{
  return this->optimal;
}

/* set dump file */
void h_set_dumpfile(pbqp *this,FILE *f)
{
  ASSERT(this != NULL);
  ASSERT(f != NULL);

  this->dump_file = f;
}

/* end of h_pbqp.h */
