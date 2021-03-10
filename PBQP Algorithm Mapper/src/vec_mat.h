/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Header file for vector/matrix objects 
 *****************************************************************************
 $Id: vec_mat.h,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $ 
 *****************************************************************************
 $Log: vec_mat.h,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/

#ifndef __VEC_MAT_H__
#define __VEC_MAT_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "support.h"

/* defines a cost number used in matrices and vecs 
 * of the PBQP problem */
typedef float num;

/* booleanean data type */
typedef enum boolean {TRUE = 1, FALSE = 0} boolean;

/* vec definition */
typedef struct vec { 
  int len;          /* length of vec */
  num data[1];      /* number of the vecs */ 
} vec;


/* mat definition */
typedef struct mat { 
  int rows;         /* number of rows */
  int cols;         /* number of columns */
  num data[1];      /* numbers of the mat */
} mat;	

#define NULL_COSTS (0)
#define INF_COSTS (INFINITY)
#define EPS (MINFLOAT)

/*-----------------*
 * Vector Routines *
 *-----------------*/

/* allocating a vector of size len */
vec *v_alloc(int len);
/* free vector */
void v_free(vec *v);
/* print vector */
void v_print(FILE *f,vec *v);
/* print vector in tex format */
void v_texprint(FILE *f,vec *v);
/* copy vector */
vec *v_copy(vec *v);
/* add vector */
void v_add(vec *dest,vec *src);
/* set vector element */
void v_set(vec *v, int i, num c);
/* get vector element */
num v_get(vec *v, int i);
/* get diagonal of a matrix */
vec *v_diagonalize(mat *m);
/* determine minimum of vector */
num v_min(vec *v);
/* determine smallest index of minimum */
int v_minidx(vec *v);
/* add row of a matrix to vector */
void v_addrow(vec *v,mat *m,int row);
/* get row of a matrix */
vec *v_row(mat *m,int row);
/* read vector */
vec *v_read(FILE *f);


/*-----------------*
 * Matrix Routines *
 *-----------------*/

/* allocate matrix */
mat *m_alloc(int rows,int cols);
/* free matrix */
void m_free(mat *m);
/* copy mat */
mat *m_copy(mat *m);
/* transpose mat */
mat *m_transpose(mat *m);
/* add two matrices */
void m_add(mat *dest,mat *src);
/* add transposed matrix */
void m_addtransposed(mat *dest,mat *src);
/* set matrix element */
void m_set(mat *m,int i,int j,num c);
/* set matrix element */
num m_get(mat *m,int i,int j);
/* subtract number from a row */
void m_subrow(mat *m,int row,num value);
/* set row */
void m_setrow(mat *m,int row,num value);
/* subtract number form a col */
void m_subcol(mat *m,int col,num value);
/* set column */
void m_setcol(mat *m,int col,num value);
/* get minimum of a row */
num m_rowmin(mat *m,int row);
/* get minimum of a column */
num m_colmin(mat *m,int col);
/* determine whether matrix is a zero matrix */
boolean m_iszero(mat *m);
/* reset matrix */
void m_reset(mat *m,num value);
/* read matrix */
mat *m_read(FILE *f);
/* print mat */
void m_print(FILE *f,mat *m);
/* print mat in tex format*/
void m_texprint(FILE *f,mat *m);


#ifndef FAST
#define M_ELEM(m,i,j) m_get(m,i,j)
#define V_ELEM(v,i) v_get(v,i)
#else 
#define M_ELEM(m,i,j) ((m)->data[(i) * (m)->cols + (j)])
#define V_ELEM(v,i) ((v)->data[(i)])
#endif

/* assert macro */
#ifndef ASSERT
#ifdef FAST
#define ASSERT(a)
#else
#define ASSERT(a) assert(a)
#endif
#endif

#endif
