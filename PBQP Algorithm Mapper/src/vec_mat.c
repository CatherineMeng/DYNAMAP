/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Vector/Matrix Routines 
 *****************************************************************************
 $Id: vec_mat.c,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $ 
 *****************************************************************************
 $Log: vec_mat.c,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "vec_mat.h"


/*-----------------*
 * Vector Routines *
 *-----------------*/

/* allocating vector of size "len" */
vec *v_alloc(int len)
{
  vec *v;
  int i;
  
  ASSERT(len > 0);
  v = (vec *)alloc_mem(sizeof(vec)+sizeof(num)*(len-1));
  ASSERT(v!=NULL);

  v->len = len;
  for(i=0;i<len;i++) {
    v->data[i] = NULL_COSTS;
  }
  
  return v;
}


/* de-allocate vector */
void v_free(vec *v)
{
  ASSERT(v != NULL);
  free_mem(v);
}

/* print vector */
void v_print(FILE *f,vec *v)
{
  int i;
  if (v!=NULL) {
    fprintf(f,"%d\n",v->len);
    for(i=0;i<v->len;i++) {
      fprintf(f," %f",v->data[i]);
    }
    fprintf(f,"\n");
  } else {
    fprintf(f,"0\n");
  }
}
/* print vector */
void v_texprint(FILE *f,vec *v)
{
  int i;
  if (v!=NULL) {
    ASSERT(v->len > 0);
    fprintf(f,"\t\\begin{pmatrix}\n");
    fprintf(f,"\t %6.4g",v->data[0]);
    for(i=1;i<v->len;i++) {
      fprintf(f,"& %6.4g",v->data[i]);
    }
    fprintf(f,"\n\t\\end{pmatrix}\n");
  }
}

/* copy vector */
vec *v_copy(vec *v)
{
  vec *dest;
  int len = v->len;
  int i;
  
  ASSERT(v != NULL);
  dest = v_alloc(len);
  for(i=0;i<len;i++) {
    dest->data[i] = v->data[i];
  }
  return dest;
}

/* add vector */
void v_add(vec *dest,vec *src)
{
  int i;
  int len;
  
  ASSERT(dest != NULL);
  ASSERT(src != NULL);
  len = dest->len;
  ASSERT(len > 0);
  ASSERT(src->len == len);
  for(i=0;i<len;i++) {
    dest->data[i] += src->data[i]; 
  }
}

/* set vector element */
void v_set(vec *v, int i, num c)
{
  ASSERT(v != NULL);
  ASSERT(i >= 0 && i < v->len);
  
  v->data[i] = c;
}

/* get vector element */
num v_get(vec *v,int i)
{
  ASSERT(v != NULL);
  ASSERT(i >= 0 && i < v->len);
  
  return v->data[i];
}

/* get diagonal of a matrix */
vec *v_diagonalize(mat *m)
{
   vec *v;
   int i;

   ASSERT(m != NULL);
   ASSERT(m->rows == m-> cols) ;
   ASSERT(m->rows > 0);

   v = v_alloc(m->rows);
   for(i=0;i<m->rows;i++) {
     v->data[i]  = m->data[i*m->cols + i]; 
   }

   return v;
}

/* determine minimum of vector */
num v_min(vec *v)
{
   int i;
   num min;

   ASSERT( v != NULL);
   ASSERT( v->len > 0);

   min = v->data[0];
   for(i=1;i<v->len;i++){
     if (v->data[i] < min) {
        min = v->data[i];
     }
   }
   return min;
}

/* determine minimum of vector */
int v_minidx(vec *v)
{
   int i;
   num min;
   int min_idx = 0;

   ASSERT( v != NULL);
   ASSERT( v->len > 0);

   min = v->data[0];
   for(i=1;i<v->len;i++){
     if (v->data[i] < min) {
        min = v->data[i];
        min_idx = i;
     }
   }
   return min_idx;
}

/* read vector from file */
vec *v_read(FILE *f)
{
   int len;
   int i;
   vec *v;

   assert(fscanf(f,"%d",&len) == 1);
   v = v_alloc(len); 
   for(i=0;i<len;i++) {
     num value;
     assert(fscanf(f,"%f",&value) == 1);
     v->data[i] = value; 
   }

   return v;
}


/* add to vector a row of a matrix */
void v_addrow(vec *v,mat *m,int row)
{
   int col;

   ASSERT(v != NULL);
   ASSERT(m != NULL);
   ASSERT(row >= 0 && row < m->rows);
   ASSERT(v->len == m->cols);

   for(col = 0; col < v->len ; col ++) {
     v->data[col] += m->data[row*m->cols + col];
   }
}

/* get row of matrix */
vec *v_row(mat *m,int row)
{
   int col;
   vec *v;

   ASSERT( m != NULL);
   ASSERT( row >= 0 && row < m->rows);
   ASSERT( m->cols > 0 );
  
   v = v_alloc(m->cols);
   for(col=0;col < m->cols; col ++) {
      v->data[col] = m->data[row*m->cols + col];  
   }
   return v;
}
 
/*-----------------*
 * Matrix Routines *
 *-----------------*/

/* allocate matrix */
mat *m_alloc(int rows,int cols)
{
  int i;
  int size = rows*cols;
  mat *p;
 
  ASSERT(rows > 0);
  ASSERT(cols > 0);
  p = (mat *)alloc_mem(sizeof(mat)+sizeof(num)*(size-1));
  ASSERT(p != NULL);

  p->rows = rows;
  p->cols = cols;
  for(i=0;i<size;i++) {
    p->data[i] = NULL_COSTS;
  }
  return p;
}

/* de-allocate matrix */
void m_free(mat *m)
{
  ASSERT(m!=NULL);
  free_mem(m); 
}

/* copy mat */
mat *m_copy(mat *m)
{
  mat *dest;
  int size;
  int i;
  
  ASSERT(m != NULL);
  size = m->rows * m->cols;

  ASSERT(size > 0);
  dest = m_alloc(m->rows,m->cols);

  for(i=0;i<size;i++){
    dest->data[i] = m->data[i];
  }

  return dest;
}

/* transpose mat */
mat *m_transpose(mat *m)
{
  int i,j;
  mat *dest;
  int rows = m->rows;
  int cols = m->cols;

  dest = m_alloc(cols,rows);

  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++) {
       dest->data[j*rows+i] = m->data[i*cols+j];
    }
  }
  return dest;
}

/* add two matrices */
void m_add(mat *dest,mat *src)
{
  int i;
  int size;

  ASSERT(src != NULL);
  ASSERT(dest != NULL);
   
  ASSERT(src->rows == dest->rows);
  ASSERT(src->cols == dest->cols);

  size = dest->rows * dest->cols;

  for(i=0;i<size;i++) {
    dest->data[i] += src->data[i];
  }
}

/* add transposed matrix */
void m_addtransposed(mat *dest,mat *src)
{
  mat *m = m_transpose(src);
  m_add(dest,m);
  m_free(m);
}

/* set a matrix element */
void m_set(mat *m,int i,int j,num c)
{
  ASSERT(m != NULL);
  ASSERT(i >= 0 && i < m->rows);
  ASSERT(j >= 0 && j < m->cols);

  m->data[i*m->cols+j] = c;
}

/* get  matrix element */
num m_get(mat *m,int i,int j)
{
  ASSERT(m != NULL);
  ASSERT(i >= 0 && i < m->rows);
  ASSERT(j >= 0 && j < m->cols);

  return m->data[i*m->cols+j];
}

/* get minimum of a row */
num m_rowmin(mat *m,int row)
{
    int col;
    num min;

    ASSERT(m != NULL);
    ASSERT(row >= 0 && row < m->rows);
   
    min=M_ELEM(m,row,0);
    for(col=1;col<m->cols;col++) {
      num v = M_ELEM(m,row,col);
      if ( v < min) {
        min = v;
      }
    }
    return min; 
}

/* get minimum of a column */
num m_colmin(mat *m,int col)
{
    int row;
    num min;

    ASSERT(m != NULL);
    ASSERT(col >= 0 && col < m->cols);
   
    min=M_ELEM(m,0,col);
    for(row=1;row<m->rows;row++) {
      num v = M_ELEM(m,row,col);
      if ( v < min) {
        min = v;
      }
    }
    return min; 
}

/* set row */
void m_setrow(mat *m,int row,num value)
{
   int col;

   ASSERT( m != NULL);
   ASSERT( row >= 0 && row < m->rows);
   for(col=0;col< m->cols;col++) {
     m->data[row * m->cols + col] = value;
   }
}

/* subtract number from a column */
void m_subrow(mat *m,int row,num value)
{
   int col;

   ASSERT( m != NULL);
   ASSERT( row >= 0 && row < m->rows);
   for(col=0;col< m->cols;col++) {
     m->data[row * m->cols + col] -= value;
   }
}

/* set column */ 
void m_setcol(mat *m,int col,num value)
{
   int row;

   ASSERT( m != NULL);
   ASSERT( col >= 0 && col < m->cols);

   for(row=0;row< m->rows;row++) {
     m->data[row * m->cols + col] = value;
   }
}
 
/* subtract number from a row */
void m_subcol(mat *m,int col,num value)
{
   int row;

   ASSERT( m != NULL);
   ASSERT( col >= 0 && col < m->cols);

   for(row=0;row< m->rows;row++) {
     m->data[row * m->cols + col] -= value;
   }
}

/* determine whether a matrix is a zero matrix */
void m_reset(mat *m,num value)
{
   int i;
   int size;

   ASSERT(m != NULL);
   size = m->cols * m->rows; 

   for(i=0;i<size;i++){
     m -> data[i] = value;
   }
} 

/* determine whether a matrix is a zero matrix */
boolean m_iszero(mat *m)
{
   int i;
   int size;

   ASSERT(m != NULL);
   size = m->cols * m->rows; 

   for(i=0;i<size;i++)
     if (fabs(m -> data[i] - NULL_COSTS) > EPS)  
        return FALSE;

   return TRUE;
} 

/* Matrix I/O */

/* read matrix */
mat *m_read(FILE *f)
{
   int rows,
       cols;
   int i,j;
   mat *m;

   assert(fscanf(f,"%d %d",&rows,&cols) == 2);

   m = m_alloc(rows,cols);
   for(i=0;i<rows;i++) {
     for(j=0;j<cols;j++) { 
        num value;
        assert(fscanf(f,"%f",&value) == 1);
        m->data[i*cols + j] = value; 
     }
   }

   return m;
}

/* print mat */
void m_print(FILE *f,mat *m)
{
  int i,j;
  num *p = m->data;

  if  (m!=NULL) {
    fprintf(f,"%d %d\n",m->rows,m->cols);
    for(i=0;i<m->rows;i++) {
      for(j=0;j<m->cols;j++) {
	fprintf(f," %f",*p++);
      }
      fprintf(f,"\n");
    }
  } else {
    fprintf(f,"0 0\n");
  }
}


/* print mat */
void m_texprint(FILE *f,mat *m)
{
  int i,j;
  num *p = m->data;

  if  (m!=NULL) {
    ASSERT(m->cols > 0);
    ASSERT(m->rows > 0);
    fprintf(f,"\t\\begin{pmatrix}\n");
    for(i=0;i<m->rows;i++) {
      fprintf(f,"\t %6.4g",*p++);
      for(j=1;j<m->cols;j++) {
	fprintf(f,"& %6.4g",*p++);
      }
      fprintf(f,"\\\\\n");
    }
    fprintf(f,"\t\\end{pmatrix}\n");
  }
}



