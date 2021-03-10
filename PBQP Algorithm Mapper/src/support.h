/*****************************************************************************
 (c) 2003 by B. Scholz
 *****************************************************************************
 Header File for support routines
 *****************************************************************************
 $Id: support.h,v 1.1.1.1 2003/12/28 20:58:54 scholz Exp $
 *****************************************************************************
 $Log: support.h,v $
 Revision 1.1.1.1  2003/12/28 20:58:54  scholz
 Imported sources

 *****************************************************************************/
#ifndef __SUPPORT_H__
#define __SUPPORT_H__

#include <stdlib.h>
#include <assert.h>

/* allocating a memory block 

   note: following function can be simply 
         replaced by something else. However,
         it should guarantee that a valid
         memory block of specified size  
         is returned.
*/

static void *alloc_mem(size_t s){
  void *p;

  p = malloc(s);
  assert(p!=NULL);

  return p;   
}

/* de-allocate memory 

   note: at the moment it is the 
   standard routine. 
*/

#define free_mem free

#endif

