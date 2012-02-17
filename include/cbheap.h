
#ifndef CBHEAP_H
#define CBHEAP_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define PARENT(x) ((x-1) >> 1)
#define LEFTCHILD(x) ((x<<1)+1)
#define RIGHTCHILD(x) ((x<<1)+2)

#define INIT_CBHEAP_SIZE		1024

typedef struct cbheap {
    uint64_t                    n;
    uint64_t                    size;
    int (*itemcmp)(const void*,const void*);
    void (*freeitem)(void*);
    void**                      A;
} cbheap_t;

cbheap_t* cbheap_create(int (*cmpfunc)(const void*,const void*),void (*dfnc)(void*));
void      cbheap_free(cbheap_t* cbh);
void*     cbheap_top(cbheap_t*);
void      cbheap_insert(cbheap_t*,void* item);
void      cbheap_delete_top(cbheap_t*);
size_t    cbheap_count(cbheap_t* ch);
size_t    cbheap_spaceusage(cbheap_t* ch);
void      cbheap_swap(void** a,void** b);
uint32_t  cbheap_height(cbheap_t* ch);
void      cbheap_resize(cbheap_t* ch);
void      cbheap_siftup(cbheap_t* ch,size_t start,size_t end);
void      cbheap_siftdown(cbheap_t* ch,size_t start,size_t end);
int       cbheap_empty(cbheap_t* ch);

/* misc */
void*     cbheap_safecalloc(size_t n);
void*     cbheap_saferealloc(void*,size_t n);


#endif

