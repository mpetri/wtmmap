
#ifndef WT_H
#define WT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "rankbv.h"

    typedef struct wt {
        uint64_t n;
        uint32_t height;
        uint32_t max_v;
        rankbv_t*  occ;
        rankbv_t** bittree;
    } wt_t;


    /* helper ops from libcds */
    static inline void*
    wt_safecalloc(size_t n)
    {
        void* mem = calloc(n,1);
        if (!mem) {
            fprintf(stderr,"ERROR: rankbv_safecalloc()");
            exit(EXIT_FAILURE);
        }
        return mem;
    }

    static inline void*
    wt_saferealloc(void* old_mem,size_t new_size)
    {
        if ((old_mem = realloc(old_mem, new_size)) == NULL) {
            fprintf(stderr, "ERROR: cbheap_saferealloc() cannot allocate"
                    "%zu blocks of memory.\n", new_size);
            exit(EXIT_FAILURE);
        }
        return (old_mem);
    }

    static inline size_t
    wt_length(wt_t* wt)
    {
        return wt->n;
    }

    static inline uint32_t
    wt_bits(size_t n)
    {
        uint32_t b = 0;
        while (n) {
            b++;
            n >>= 1;
        }
        return b;
    }

#define wt_max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define wt_min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

    static inline int
    wt_marked(uint32_t v,uint32_t h,uint32_t bit)
    {
        return v & (1<<(h-bit-1));
    }

    static inline uint32_t
    wt_mark(uint32_t v,uint32_t h,uint32_t bit)
    {
        return v | (1<<(h-bit-1));
    }


    static inline uint32_t
    wt_getsym(uint64_t* B,register size_t bits,register size_t pos)
    {
        uint32_t* A = (uint32_t*)B;
        register size_t i=pos*bits/RBVW32, j=pos*bits-RBVW32*i;
        uint32_t result;
        if (j+bits <= RBVW32)
            result = (A[i] << (RBVW32-j-bits)) >> (RBVW32-bits);
        else {
            result = A[i] >> j;
            result = result | (A[i+1] << (RBVW-j-bits)) >> (RBVW32-bits);
        }
        return result;
    }

    static inline void
    wt_setsym(uint64_t* B,register size_t bits,register size_t pos,register uint32_t sym)
    {
        uint32_t* A = (uint32_t*)B;
        size_t i=pos*bits/RBVW32, j=pos*bits-i*RBVW32;
        uint32_t mask = ((j+bits) < RBVW32 ? ~0u << (j+bits) : 0)
                        | ((RBVW32-j) < RBVW32 ? ~0u >> (RBVW32-j) : 0);
        A[i] = (A[i] & mask) | sym << j;
        if (j+bits>RBVW32) {
            mask = ((~0u) << (bits+j-RBVW32));
            A[i+1] = (A[i+1] & mask)| sym >> (RBVW32-j);
        }
    }

    /* result structs */

    typedef struct wt_quant {
        uint32_t sym;
        size_t freq;
    } wt_quant_t;

    typedef struct wt_range {
        size_t sp;
        size_t ep;
    } wt_range_t;

    typedef struct wt_item {
        uint32_t sym;
        size_t freq;
        size_t weight;
    } wt_item_t;

    typedef struct wt_result {
        size_t size;
        size_t m;
        wt_item_t* items;
    } wt_result_t;

    typedef struct wt_topkrange {
        size_t left;
        size_t right;
        size_t start;
        size_t end;
        uint32_t sym;
        size_t freq;
        uint32_t lvl;
    } wt_topkrange_t;

    static int
    wt_item_cmp(const void* a,const void* b)
    {
        wt_item_t* wa = (wt_item_t*)a;
        wt_item_t* wb = (wt_item_t*)b;
        float wfa = wa->freq/(float)(wa->weight+1);
        float wfb = wb->freq/(float)(wb->weight+1);
        if (wfa < wfb) return -1;
        if (wfb > wfa) return 1;
        return 0;
    }

    static inline wt_item_t*
    wt_new_item(uint32_t sym,size_t freq,size_t weight)
    {
        wt_item_t* item = (wt_item_t*) wt_safecalloc(sizeof(wt_item_t));
        item->sym = sym;
        item->freq = freq;
        item->weight = weight;
        return item;
    }

    static inline wt_result_t*
    wt_newresult()
    {
        wt_result_t* res = (wt_result_t*) wt_safecalloc(sizeof(wt_result_t));
        res->size = 8192;
        res->items = (wt_item_t*) wt_safecalloc(res->size*sizeof(wt_item_t));
        res->m = 0;

        return res;
    }

    static inline void
    wt_addresult(wt_result_t* res,uint32_t sym,size_t freq,size_t weight)
    {
        if (res->m == res->size) {
            res->items = (wt_item_t*) wt_saferealloc(res->items,2*res->size*sizeof(wt_item_t));
            res->size *= 2;
        }
        res->items[res->m].sym = sym;
        res->items[res->m].freq = freq;
        res->items[res->m].weight = weight;
        res->m++;
    }

    static inline void
    wt_freeresult(wt_result_t* res)
    {
        free(res->items);
        free(res);
    }

    static inline void
    wt_sortresult(wt_result_t* res)
    {
        qsort(res->items,res->m,sizeof(wt_item_t),wt_item_cmp);
    }


    /* rankbv functions */
    wt_t*        wt_init(size_t n);
    wt_t*        wt_create(uint64_t* A,size_t bits,size_t n,uint32_t f);
    void         wt_free(wt_t* wt);
    void         wt_build(wt_t* wt,uint64_t* A,size_t bits,size_t n);
    void         wt_buildlvl(wt_t* wt,uint64_t* A,size_t bits,uint32_t lvl,size_t n,size_t offset);
    uint32_t     wt_access(wt_t* wt,size_t i);
    size_t       wt_rank(wt_t* wt,uint32_t sym,size_t i);
    size_t       wt_select(wt_t* wt,uint32_t sym,size_t x);
    size_t       wt_count(wt_t* wt,uint32_t sym);
    void         wt_print(wt_t* wt);
    uint32_t     wt_quantile(wt_t* wt,size_t left,size_t right,size_t quantile);
    wt_quant_t   wt_quantile_freq(wt_t* wt,size_t left,size_t right,size_t quantile);
    wt_result_t* wt_mostfrequent(wt_t* wt,size_t left,size_t right,size_t k);
    wt_result_t* wt_intersect(wt_t* wt,wt_range_t* ranges,size_t m,size_t threshold);

    /* save/load */
    size_t    wt_spaceusage(wt_t* rbv);
    wt_t*     wt_load(FILE* f);
    void      wt_save(wt_t* rbv,FILE* f);


#ifdef __cplusplus
}
#endif

#endif

