
#include "rankbv.h"
#include "wt.h"
#include "cbheap.h"

/*#define _WT_DEBUG_*/
#include <time.h>

wt_t*
wt_init(size_t n)
{
    wt_t* wt = (wt_t*) wt_safecalloc(sizeof(wt_t));

    wt->n = n;
    wt->height = 0;
    wt->occ = NULL;
    wt->bittree = NULL;
    wt->max_v = 0;

    return wt;
}

wt_t*
wt_create(uint64_t* A,size_t bits,size_t n,uint32_t f)
{
    size_t i;
    wt_t* wt = wt_init(n);

    /* calc height */
    for (i=0; i<n; i++) wt->max_v = wt_max(wt_getsym(A,bits,i),wt->max_v);
    wt->height = wt_bits(wt->max_v);

#ifdef _WT_DEBUG_
    fprintf(stdout,"wt::create() height = %u max_v %u\n",wt->height,wt->max_v);
#endif

    /* count occs */
    uint64_t* occs = (uint64_t*) wt_safecalloc((wt->max_v+1)*sizeof(uint64_t));
    for (i=0; i<n; i++) occs[wt_getsym(A,bits,i)]++;
    /* cummulative counts */
    for (i=1; i<=wt->max_v; i++) occs[i] += occs[i-1];

    /* build select() struct */
    wt->occ = rankbv_init(n+1,f);
    for (i=0; i<=wt->max_v; i++) rankbv_setbit(wt->occ,occs[i]-1);
    rankbv_setbit(wt->occ,n);
    rankbv_build(wt->occ);
    free(occs);

    /* create tree */
    wt->bittree = (rankbv_t**) wt_safecalloc(wt->height*sizeof(rankbv_t*));
    for (i=0; i<wt->height; i++) wt->bittree[i] = rankbv_init(n,f);

    wt_build(wt,A,bits,n);

    return wt;
}

void
wt_free(wt_t* wt)
{
    size_t i;
    if (wt) {
        if (wt->occ) rankbv_free(wt->occ);
        if (wt->bittree) {
            for (i=0; i<wt->height; i++) rankbv_free(wt->bittree[i]);
            free(wt->bittree);
        }
        free(wt);
    }
}

void
wt_build(wt_t* wt,uint64_t* A,size_t bits,size_t n)
{
    size_t i;

    /* build levels */
    wt_buildlvl(wt,A,bits,0,n,0);

    /* build rank */
#ifdef _WT_DEBUG_
    fprintf(stdout,"\nWT::BuildRank(%zu)\n",i);
#endif
    for (i=0; i<wt->height; i++) {
        rankbv_build(wt->bittree[i]);
    }
}

void
wt_buildlvl(wt_t* wt,uint64_t* A,size_t bits,uint32_t lvl,size_t n,size_t offset)
{
    size_t i;
    if (lvl == wt->height) {
        free(A);
        return;
    }

    /* count symbols */
    size_t cone = 0;
    for (i=0; i<n; i++) {
        if (wt_marked(wt_getsym(A,bits,i) ,wt->height,lvl)) cone++;
    }
    size_t czero = n - cone;

    /* move symbols to left/right subarray */
    size_t cleft,cright;
    cleft = cright = 0;

    uint64_t* left = (uint64_t*)  wt_safecalloc(((czero*bits)/RBVW+1)*sizeof(uint64_t));
    uint64_t* right = (uint64_t*) wt_safecalloc(((cone*bits)/RBVW+1)*sizeof(uint64_t));

    for (i=0; i<n; i++) {
        uint32_t sym = wt_getsym(A,bits,i);
        if (wt_marked(sym,wt->height,lvl)) {
            wt_setsym(right,bits,cright++,sym);
            rankbv_setbit(wt->bittree[lvl],offset+i);
        } else {
            wt_setsym(left,bits,cleft++,sym);
        }
    }
    free(A);

    wt_buildlvl(wt,left,bits,lvl+1,cleft,offset);
    wt_buildlvl(wt,right,bits,lvl+1,cright,offset+cleft);
}

size_t
wt_count(wt_t* wt,uint32_t sym)
{
    if (sym) return rankbv_select1(wt->occ,sym) - rankbv_select1(wt->occ,sym-1) + 1;
    return rankbv_select1(wt->occ,sym);
}


uint32_t
wt_access(wt_t* wt,size_t pos)
{
    uint32_t lvl = 0;
    uint32_t ret = 0;
    size_t start = 0;
    size_t end = wt->n - 1;
    size_t before;
    rankbv_t* bs;

    while (lvl<wt->height) {
        bs = wt->bittree[lvl];

        if (start==0) before = 0;
        else before = rankbv_rank1(bs,start-1);

        if (rankbv_access(bs,pos)) {
            ret = wt_mark(ret,wt->height,lvl);
            pos = rankbv_rank1(bs,pos-1)-before;
            start = rankbv_rank1(bs,end) - before;
            start = end-start+1;
            pos += start;
        } else {
            pos = pos - start - (rankbv_rank1(bs,pos) - before);
            end = end - start - (rankbv_rank1(bs,end) - before);
            end += start;
            pos += start;
        }
        lvl++;
    }
    return ret;
}

size_t
wt_select(wt_t* wt,uint32_t sym,size_t j)
{
    uint32_t mask = (1<<wt->height)-2;
    uint32_t sum = 2;
    uint32_t lvl = wt->height-1;
    size_t pos=j;
    size_t start = 0;
    size_t end = 0;
    rankbv_t* bs;

    while (1) {
        bs = wt->bittree[lvl];
        start = sym&mask; /* get start */
        end = wt_min(wt->max_v+1,start+sum);

        if (start!=0) start = rankbv_select1(wt->occ,start)+1;
        end = rankbv_select1(wt->occ,end+1)-1;

        size_t ones_start = rankbv_rank1(bs,start-1);
        if (wt_marked(sym,wt->height,lvl))
            pos = rankbv_select1(bs,ones_start+pos)-start+1;
        else
            pos = rankbv_select0(bs,start-ones_start+pos)-start+1;

        mask <<=1;
        sum <<=1;
        if (!lvl) break;
        lvl--;
    }

    return pos-1;
}

size_t
wt_rank(wt_t* wt,uint32_t sym,size_t pos)
{
    uint32_t lvl = 0;
    size_t start = 0;
    size_t end = wt->n-1;
    size_t count=0;
    size_t before;
    rankbv_t* bs;

    while (lvl<wt->height) {
        bs = wt->bittree[lvl];

        if (start==0) before = 0;
        else before = rankbv_rank1(bs,start-1);

        if (wt_marked(sym,wt->height,lvl)) {
            pos=rankbv_rank1(bs,pos)-before-1;
            count=pos+1;
            start=(rankbv_rank1(bs,end)-before);
            start=end-start+1;
            pos+=start;
        } else {
            pos=pos-start+before-rankbv_rank1(bs,pos);
            count=pos+1;
            end=end-start-(rankbv_rank1(bs,end)-rankbv_rank1(bs,start-1));
            end+=start;
            pos+=start;
        }
        lvl++;
        if (count==0) return 0;
    }
    return count;
}

void
wt_print(wt_t* wt)
{
    size_t i;
    if (wt) {
        for (i=0; i<wt->height; i++) {
            fprintf(stdout,"(%zu) ",i);
            rankbv_print(wt->bittree[i]);
        }
    }
}

size_t
wt_spaceusage(wt_t* wt)
{
    size_t i=0;
    size_t treespace = 0;
    for (i=0; i<wt->height; i++) {
        treespace += rankbv_spaceusage(wt->bittree[i]);
    }
    return sizeof(wt) +
           rankbv_spaceusage(wt->occ) +
           treespace;
}

wt_t*
wt_load(FILE* f)
{
    size_t i;
    wt_t* wtl = wt_init(0);
    if (fread(&wtl->n,sizeof(uint64_t),1,f)!=1) {
        fprintf(stdout,"error reading wt->n\n");
        exit(EXIT_FAILURE);
    }
    if (fread(&wtl->height,sizeof(uint32_t),1,f)!=1) {
        fprintf(stdout,"error reading wt->height\n");
        exit(EXIT_FAILURE);
    }
    if (fread(&wtl->max_v,sizeof(uint32_t),1,f)!=1) {
        fprintf(stdout,"error reading wt->max_v\n");
        exit(EXIT_FAILURE);
    }

#ifdef _WT_DEBUG_
    fprintf(stdout,"WT::Load() n=%zu height=%u max_v=%u\n",wtl->n,wtl->height,wtl->max_v);
#endif

#ifdef _WT_DEBUG_
    fprintf(stdout,"WT::Load() rankbv_load(wt->occ)\n");
#endif
    wtl->occ = rankbv_load(f);
    wtl->bittree = (rankbv_t**) wt_safecalloc(wtl->height*sizeof(rankbv_t*));
    for (i=0; i<wtl->height; i++) {
#ifdef _WT_DEBUG_
        fprintf(stdout,"WT::Load() rankbv_load(wt->bittree[%zu])\n",i);
#endif
        wtl->bittree[i] = rankbv_load(f);
    }
    return wtl;
}

void
wt_save(wt_t* wt,FILE* f)
{
#ifdef _WT_DEBUG_
    fprintf(stdout,"WT::Write() n=%zu height=%u max_v=%u\n",wt->n,wt->height,wt->max_v);
#endif
    size_t i;
    if (fwrite(&(wt->n),sizeof(uint64_t),1,f)!=1) {
        fprintf(stdout,"error writing wt->n\n");
        exit(EXIT_FAILURE);
    }
    uint32_t height = wt->height;
    if (fwrite(&height,sizeof(uint32_t),1,f)!=1) {
        fprintf(stdout,"error writing wt->height\n");
        exit(EXIT_FAILURE);
    }
    if (fwrite(&(wt->max_v),sizeof(uint32_t),1,f)!=1) {
        fprintf(stdout,"error writing wt->max_v\n");
        exit(EXIT_FAILURE);
    }

#ifdef _WT_DEBUG_
    fprintf(stdout,"WT::Write() rankbv_save(wt->occ)\n");
#endif
    rankbv_save(wt->occ,f);
    for (i=0; i<wt->height; i++) {
#ifdef _WT_DEBUG_
        fprintf(stdout,"WT::Write() rankbv_save(wt->bittree[%zu])\n",i);
#endif
        rankbv_save(wt->bittree[i],f);
    }
}

wt_quant_t
wt_quantile_freq(wt_t* wt,size_t left,size_t right,size_t q)
{
    /* decrease q as the smallest element q=1 is
     * found by searching for 0 */
    q--;
    uint32_t sym = 0;
    size_t freq = 0;
    uint32_t lvl = 0;
    size_t start = 0, end = wt->n-1;
    size_t before;
    rankbv_t* bs;

    while (lvl<wt->height) {
        bs = wt->bittree[lvl];

        /* calc start of level bound */
        if (start == 0) before = 0;
        else before = rankbv_rank1(bs,start-1);

        /* number of 1s before T[l..r] */
        size_t rank_before_left = rankbv_rank1(bs,start+left-1);
        /* number of 1s before T[r] */
        size_t rank_before_right = rankbv_rank1(bs,start+right);
        /* number of 1s in T[l..r] */
        size_t num_ones = rank_before_right - rank_before_left;
        /* number of 0s in T[l..r] */
        size_t num_zeros = (right-left+1) - num_ones;

        /* if there are more than q 0s we go right. left otherwise */
        if (q >= num_zeros) { /* go right */
            freq = num_ones; /* calc freq */
            /* set bit to 1 in sym */
            sym = wt_mark(sym,wt->height,lvl);
            /* number of 1s before T[l..r] within the current node */
            left = rank_before_left - before;
            /* number of 1s in T[l..r] */
            right = rank_before_right - before - 1;
            q = q - num_zeros;
            /* calc starting pos of right childnode */
            start = end - (rankbv_rank1(bs,end)-before) + 1;
        } else { /* go left q = q // sym == sym */
            freq = num_zeros; /* calc freq */
            /* number of zeros before T[l..r] within the current node */
            left = left - (rank_before_left - before);
            /* number of zeros in T[l..r] + left bound */
            right = right - (rank_before_right - before);
            /* calc end pos of left childnode */
            end = end - (rankbv_rank1(bs,end) - before);
        }
        lvl++;
    }
    wt_quant_t qf;
    qf.sym = sym;
    qf.freq = freq;
    return qf;
}

uint32_t
wt_quantile(wt_t* wt,size_t left,size_t right,size_t quantile)
{
    wt_quant_t q = wt_quantile_freq(wt,left,right,quantile);
    return q.sym;
}

int topkrange_cmp(const void* a,const void* b)
{
    wt_topkrange_t* wa = (wt_topkrange_t*)a;
    wt_topkrange_t* wb = (wt_topkrange_t*)b;
    return (wb->right-wb->left+1) - (wa->right-wa->left+1);
}

wt_topkrange_t*
new_topk_range(size_t l,size_t r,size_t s,size_t e,uint32_t lvl,uint32_t sym,size_t f)
{
    wt_topkrange_t* res = (wt_topkrange_t*) wt_safecalloc(sizeof(wt_topkrange_t));
    res->left = l;
    res->right = r;
    res->start = s;
    res->end = e;
    res->lvl = lvl;
    res->sym = sym;
    res->freq = f;
    return res;
}

wt_result_t*
wt_mostfrequent(wt_t* wt,size_t left,size_t right,size_t k)
{
    size_t before,newleft,newright;

    wt_result_t* res = wt_newresult();
    cbheap_t* h = cbheap_create(topkrange_cmp,free);
    cbheap_insert(h,new_topk_range(left,right,0,wt->n-1,0,0,0));

    while (h->n > 0) {
        wt_topkrange_t* cr = (wt_topkrange_t*) cbheap_top(h);
        size_t rleft = cr->left;
        size_t rright = cr->right;
        size_t rstart = cr->start;
        size_t rend = cr->end;
        uint32_t rlvl = cr->lvl;
        uint32_t rsym = cr->sym;
        size_t rfreq = cr->freq;

        /* remove old range */
        cbheap_delete_top(h); /* deletes cr */

        rankbv_t* bs = wt->bittree[cr->lvl];
        if (cr->lvl == wt->height) {  /* leaf node */
            wt_addresult(res,rsym,rfreq,0);
            if (res->m == k) {  /* we got the top-k */
                break;
            }
            continue;
        }
        /* calc start of level bound */
        if (cr->start == 0) before = 0;
        else before = rankbv_rank1(bs,rstart-1);

        /* number of 1s before T[l..r] */
        size_t rank_before_left = rankbv_rank1(bs,rstart+rleft-1);
        /* number of 1s before T[r] */
        size_t rank_before_right = rankbv_rank1(bs,rstart+rright);
        /* number of 1s in T[l..r] */
        size_t num_ones = rank_before_right - rank_before_left;
        /* number of 0s in T[l..r] */
        size_t num_zeros = (rright-rleft+1) - num_ones;
        /* new interval bounds */
        size_t newintbound = rend - (rankbv_rank1(bs,rend)-before);

        if (num_ones) {  /* right child */
            size_t sym = wt_mark(rsym,wt->height,rlvl);
            /* number of 1s before T[l..r] within the current node */
            newleft = rank_before_left - before;
            /* number of 1s in T[l..r] */
            newright = rank_before_right - before - 1;
            /* calc starting pos of right childnode */
            size_t start = newintbound + 1;
            cbheap_insert(h,new_topk_range(newleft,newright,start,rend,rlvl+1,sym,num_ones));
        }
        if (num_zeros) {  /* left child */
            /* number of zeros before T[l..r] within the current node */
            newleft = rleft - (rank_before_left - before);
            /* number of zeros in T[l..r] + left bound */
            newright = rright - (rank_before_right - before);
            /* calc end pos of left childnode */
            size_t end = newintbound;
            cbheap_insert(h,new_topk_range(newleft,newright,rstart,end,rlvl+1,rsym,num_zeros));
        }

    }
    cbheap_free(h);
    return res;
}

