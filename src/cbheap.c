
#include "cbheap.h"


cbheap_t*
cbheap_create(int (*cmpfunc)(const void*,const void*),void (*del)(void*))
{
    cbheap_t* cbh = (cbheap_t*) cbheap_safecalloc(sizeof(cbheap_t));
    cbh->n = 0;
    cbh->itemcmp = cmpfunc;
    cbh->freeitem = del;
    cbh->A = (void**) cbheap_safecalloc(INIT_CBHEAP_SIZE*sizeof(void*));
    cbh->size = INIT_CBHEAP_SIZE;

    return cbh;
}

void
cbheap_free(cbheap_t* cbh)
{
    size_t i;
    for (i=0; i<cbh->n; i++) cbh->freeitem(cbh->A[i]);
    free(cbh->A);
    free(cbh);
}

void*
cbheap_top(cbheap_t* ch)
{
    if (ch->n) return ch->A[0];
    return NULL;
}

void
cbheap_insert(cbheap_t* ch,void* item)
{
    if (ch->size == ch->n) cbheap_resize(ch);

    /** add item to the end of array */
    ch->A[ch->n] = item;
    ch->n++;

    /* heapify */
    cbheap_siftup(ch,0,ch->n-1);
}

void
cbheap_delete_top(cbheap_t* ch)
{
    void* del = NULL;
    /* empty heap */
    if (ch && ch->n == 0) return;

    del = ch->A[0];

    /* put last item in top pos */
    ch->A[0] = ch->A[ch->n-1];

    /* make heap */
    ch->n--;
    cbheap_siftdown(ch,0,ch->n);

    ch->freeitem(del);
}

void
cbheap_siftup(cbheap_t* ch,size_t start,size_t end)
{
    size_t child = end;
    size_t parent;

    while (child > start) {
        parent = PARENT(child);
        if (ch->itemcmp(ch->A[parent],ch->A[child]) > 0) {
            cbheap_swap(&(ch->A[parent]),&(ch->A[child]));
            child = parent;
        } else {
            return;
        }
    }
}

void
cbheap_siftdown(cbheap_t* ch,size_t start,size_t end)
{
    size_t swap,left,right;
    size_t root = start;
    while (LEFTCHILD(root) < end) {
        swap = root;
        left = LEFTCHILD(root);
        right = RIGHTCHILD(root);
        if (ch->itemcmp(ch->A[swap],ch->A[left]) > 0)
            swap = left;
        if (right < end && ch->itemcmp(ch->A[swap],ch->A[right]) > 0)
            swap = right;
        if (swap != root) {
            cbheap_swap(&(ch->A[swap]),&(ch->A[root]));
            root = swap;
        } else
            return;
    }
}

size_t
cbheap_count(cbheap_t* ch)
{
    return ch->n;
}

int
cbheap_empty(cbheap_t* ch)
{
    if (ch->n==0) return 1;
    else return 0;
}

size_t
cbheap_spaceusage(cbheap_t* ch)
{
    return sizeof(cbheap_t) + ch->n*sizeof(void*);
}

void
cbheap_swap(void** a,void** b)
{
    void* c = *a;
    *a = *b;
    *b = c;
}

uint32_t
cbheap_height(cbheap_t* ch)
{
    size_t n = ch->n;
    size_t height = 0;
    while (n) {
        height++;
        n = n>>1;
    }
    return height;
}

void
cbheap_resize(cbheap_t* ch)
{
    size_t new_size = 2*ch->size;
    ch->A = (void**) cbheap_saferealloc(ch->A,new_size*sizeof(void*));
    ch->size = new_size;
}

/* misc */
void*
cbheap_safecalloc(size_t n)
{
    void* mem = calloc(n,1);
    if (!mem) {
        fprintf(stderr,"ERROR: cbheap_safecalloc()");
        exit(EXIT_FAILURE);
    }
    return mem;
}

void*
cbheap_saferealloc(void* old_mem,size_t new_size)
{
    if ((old_mem = realloc(old_mem, new_size)) == NULL) {
        fprintf(stderr, "ERROR: cbheap_saferealloc() cannot allocate"
                "%zu blocks of memory.\n", new_size);
        exit(EXIT_FAILURE);
    }
    return (old_mem);
}
