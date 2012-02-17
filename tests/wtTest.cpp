#include "TestHarness.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#include "wt.h"

uint8_t* init_T(size_t* n)
{
    size_t i;
    uint8_t* T = (uint8_t*) malloc(50);

    for (i=0; i<50; i++) {
        T[i] = i % 8;
    }
    *n = 50;
    return T;
}

uint8_t* init_TRand(size_t* n)
{
    size_t i;
    uint8_t* T = (uint8_t*) malloc(50000);

    for (i=0; i<50000; i++) {
        T[i] = rand() % 256;
    }
    *n = 50000;
    return T;
}

TEST(wt , saveload)
{
    size_t n,i,j;
    uint8_t* T = init_TRand(&n);
    uint8_t* Tcopy = (uint8_t*) malloc(n);
    memcpy(Tcopy,T,n);

    FILE* f = fopen("wt.test","w");

    wt_t* wt = wt_create((uint64_t*)T,8,n,4);
    wt_save(wt,f);
    fclose(f);

    f = fopen("wt.test","r");
    wt_t* wtl = wt_load(f);
    fclose(f);

    /* select test */
    for (i=0; i<200; i++) {
        size_t pos = rand() % n;
        size_t sym = Tcopy[pos];
        size_t cnt = 0;
        for (j=0; j<=pos; j++) if (Tcopy[j]==sym) cnt++;
        CHECK(wt_select(wt,sym,cnt) == wt_select(wtl,sym,cnt));
    }

    /* access */
    for (i=0; i<n; i++) {
        CHECK(wt_access(wt,i) == Tcopy[i]);
    }

    /* rank */
    for (i=0; i<200; i++) {
        size_t pos = rand() % n;
        for (j=0; j<255; j++) {
            CHECK(wt_rank(wt,j,pos) == wt_rank(wtl,j,pos));
        }
    }

    remove("wt.test");
    wt_free(wt);
    wt_free(wtl);
    free(Tcopy);
}


TEST(wt , select)
{
    size_t n,i,j;
    uint8_t* T = init_TRand(&n);
    uint8_t* Tcopy = (uint8_t*) malloc(n);
    memcpy(Tcopy,T,n);

    wt_t* wt = wt_create((uint64_t*)T,8,n,4);

    for (i=0; i<200; i++) {
        size_t pos = rand() % n;
        size_t sym = Tcopy[pos];
        size_t cnt = 0;
        for (j=0; j<=pos; j++) if (Tcopy[j]==sym) cnt++;
        CHECK(wt_select(wt,sym,cnt) == pos);
    }

    wt_free(wt);
    free(Tcopy);
}

TEST(wt , access)
{
    size_t n,i;
    uint8_t* T = init_T(&n);
    wt_t* wt = wt_create((uint64_t*)T,8,n,4);

    for (i=0; i<n; i++) {
        CHECK(wt_access(wt,i) == (i%8));
    }

    wt_free(wt);

    T = init_TRand(&n);
    uint8_t* Tcopy = (uint8_t*) malloc(n);
    memcpy(Tcopy,T,n);
    wt = wt_create((uint64_t*)T,8,n,4);

    for (i=0; i<n; i++) {
        CHECK(wt_access(wt,i) == Tcopy[i]);
    }

    wt_free(wt);
    free(Tcopy);
}

TEST(wt , rank)
{
    size_t n,i,j;
    uint8_t* T = init_TRand(&n);
    uint8_t* Tcopy = (uint8_t*) malloc(n);
    memcpy(Tcopy,T,n);

    wt_t* wt = wt_create((uint64_t*)T,8,n,4);

    for (i=0; i<200; i++) {
        size_t pos = rand() % n;
        size_t count[256] = {0};
        for (j=0; j<=pos; j++) count[ Tcopy[j] ]++;
        for (j=0; j<255; j++) {
            CHECK(wt_rank(wt,j,pos) == count[j]);
        }
    }

    wt_free(wt);
    free(Tcopy);
}

TEST(wt , setgetsym)
{
    uint64_t* T = (uint64_t*) malloc(8*sizeof(uint64_t));

    for (size_t i=0; i<64; i++) {
        wt_setsym(T,8,i,(i*2)%256);
    }

    for (size_t i=0; i<64; i++) {
        CHECK(wt_getsym(T,8,i) == (i*2)%256);
    }

    free(T);
}

TEST(wt , quantile)
{
    uint64_t* T = (uint64_t*) malloc(4*sizeof(uint64_t));
    wt_setsym(T,6,0,7);
    wt_setsym(T,6,1,2);
    wt_setsym(T,6,2,13);
    wt_setsym(T,6,3,12);
    wt_setsym(T,6,4,9);
    wt_setsym(T,6,5,11);
    wt_setsym(T,6,6,15);
    wt_setsym(T,6,7,14);
    wt_setsym(T,6,8,17);
    wt_setsym(T,6,9,3);
    wt_setsym(T,6,10,8);
    wt_setsym(T,6,11,10);
    wt_setsym(T,6,12,11);
    wt_setsym(T,6,13,1);
    wt_setsym(T,6,14,18);
    wt_setsym(T,6,15,5);
    wt_setsym(T,6,16,6);
    wt_setsym(T,6,17,4);
    wt_setsym(T,6,18,16);

    wt_t* wt = wt_create(T,6,20,4);

    CHECK(wt_quantile(wt,8,13,1) == 1);
    CHECK(wt_quantile(wt,8,13,2) == 3);
    CHECK(wt_quantile(wt,8,13,3) == 8);
    CHECK(wt_quantile(wt,8,13,4) == 10);
    CHECK(wt_quantile(wt,8,13,5) == 11);
    CHECK(wt_quantile(wt,8,13,6) == 17);

    wt_free(wt);
}

TEST(wt , mostfrequent)
{
    uint64_t* T = (uint64_t*) malloc(4*sizeof(uint64_t));
    wt_setsym(T,6,0,7);
    wt_setsym(T,6,1,2);
    wt_setsym(T,6,2,2);
    wt_setsym(T,6,3,4);
    wt_setsym(T,6,4,7);
    wt_setsym(T,6,5,4);
    wt_setsym(T,6,6,6);
    wt_setsym(T,6,7,7);
    wt_setsym(T,6,8,3);
    wt_setsym(T,6,9,6);
    wt_setsym(T,6,10,7);
    wt_setsym(T,6,11,4);
    wt_setsym(T,6,12,6);
    wt_setsym(T,6,13,5);
    wt_setsym(T,6,14,18);
    wt_setsym(T,6,15,7);
    wt_setsym(T,6,16,4);
    wt_setsym(T,6,17,2);
    wt_setsym(T,6,18,1);

    wt_t* wt = wt_create(T,6,20,4);

    wt_result_t* res = wt_mostfrequent(wt,0,18,2);

    CHECK(res->m == 2);
    CHECK(res->items[0].sym == 7);
    CHECK(res->items[1].sym == 4);
    CHECK(res->items[0].freq == 5);
    CHECK(res->items[1].freq == 4);

    wt_freeresult(res);

    wt_free(wt);
}
