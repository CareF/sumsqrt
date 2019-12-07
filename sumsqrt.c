#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 50
#define HALFIDX (((idxbits)1<<(NMAX/2)))
#define REMAIN (((idxbits)1<<((NMAX+1)/2)))
typedef long long idxbits;
typedef long double sfloat;
sfloat sqrts[NMAX];

sfloat idxToSum(idxbits idx, int nmin, int nmax) {
    /*sum over bit in idx that's1 from nmin < x <= nmax */
    int i;
    sfloat sum = 0;
    for(i=nmin; i<nmax; i++) {
        if(idx%2)
            sum += sqrts[i];
        idx >>= 1;
    }
    return sum;
}

/*RAND_MAX = 2^31 - 1 here for gcc 9.2 on x86_64*/
#define idxRand (((idxbits)rand())*((idxbits)RAND_MAX+1) + (idxbits)rand())
void idxQSort(const sfloat *nums, idxbits *idxs, idxbits l, idxbits r) {
    /*Quick sort indexes for nums*/
    if(l < r-1) {
        idxbits i = l;
        idxbits j = r-1;
        idxbits m = l + idxRand%(r-l);
        sfloat x = nums[idxs[m]];
        m = i;
        while(m <= j) {
            idxbits tmp;
            if(nums[idxs[m]] < x){
                tmp = idxs[i];
                idxs[i++] = idxs[m];
                idxs[m++] = tmp;
            }
            else if(nums[idxs[m]] > x){
                tmp = idxs[j];
                idxs[j--] = idxs[m];
                idxs[m] = tmp;
            }
            else{
                m++;
            }
        }
        #pragma omp task
        idxQSort(nums, idxs, l, i);
        #pragma omp task
        idxQSort(nums, idxs, j+1, r);
    }
}

idxbits idxBiSect(const sfloat *nums, const idxbits *idxs, sfloat target) {
    /*return min idxs s.t. nums[idx] >= target*/
    idxbits lo=0;
    idxbits hi=HALFIDX;
    while(lo<hi){
        idxbits mid = (lo+hi)/2;
        if(nums[idxs[mid]] < target)
            lo = mid+1;
        else
            hi = mid;
    }
    if(lo < HALFIDX)
        return idxs[lo];
    return HALFIDX;
}

int main() {
    /*Have to malloc or segmentation fault*/
    sfloat *halfSum = malloc(HALFIDX*sizeof(sfloat));
    idxbits *halfIdx = malloc(HALFIDX*sizeof(idxbits));
    idxbits n;
    sfloat target;
    sfloat res;
    idxbits residx;
    for(n=0; n<NMAX; n++) {
        sqrts[n] = sqrtl(n+1);
    }
    target = idxToSum(((idxbits)1<<NMAX)-1, 0, NMAX)/2;
    res = target;
    printf("Half of sum sqrt %.16llf\n", target);

    for(n=0; n<HALFIDX; n++){
        halfIdx[n] = n;
        halfSum[n] = idxToSum(n, 0, NMAX/2);
    }
    #pragma omp parallel
    #pragma omp single
    idxQSort(halfSum, halfIdx, 0, HALFIDX);
    /* for(n=0; n<HALFIDX; n++) */
    /*     printf("%llf\n", n, halfSum[halfIdx[n]]); */

    #pragma omp parallel
    {
        idxbits n;
        sfloat res_local = res;
        idxbits residx_local;
		#pragma omp for
        for(n=0; n<REMAIN; n++){
            sfloat sum = idxToSum(n, NMAX/2, NMAX);
            idxbits l = idxBiSect(halfSum, halfIdx, target-sum);
            if(l < HALFIDX && sum + halfSum[l] - target < res_local) {
                res_local = sum + halfSum[l] - target;
                residx_local = l + (n << (NMAX/2));
            }
        }
    #pragma omp critical
        {
            if(res_local < res) {
                res = res_local;
                residx = residx_local;
            }
        }
    }

    printf("min sum of sqrt (%lld): ", residx);
    for(n=0; n<NMAX; n++){
        if(residx & ((idxbits)1<<n))
            printf("%d ", n+1);
    }
    printf(", \ndiff=%.16lle\n", res);
    free(halfSum);
    free(halfIdx);
    return 0;
}

/*
 * Output with double:
 * Half of sum sqrt 119.5179003017603776
 * min sum of sqrt (566001980646974): 2 3 4 5 6 10 11 14 16 17 18 19 20 21 22 23 24 26 28 29 32 34 35 39 40 42 50 ,
 * diff=9.9475983006414026e-14
 *
 * with long double:
 * Half of sum sqrt 119.5179003017603917
 * min sum of sqrt (566001980646974): 2 3 4 5 6 10 11 14 16 17 18 19 20 21 22 23 24 26 28 29 32 34 35 39 40 42 50 ,
 * diff=7.2719608112947753e-14
 *
 * with long double and sqrtl:
 * Half of sum sqrt 119.5179003017603922
 * min sum of sqrt (566001980614463): 1 2 3 4 5 6 9 10 11 14 17 18 19 20 21 22 23 24 26 28 29 32 34 35 39 40 42 50 ,
 * diff=7.1491423891956174e-14
 */
