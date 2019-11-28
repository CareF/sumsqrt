#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 50
#define HALFIDX (((idxbits)1<<(NMAX/2)))
#define REMAIN (((idxbits)1<<((NMAX+1)/2)))
typedef long long idxbits;

double idxToSum(idxbits idx, int nmin, int nmax) {
    int i;
    double sum = 0;
    for(i=nmin; i<=nmax; i++) {
        if(idx%2)
            sum += sqrt(i);
        idx >>= 1;
    }
    return sum;
}

/*RAND_MAX = 2^31 - 1 here*/
#define idxRand (((idxbits)rand())*((idxbits)RAND_MAX+1) + (idxbits)rand())
void idxQSort(const double *nums, idxbits *idxs, idxbits l, idxbits r) {
    /*Quick sort indexes for nums*/
    if(l < r-1) {
        idxbits i = l;
        idxbits j = r-1;
        idxbits m = l + idxRand%(r-l);
        double x = nums[idxs[m]];
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

idxbits idxBiSect(const double *nums, const idxbits *idxs, double target) {
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
    double *halfSum = malloc(HALFIDX*sizeof(double));
    idxbits *halfIdx = malloc(HALFIDX*sizeof(idxbits));
    idxbits n;
    double target = idxToSum(((idxbits)1<<NMAX)-1, 1, NMAX)/2;
    double res = target;
    idxbits residx;
    printf("Half of sum sqrt %.16f\n", target);

    for(n=0; n<HALFIDX; n++){
        halfIdx[n] = n;
        halfSum[n] = idxToSum(n, 1, NMAX/2);
    }
    #pragma omp parallel
    #pragma omp single
    idxQSort(halfSum, halfIdx, 0, HALFIDX);
    /* for(n=0; n<HALFIDX; n++) */
    /*     printf("%f\n", n, halfSum[halfIdx[n]]); */

    #pragma omp parallel
    {
        double res_local = res;
        idxbits residx_local;
        for(n=0; n<REMAIN; n++){
            double sum = idxToSum(n, NMAX/2+1, NMAX);
            idxbits l = idxBiSect(halfSum, halfIdx, target-sum);
            if(l < HALFIDX && sum + halfSum[l] - target < res_local) {
                res_local = sum + halfSum[l] - target;
                residx_local = l + (n << NMAX/2);
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
    for(n=1; n<=NMAX; n++){
        if(residx & ((idxbits)1<<(n-1)))
            printf("%d ", n);
    }
    printf(", \ndiff=%.16e\n", res);
    return 0;
}

/*
 * Output:
 * Half of sum sqrt 119.5179003017603776
 * min sum of sqrt (566001980646974): 2 3 4 5 6 10 11 14 16 17 18 19 20 21 22 23 24 26 28 29 32 34 35 39 40 42 50 ,
 * diff=9.2370555648813024e-14
 */
