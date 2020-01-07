#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NMAX 50
#define HALFIDX ((idxbits)1<<(NMAX/2))
#define REMAIN ((idxbits)1<<((NMAX+1)/2))
typedef long long idxbits;
typedef long double sfloat;
typedef struct {
    idxbits idx;
    sfloat val;
} sqsum;

/*RAND_MAX = 2^31 - 1 here for gcc 9.2 on x86_64*/
#define idxRand (((idxbits)rand())*((idxbits)RAND_MAX+1) + (idxbits)rand())
void QSort(sqsum *nums, idxbits l, idxbits r) {
    /*Quick sort indexes for nums*/
    if(l < r-1) {
        idxbits i = l;
        idxbits j = r-1;
        idxbits m = l + idxRand%(r-l);
        sfloat x = nums[m].val;
        m = i;
        while(m <= j) {
            sqsum tmp;
            if(nums[m].val < x){
                tmp = nums[i];
                nums[i++] = nums[m];
                nums[m++] = tmp;
            }
            else if(nums[m].val > x){
                tmp = nums[j];
                nums[j--] = nums[m];
                nums[m] = tmp;
            }
            else{
                m++;
            }
        }
        #pragma omp task
        QSort(nums, l, i);
        #pragma omp task
        QSort(nums, j+1, r);
    }
}

void merge(sqsum *nums, idxbits l, idxbits m, idxbits r) {
    sqsum *ary = malloc((r-l)*sizeof(sqsum));
    idxbits i = l, j = m, t = 0;
    while(i < m && j < r) {
        if(nums[i].val <= nums[j].val)
            ary[t++] = nums[i++];
        else
            ary[t++] = nums[j++];
    }
    while(i < m)
        ary[t++] = nums[i++];
    while(j < r) 
        ary[t++] = nums[j++];

    t = 0;
    while(l < r)
        nums[l++] = ary[t++];
    free(ary);
}

void MergeSort_s(sqsum *nums, idxbits l, idxbits r) {
    if(l < r-1) {
        idxbits m = (l + r)/2;
        MergeSort_s(nums, l, m);
        MergeSort_s(nums, m, r);
        merge(nums, l, m, r);
    }
    return;
}

void MergeSort(sqsum *nums, idxbits l, idxbits r) {
    if(l < r-32) {
        idxbits m = (l + r)/2;
        #pragma omp task
        MergeSort(nums, l, m);
        #pragma omp task
        MergeSort(nums, m, r);
        #pragma omp taskwait
        merge(nums, l, m, r);
    }
    else 
        MergeSort_s(nums, l, r);
    return;
}

int main() {
    sqsum *half = malloc(HALFIDX*sizeof(sqsum));
    sqsum *remain = malloc(REMAIN*sizeof(sqsum));
    int n;
    sfloat target;
    sqsum res;

    half[0].val = 0.0;
    half[0].idx = 0;
    for(n = 1; n <= NMAX/2; n++){
        /* halfSum[1<<(n-1):1<<n] = sqrt(n) + halfSum[:1<<(n-1)] */
        idxbits t;
        const sfloat sqrtn = sqrtl(n);
        const idxbits offcut = (idxbits)1<<(n-1);
        sqsum * const newarray = half + offcut;
        /* #pragma omp parallel for */
        for(t = 0; t < offcut; t++) {
            newarray[t].idx = half[t].idx + offcut;
            newarray[t].val = half[t].val + sqrtn;
        }
    }
    #pragma omp parallel
    #pragma omp single
    /* MergeSort(half, 0, HALFIDX); */
    QSort(half, 0, HALFIDX);

    remain[0].val = 0.0;
    remain[0].idx = 0;
    for(n = NMAX/2+1; n <= NMAX; n++){
        idxbits t;
        const sfloat sqrtn = sqrtl(n);
        const idxbits offcut = (idxbits)1<<(n-NMAX/2-1);
        sqsum * const newarray = remain + offcut;
        /* #pragma omp parallel for */
        for(t = 0; t < offcut; t++) {
            newarray[t].idx = remain[t].idx + ((idxbits)1<<(n-1));
            newarray[t].val = remain[t].val + sqrtn;
        }
    }
    target = (half[HALFIDX-1].val + remain[REMAIN-1].val)/2.0;
    printf("Looking for sum close to %.16Lf\n", target);

    res.val = target;
    #pragma omp parallel
    {
        idxbits t;
        sqsum res_local = res;
        #pragma omp for
        for(t = 0; t < REMAIN; t++){
            sfloat m = target - remain[t].val;
            idxbits lo = 0;
            idxbits hi = HALFIDX;
            while(lo < hi){
                idxbits mid = (lo + hi)/2;
                if(half[mid].val < m)
                    lo = mid + 1;
                else 
                    hi = mid;
            }
            if(lo < HALFIDX && half[lo].val - m < res_local.val) {
                res_local.val = half[lo].val - m;
                res_local.idx = half[lo].idx + remain[t].idx;
            }
        }
        #pragma omp critical
        {
            if(res_local.val < res.val) {
                res = res_local;
            }
        }
    }

    printf("min sum of sqrt (%lld), delta %.16Le\n", res.idx, res.val);

    free(half);
    free(remain);
    return 0;
}
