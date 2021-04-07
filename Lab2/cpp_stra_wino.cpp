#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

using namespace std;

typedef unsigned long long Size_T;
typedef double Data_T;
typedef long Addr_T;

const double eps = 1e-6;
const int TIMES = 3;
const double INF = 1e60;

#define min(a, b) ((a) < (b) ? (a) : (b))
#define del_arr(obj) \
    do { \
        if ((obj) != nullptr) { \
            free (obj); \
            (obj) = nullptr; \
        } \
    }while(0)

unsigned int log2_ceil(unsigned long long a){
    unsigned int nbits = __builtin_popcountl(a);
    unsigned int res = nbits > 1 ? 64 - __builtin_clzl(a) : 63 - __builtin_clzl(a);
    return res;
}




inline double clock_tdiff(struct timespec const & start, struct timespec const & end){
    return (end.tv_sec - start.tv_sec) + 1e-9 * (end.tv_nsec - start.tv_nsec);
}

void check(Size_T n, Data_T * res, Data_T * & ans){
    bool check_ok = true;
    if(ans == nullptr){
        ans = (Data_T *) calloc(n * n, sizeof(Data_T));
        for(Addr_T i = 0; i < n; i++){
            for(Addr_T j = 0;j < n; j++){
                ans[i * n + j] = res[i * n + j];
            }
        }
        printf("assign ans finished!\n");
    }else{
        for(Addr_T i = 0; i < n; i++){
            for(Addr_T j = 0;j < n; j++ ){
                if(abs(res[i * n + j] - ans[i * n + j]) > eps){
                    printf("Wrong in (%ld,%ld),res = %.6f, ans = %.6f\n",i, j, res[i * n + j], ans[i * n + j]);
                    check_ok = false;
                    break;
                }
            }
        }
        printf("check finished!\n");
    }
    return;
}

/*the matrices in mm_loop are row-majored*/
void mm_loop(Data_T * A, Data_T * B, Size_T n, Data_T * res){
    for(Addr_T i = 0; i < n; i++){
        for (Addr_T j = 0; j < n;j++){
            for(Addr_T k = 0;k < n;k++){
                res[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    return;
}

void stra_wino(Data_T * A, Data_T * B, Size_T n, Data_T * res){
    //to do...
    for(Addr_T i = 0; i < n; i++){
        for (Addr_T j = 0; j < n;j++){
            for(Addr_T k = 0;k < n;k++){
                res[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    return;
}

int main(int argc,char *argv[]){
    if(argc < 2){
        printf("usage: stra-wino n(order of matrix)");
        exit(0);
    }
    Size_T n = atoll(argv[1]);
    n = 1 << (log2_ceil(n));
    printf("round up to %lld\n",n);
    Data_T * A = (Data_T *)calloc(n*n, sizeof(Data_T));
    Data_T * B = (Data_T *)calloc(n*n, sizeof(Data_T));
    
    /*initialize A and B*/
    for(Addr_T i = 0;i < n;i++){
        timespec t;
        clock_gettime(CLOCK_MONOTONIC, &t);
        struct drand48_data rand_state;
        srand48_r(t.tv_nsec * i, &rand_state);
        for(Addr_T ii = i * n; ii < (i+1) * n; ++ii){
            drand48_r(&rand_state, &A[ii]);
        }
        for(Addr_T ii = i * n; ii < (i+1) * n;++ii){
            drand48_r(&rand_state, &B[ii]);
        }
    }

    Data_T * ans = nullptr;
    Data_T * res = (Data_T *) calloc(n * n, sizeof(Data_T));
    
    mm_loop(A, B, n, res);
    check(n,res,ans);

    double tmin;
    tmin = INF;
    timespec start, end;
    for(int times = 0; times < TIMES; times++){
        memset(res, 0, n * n * sizeof(Data_T));
        clock_gettime(CLOCK_MONOTONIC, &start);
        stra_wino(A, B, n, res);
        clock_gettime(CLOCK_MONOTONIC, &end);
        tmin = min(tmin, clock_tdiff(start, end));
    }
    check(n, res, ans);
    printf("stra_wino takes %.6f seconds\n",tmin);
 
/*
    for(int i = 0; i < n; i++){
        for(int j = 0;j < n;j++){
            printf("A[i * n + j] = %.6f, B[i * n + j] = %.6f\n",A[i * n + j], B[i * n + j]);
        }
    }

    for(int i = 0; i < n; i++){
        for(int j = 0;j < n;j++){
            printf("ans[i * n + j] = %.6f, res[i * n + j] = %.6f\n", ans[i * n + j], res[i * n + j]);
        }
    }
*/

    del_arr(A);
    del_arr(B);
    del_arr(ans);
    del_arr(res);
    return 0;
}
    


