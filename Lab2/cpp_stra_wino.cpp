//19302010011 王海伟

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
                    //break;
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
    if(n <= 32){
        for(Addr_T i = 0; i < n; i++){
            for (Addr_T j = 0; j < n;j++){
                for(Addr_T k = 0;k < n;k++){
                    res[i * n + j] += A[i * n + k] * B[k * n + j];
                }
            }
        }
        return;
    }

    Data_T * A11 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * A12 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * A21 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * A22 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));

    Data_T * B11 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * B12 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * B21 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * B22 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));

    Data_T * C11 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * C12 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * C21 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));
    Data_T * C22 = (Data_T *) calloc((n / 2) * (n / 2), sizeof(Data_T));

    Addr_T k1 = 0;
    Addr_T k2 = 0;

    for(Addr_T j = 0; j < n/2; j++) {
        k1 = j * (n/2);
        k2 = j * n;
        for (Addr_T i = 0; i < n / 2; i+=2) {
            A11[i + k1] = A[i + k2];
            A11[i + k1 +1] = A[i + k2+1];
            B11[i + k1] = B[i + k2];
            B11[i + k1 +1] = B[i + k2 +1];
        }

        k1 = j * (n/2) - (n/2);
        for (Addr_T i = n / 2; i < n; i+=2) {
            A12[i + k1] = A[i + k2];
            B12[i + k1] = B[i + k2];
            A12[i + k1 +1] = A[i + k2 +1];
            B12[i + k1 +1] = B[i + k2 +1];
        }
    }

    for(Addr_T j = n/2; j < n; j++) {
        k1 = (j-n/2) * (n/2);
        k2 = j * n;
        for (Addr_T i = 0; i < n / 2; i+=2) {
            A21[i + k1] = A[i + k2];
            B21[i + k1] = B[i + k2];
            A21[i + k1 +1] = A[i + k2 +1];
            B21[i + k1 +1] = B[i + k2 +1];
        }
        k1 = (j-n/2) * (n/2) - (n/2);
        for (Addr_T i = n / 2; i < n; i+=2) {
            A22[i + k1] = A[i + k2];
            B22[i + k1] = B[i + k2];
            A22[i + k1+1] = A[i + k2+1];
            B22[i + k1+1] = B[i + k2+1];
        }
    }


    Data_T * Q0 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));

    Data_T * Q1 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    
    Data_T * Q2 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * T2 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * S2 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));

    Data_T * Q3 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * T3 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * S3 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));

    Data_T * Q4 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * T4 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * S4 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));

    Data_T * Q5 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * T5 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));

    Data_T * Q6 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    Data_T * S6 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
    

    for (Addr_T i = 0; i < (n/2)*(n/2); i+=2) {
        T2[i] = A21[i] + A22[i];
        S2[i] = B12[i] - B11[i];

        T4[i] = A11[i] - A21[i];
        S4[i] = B22[i] - B12[i];

        T2[i+1] = A21[i+1] + A22[i+1];
        S2[i+1] = B12[i+1] - B11[i+1];

        T4[i+1] = A11[i+1] - A21[i+1];
        S4[i+1] = B22[i+1] - B12[i+1];
    }

    for (Addr_T i = 0; i < (n/2)*(n/2); i+=2) {
        T3[i] = T2[i] - A11[i];
        S3[i] = B22[i] - S2[i];

        T3[i+1] = T2[i+1] - A11[i+1];
        S3[i+1] = B22[i+1] - S2[i+1];
    }

    for (Addr_T i = 0; i < (n/2)*(n/2); i+=2) {
        T5[i] = A12[i] - T3[i];
        S6[i] = S3[i] - B21[i];

        T5[i+1] = A12[i+1] - T3[i+1];
        S6[i+1] = S3[i+1] - B21[i+1];
    }


    stra_wino(A11, B11, n / 2, Q0);
    stra_wino(A12, B21, n / 2, Q1);
    stra_wino(T2, S2, n / 2, Q2);
    stra_wino(T3, S3, n / 2, Q3);
    stra_wino(T4, S4, n / 2, Q4);
    stra_wino(T5, B22, n / 2, Q5);
    stra_wino(A22, S6, n / 2, Q6);

    free(T2);
    free(S2);
    free(T3);
    free(S3);
    free(T4);
    free(S4);
    free(T5);
    free(S6);

    Data_T * U1 = (Data_T *) calloc((n/2) * (n/2), sizeof(Data_T));
   
    for (Addr_T i = 0; i < (n/2)*(n/2); i++) {
        U1[i] = Q0[i] + Q3[i];
    }

    for (Addr_T i = 0; i < (n/2)*(n/2); i++) {
        C11[i] = Q0[i] + Q1[i];
        C12[i] = U1[i] + Q2[i] + Q5[i];
        C21[i] = U1[i] + Q4[i] - Q6[i];
        C22[i] = U1[i] + Q4[i] + Q2[i];
    }

    for(Addr_T j = 0; j < n/2; j++) {
        for (Addr_T i = 0; i < n / 2; i++) {
            res[i + j * n] = C11[i + j * (n/2)];
        }
        for (Addr_T i = n / 2; i < n; i++) {
            res[i + j * n] = C12[i - (n/2) + j * (n/2)];
        }
    }

    for(Addr_T j = n/2; j < n; j++) {
        for (Addr_T i = 0; i < n / 2; i++) {
            res[i + j * n] = C21[i + (j-n/2) * (n/2)];
        }
        for (Addr_T i = n / 2; i < n; i++) {
            res[i + j * n] = C22[i - (n/2) + (j-n/2) * (n/2)];
        }
    }

    free(A11);
    free(A12);
    free(A21);
    free(A22);
    free(B11);
    free(B12);
    free(B21);
    free(B22);
    free(C11);
    free(C12);
    free(C21);
    free(C22);
    free(Q1);
    free(Q2);
    free(Q3);
    free(Q0);
    free(Q4);
    free(Q5);
    free(Q6);
    free(U1);


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



