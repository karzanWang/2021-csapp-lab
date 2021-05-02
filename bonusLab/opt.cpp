#include <memory>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include <type_traits>
#include <string>
#include <iomanip>
#include <ctime>
#include<cmath>

template <int i>
double get_multi(double* a) {
    return a[0] *  get_multi<i - 1>(a + 1);
}

template <>
double get_multi<1>(double* a) {
    return a[0];
}

template<int size,int accNum,int uroll,int accCount>
class unrollTest {
private:
    double vecs[size];
    double accs[accNum];
    void initVec(){
        srand(time(0));//设置种子
        const int n = 1000;
        for (double &v : vecs) {
            v = rand() % 1000;
        }
        for (double &a : accs) {
            a = 1;
        }
    };
public:
    unrollTest(){
        initVec();
    };
    void loop(){
        double acc = 1;
        int i = 0;
        int j = 0;
        for(i = 0;i < size;i+=uroll){
            for(j = 0;j < accNum;++j){
                accs[j] = accs[j] * get_multi<accCount>(vecs + i);
            }
        }

        for(i;i < size;++i){
            accs[0] = accs[j] * vecs[i];
        }

        for(j = 0;j < accNum;++j){
            acc *= accs[j];
        }

    }
};


// ***********************************
// This is for measuring CPU clocks
#if defined(__i386__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#elif defined(__x86_64__)
static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif
// ***********************************


static void LeastSquaresFitting(int nData[],int nLen,double &a,double &b,double &r)
{
    double av_x,av_y; //声明变量
    double L_xx,L_yy,L_xy;
    double *fData = new double[nLen];
    //变量初始化
    av_x = 0; //X的平均值
    av_y = 0; //Y的平均值
    L_xx = 0; //Lxx
    L_yy = 0; //Lyy
    L_xy = 0; //Lxy
    int i = 0;
    for(i = 0; i < nLen; i++) //计算X、Y的平均值
    {
        fData[i] = log((double)nData[i]);
        av_x += i;
        av_y += fData[i];
    }
    av_x = av_x/nLen;
    av_y = av_y/nLen;

    for(i = 0; i < nLen; i+=100) //计算Lxx、Lyy和Lxy
    {
        L_xx += (i-av_x)*(i-av_x);
        L_yy += (fData[i]-av_y)*(fData[i]-av_y);
        L_xy += (i-av_x)*(fData[i]-av_y);
    }
    a = L_xy/L_xx; //斜率
    b = av_y-L_xy*av_x/L_xx; //截距
    r = double(L_xy/sqrt(L_xx*L_yy)); //相关系数r
    r *= r;
    delete fData;
}


int main() {

    int cpu_checkpoint_start, cpu_checkpoint_finish;

    int nData[4];
    int nLen = 4;
    double a,b,c;

    const int accNum = 1;
    const int uroll = 1;
    const int accCount = uroll/accNum;

    unrollTest<100, accNum, uroll, accCount> unrollTest1;
    cpu_checkpoint_start = rdtsc();
    unrollTest1.loop();
    cpu_checkpoint_finish = rdtsc();
    int avg_cpu_clocks = (cpu_checkpoint_finish - cpu_checkpoint_start);
    nData[0] = avg_cpu_clocks;

    unrollTest<200, accNum, uroll, accCount> unrollTest2;
    cpu_checkpoint_start = rdtsc();
    unrollTest2.loop();
    cpu_checkpoint_finish = rdtsc();
    avg_cpu_clocks = (cpu_checkpoint_finish - cpu_checkpoint_start);
    nData[1] = avg_cpu_clocks;

    unrollTest<300, accNum, uroll, accCount> unrollTest3;
    cpu_checkpoint_start = rdtsc();
    unrollTest3.loop();
    cpu_checkpoint_finish = rdtsc();
    avg_cpu_clocks = (cpu_checkpoint_finish - cpu_checkpoint_start);
    nData[2] = avg_cpu_clocks;

    unrollTest<400, accNum, uroll, accCount> unrollTest4;
    cpu_checkpoint_start = rdtsc();
    unrollTest4.loop();
    cpu_checkpoint_finish = rdtsc();
    avg_cpu_clocks = (cpu_checkpoint_finish - cpu_checkpoint_start);
    nData[3] = avg_cpu_clocks;


    LeastSquaresFitting(nData,nLen,a,b,c);
    std::cout<<"lan:"<<a;
    return 0;
}
