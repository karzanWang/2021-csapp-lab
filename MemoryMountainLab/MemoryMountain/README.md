### 环境

linux,macOS.matlab

### 说明

文件夹中有`mountain.c`和`mountain1.c`，`mountain.c`是单个矩阵自乘的版本，另外一个是两个矩阵互乘，两个矩阵互乘的MAXBYTES如果比较大的话程序会慢到无法接受，事实上对于单个矩阵自乘可以接受的大小内所生成的图和两个矩阵互乘得到的相似。

在linux环境下在代码文件夹中执行`./moutain`或者`make clean`加`make`重新编译（这里编译的是`mountain.c`，如果需要测试另一个版本可以修改`makefile`）。

然后把输出的数据放在`data.txt`中，需要复制的内容和格式在`data.txt`中已经有一份示例了。

然后在mac上终端执行`export PATH=/Applications/MATLAB_R2021a.app/bin:$PATH`(这里取决于matlab的版本)，在有data.txt和matlab代码的文件夹下执行`matlab -nodesktop -nosplash -r "display_mountain data.txt"`就可以得到图像



![mountain](/Users/karz4n/Documents/GitHub/2021-csapp-lab/MemoryMountainLab/MemoryMountain/mountain.png)