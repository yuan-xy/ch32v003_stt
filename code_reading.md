> 解释代码：
  void simple_int_fft(int size)
  {
      unsigned int even, odd, span, log=0, rootindex;    // indexes
      int
  temp;
      log=0;
      int i;

      for(span=size>>1; span; span>>=1, log++)
      {

          for(odd=span;
  odd<size; odd++)         // iterate over the dual nodes
          {

              odd |= span;
  // iterate over odd blocks only
              even = odd ^ span;              // even part of the dual node pair
  //printf("even=%i,odd=%i\n",even,odd);

              temp = re[even] + re[odd];
  re[odd] = re[even] - re[odd];
              re[even] = temp;

              temp = im[even] +
  im[odd];
              im[odd] = im[even] - im[odd];
              im[even] = temp;
  rootindex = (even<<log) & (size-1); // find root of unity index
              if(rootindex)                    //
  skip rootindex[0] (has an identity)
              {
              //printf("rootindex=%i\n",rootindex);
          //
  sanity checks....
          //if (re[odd]*tr[rootindex]>(1<<30)) {printf("oops fft\n"); exit(0);}
          //if
  (re[odd]*ti[rootindex]>(1<<30)) {printf("oops fft\n"); exit(0);}
          //if (im[odd]*tr[rootindex]>(1<<30))
  {printf("oops fft\n"); exit(0);}
          //if (im[odd]*ti[rootindex]>(1<<30)) {printf("oops fft\n"); exit(0);}
  temp=re[odd]*tr[rootindex]/RES+im[odd]*ti[rootindex]/RES;
  im[odd]=im[odd]*tr[rootindex]/RES-re[odd]*ti[rootindex]/RES;
                  re[odd] = temp;
              }
  } // end of loop over n

       } // end of loop over FFT stages
       //exit(0);
  } //end of function

● 这是一个**128点基数-2整数FFT（快速傅里叶变换）**的实现。让我详细解释这个算法的每个部分：

  算法概述

  这是一个原位计算的Cooley-Tukey FFT算法，使用蝶形运算和旋转因子进行频率域变换。

  变量含义

  unsigned int even, odd, span, log=0, rootindex;  // 索引变量
  int temp;  // 临时变量

  - even/odd: 蝶形运算的偶数和奇数索引
  - span: 当前阶段的蝶形间距
  - log: 阶段计数器（对数）
  - rootindex: 旋转因子索引
  - temp: 蝶形运算的临时存储

  算法结构分析

  1. 外层循环 - FFT阶段控制

  for(span=size>>1; span; span>>=1, log++)
  - 初始: span = size/2 (对于128点FFT，span=64)
  - 每次迭代: span减半，log递增
  - 阶段数: 对于128点FFT，需要7个阶段 (2^7=128)

  2. 内层循环 - 蝶形运算

  for(odd=span; odd<size; odd++)
  处理每个阶段中的所有蝶形对。

  核心蝶形运算

  基本蝶形运算（无旋转因子）

  temp = re[even] + re[odd];       // 实部加法
  re[odd] = re[even] - re[odd];    // 实部减法
  re[even] = temp;

  temp = im[even] + im[odd];        // 虚部加法
  im[odd] = im[even] - im[odd];     // 虚部减法
  im[even] = temp;

  这相当于复数运算：
  Even' = Even + Odd
  Odd' = Even - Odd

  旋转因子乘法（带旋转因子）

  temp = re[odd]*tr[rootindex]/RES + im[odd]*ti[rootindex]/RES;
  im[odd] = im[odd]*tr[rootindex]/RES - re[odd]*ti[rootindex]/RES;
  re[odd] = temp;

  这相当于复数乘法：
  Odd' = Odd × W
  其中 W = tr[rootindex]/RES + j×ti[rootindex]/RES

  索引计算关键技巧

  1. 奇偶索引生成

  odd |= span;                    // 确保odd索引在正确的块中
  even = odd ^ span;              // 计算对应的even索引

  举例（第一阶段，span=64）：
  - odd = 64 → odd|64=64, even=64^64=0
  - odd = 65 → odd|64=65, even=65^64=1
  - 处理索引对: (0,64), (1,65), (2,66)...

  2. 旋转因子索引计算

  rootindex = (even<<log) & (size-1);

  原理: 利用位反转来找到正确的旋转因子。

  举例（第二阶段，span=32, log=1）：
  - even=0: rootindex = (0<<1) & 127 = 0 → 不使用旋转因子
  - even=2: rootindex = (2<<1) & 127 = 4 → 使用W[4]

  128点FFT的具体阶段

  ┌──────┬──────┬──────────┬────────────┬──────────────────────────┐
  │ 阶段 │ span │ 蝶形组数 │ 每组蝶形数 │ 旋转因子                 │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 1    │ 64   │ 1组      │ 64个蝶形   │ W[0]                     │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 2    │ 32   │ 2组      │ 32个蝶形   │ W[0,32]                  │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 3    │ 16   │ 4组      │ 16个蝶形   │ W[0,16,32,48]            │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 4    │ 8    │ 8组      │ 8个蝶形    │ W[0,8,16,24,32,40,48,56] │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 5    │ 4    │ 16组     │ 4个蝶形    │ W[0-60,步长4]            │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 6    │ 2    │ 32组     │ 2个蝶形    │ W[0-62,步长2]            │
  ├──────┼──────┼──────────┼────────────┼──────────────────────────┤
  │ 7    │ 1    │ 64组     │ 1个蝶形    │ W[0-63,步长1]            │
  └──────┴──────┴──────────┴────────────┴──────────────────────────┘

  整数运算优化

  分辨率处理

  /RES  // RES=8192，提供13位定点数分辨率

  优势：
  - 避免浮点运算（CH32V003无FPU）
  - 使用预计算的旋转因子表tr[]和ti[]
  - 定点数运算效率高

  算法特点

  1. 原位计算

  - 输入输出使用相同的re[]和im[]数组
  - 节省内存（重要约束：只有2KB RAM）

  2. 位反转输出

  - 输出频域数据是位反转顺序的
  - 后续处理需要重新排序（在mel_mx.h中处理）

  3. 实数输入优化

  - 输入是实信号（虚部全为0）
  - 利用实信号的对称性减少计算量

  性能考虑

  对于128点FFT：
  - 理论复杂度: O(N log N) = 128 × 7 = 896次操作
  - 实际计算: 约7个阶段，每个阶段64次蝶形运算
  - 适合实时处理: 在48MHz主频下能实现100fps处理