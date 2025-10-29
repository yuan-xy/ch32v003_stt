# 制作中文0-9语音codebook

## 一、清除原有样本数据

    rm ch32v003_stt/training/desktop/cb.txt
    rm ch32v003_stt/misc/cb.txt
    rm ch32v003_stt/codebook.h

## 二、重新采集中文数字音频样本
### 1. 准备工作
确保 CH32V003 开发板、麦克风模块连接正常，烧录dump_raw_audio程序（用于输出原始音频）：


    cd ch32v003_stt/training/ch32v003/dump_raw_audio
    make flash  # 需配置好WCH-Link烧录工具


程序烧录后，ch32v003会不断通过串口输出采集到的音频数据。

使用screen或picocom监听串口， 确认串口通信正常（波特率 230400）：

    sudo picocom -b 230400 -f h -p n -d 8 -s 1 /dev/ttyUSB0

### 2. 录制中文数字样本
针对汉语数字 “零、一、二、…、九”，每个数字录制几条不同样本（覆盖不同语速、音调、音量）：
每次录制流程：
串口工具显示连续数据（表示麦克风正常工作）。
对着麦克风清晰说出目标数字（如 “三”），发音后停顿 1-2 秒（避免样本重叠）。




## 三、提取中文数字特征并生成cb.txt
使用go.c程序处理原始音频，提取 MFCC 特征并标注汉语数字标签：

    cd ch32v003_stt/training/desktop/
    make
    ./go

### 标注特征样本
程序运行后会实时处理音频，当检测到语音能量超过阈值时开始记录特征，流程如下：
程序提示 “listening again....” 时，对着麦克风说一个中文数字（如 “五”）。
发音结束后，程序会检测到静音帧，此时按键盘对应数字键（如按5），将特征与标签 “5” 关联。
重复操作，为每个中文数字（0-9）标注足够样本。
标注完成后，程序会自动生成新的cb.txt，存储中文数字的特征数据。




## 四、生成汉语数字的codebook.h
使用make_codebook_h.c将cb.txt转换为嵌入式设备可直接使用的头文件：

    cd ch32v003_stt/misc
    make
    cp ../training/desktop/cb.txt .
    ./make_codebook_h

生成的codebook.h会包含中文数字 “零” 到 “九” 的特征模板，格式与原文件一致。


