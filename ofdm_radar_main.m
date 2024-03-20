clc; 
clear; 
addpath('fig'); % 添加路径

global c0 fc lambda M N delta_f Ts CPsize; % 定义全局变量

%% ISAC Transmitter
% System parameters
c0 = 3e+8;  % 光速(m/s)
fc = 30e+9;  % 载波频率，30GHz
lambda = c0 / fc; % 波长
M = 1024;     % 子载波数（每个subframe的）
N = 15;       % subframe长度=N*符号长度
delta_f = 120e+3; % 子载波间隔，120kHz
T = 1 / delta_f; % 符号长度，8.3us
Tcp = T / 4;    % 保护间隔取符号长度的1/4, 2.1us
Ts = T + Tcp;  % 总符号时延=符号持续时间+循环前缀
CPsize = M / 4; % 循环前缀码的长度，子载波数量的1/4
bitsPerSymbol = 4; % 每符号比特数
qam = 2^(bitsPerSymbol); % QAM调制，每正交方向4个幅度取值，使用16-QAM

% Transmit data
data = randi([0 qam-1], M, N); % M行N列整数矩阵，使用M个子载波的N个OFDM符号
                                                    % randi([0 qam-1], ___) 包含从区间 [0,qam-1] 的均匀离散分布中得到的整数
TxData = qammod(data, qam, 'gray'); % bit->symbol，使用16-QAM，Gray码
%y = qammod(data, 16, 'UnitAveragePower', true, 'PlotConstellation', true);
%% OFDM modulator
% 此处发送包含N个subframe的帧，即N*M长的符号序列，串并转换后通过M个平行信道发送
% M长的数据与L长的信道响应进行线性卷积，避免ISI，要加上保护间隔
% 同时由于要使用循环卷积实现，保护间隔为循环前缀
TxSignal = ifft(TxData, M); % 对每一列进行M点IFFT变换（生成时域信号）
TxSignal_cp = [TxSignal(M - CPsize + 1: M, :); TxSignal]; % 对每个OFDM符号添加循环前缀（将末尾的CPsize个采样拿到前方）
TxSignal_cp = reshape(TxSignal_cp, [], 1); % 将添加循环前缀后的信号重新整形为一维时域信号

%% Channel
% Sensing Data Generation
SNR = 30; % 设置信噪比
r = [30]; % 设置目标距离
v = [20]; % 设置目标速度
RxSignal = sensingSignalGen(TxSignal_cp, r, v, SNR); % 接收信号信号
k = length(r); % 获取目标数量

%% OFDM Radar Receivers

% 1. 2DFFT based (classical method, reference is omitted here)
Rx = RxSignal(1:size(TxSignal_cp,1),:); % 从感知信号中提取接收信号
Rx = reshape(Rx, [], N); % 将接收信号重新整形为N个OFDM符号的矩阵
Rx = Rx(CPsize + 1 : M + CPsize,:); % 去除循环前缀并重新整形为矩阵
Rx_dem = fft(Rx,M); % 对接收信号进行M点FFT
CIM_2dfft = Rx_dem .* conj(TxData); % 计算匹配滤波后的复数矩阵
RDM_2dfft = fft(ifft(CIM_2dfft,M).',10*N); % 计算距离多普勒映射

% plot the range doppler map
figure(1); % 创建一个新的图形窗口
range_2dfft = linspace(0,c0/(2*delta_f),M+1); % 计算距离向量
range_2dfft = range_2dfft(1:M); % 截取距离向量的前M个值
velocity_2dfft = linspace(0,lambda/2/Ts,10*N+1); % 计算多普勒向量
velocity_2dfft = velocity_2dfft(1:10*N); % 截取多普勒向量的前10*N个值
[X,Y] = meshgrid(range_2dfft,velocity_2dfft); % 创建网格矩阵
RDM_2dfft_norm = 10*log10(abs(RDM_2dfft)/max(abs(RDM_2dfft),[],'all')); % 归一化距离多普勒映射
mesh(X,Y,(RDM_2dfft_norm)); % 绘制三维网格图
title('2D-FFT based method'); % 设置标题
xlabel('range(m)'); % 设置x轴标签
ylabel('velocity(m/s)'); % 设置y轴标签
savefig('fig/figure1.fig'); % 保存图像到指定文件夹

% 2. CCC-based (Method proposed by Kai Wu et al.)
figure(2); % 创建第二个图形窗口
% setting parameters for CCC-based sensing method
% please refer to the paper for the meaning of these parameters
mildM = 512; % 设置参数mildM
Qbar = 64; % 设置参数Qbar
mildQ = 128; % 设置参数mildQ
% CCC
[r_cc,RDM] = cccSensing(RxSignal,TxSignal_cp,mildM,Qbar,mildQ); % 执行CCC感知方法
% plot the range doppler map
Tsa = 1/delta_f/M; % 计算时间分辨率
mildN = floor((length(TxSignal_cp)-Qbar-mildQ)/(mildM - Qbar)); % 计算参数mildN
range_ccc = linspace(0,c0/2*Tsa*mildM, mildM+1); % 计算距离向量
range_ccc = range_ccc(1:mildM); % 截取距离向量的前mildM个值
doppler_ccc = linspace(0,lambda/(mildM-Qbar)/Tsa/2,10*mildN+1); % 计算多普勒向量
range_ccc = range_ccc(1:mildM); % 截取距离向量的前mildM个值
doppler_ccc = doppler_ccc(1:10*mildN); % 截取多普勒向量的前10*mildN个值
RDM_norm = 10*log10(abs(RDM)/max(abs(RDM),[],'all')); % 归一化RDM
[X,Y] = meshgrid(range_ccc,doppler_ccc); % 创建网格矩阵
mesh(X,Y,(RDM_norm)); % 绘制三维网格图
title('CCC based method'); % 设置标题
xlabel('range(m)'); % 设置x轴标签
ylabel('velocity(m/s)'); % 设置y轴标签
savefig('fig/figure2.fig'); % 保存图像到指定文件夹

% 3. Super resolution sensing method
% 3.1 MUSIC based (a time consuming but precise method)
CIM = Rx_dem .*conj(TxData); % 计算CIM矩阵
[P_music_range,P_music_velo] = MUSICforOFDMsensing(CIM,k); % 调用MUSIC算法进行感知
% plot the MUSIC power spectrum
figure(3); % 创建第三个图形窗口
title('MUSIC for OFDM sensing'); % 设置标题
subplot(1,2,1); % 创建1行2列的子图，并定位到第一个子图
plot(linspace(0,100,length(P_music_range)),abs(P_music_range)/max(abs(P_music_range))); % 绘制距离估计的MUSIC功率谱
ylabel('Pmusic'); % 设置y轴标签
xlabel('range(m)'); % 设置x轴标签
ylim([10^-3,1]); % 设置y轴的显示范围
subplot(1,2,2); % 定位到第二个子图
plot(linspace(0,100,M),abs(P_music_velo)/max(abs(P_music_velo))); % 绘制速度估计的MUSIC功率谱
ylabel('Pmusic'); % 设置y轴标签
xlabel('velocity(m/s)'); % 设置x轴标签
ylim([10^-3,1]); % 设置y轴的显示范围
title('MUSIC for velocity estimation'); % 设置子图标题
savefig('fig/figure3.fig'); % 保存图像到指定文件夹

% 3.2 ESPRIT based method
[range,velocity] = ESPRITforOFDMsensing(CIM,k); % 调用ESPRIT算法进行感知
fprintf('The estimation result of TLS-ESPRIT is :\n'); % 打印估计结果的标题
fprintf('Range = %f\n',range); % 打印距离估计结果
fprintf('Velocity = %f\n',velocity); % 打印速度估计结果

