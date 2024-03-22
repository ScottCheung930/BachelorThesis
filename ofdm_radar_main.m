clc; 
clear; 
addpath('fig'); % 添加路径

global c0 fc lambda M N delta_f Ts CPsize NR NT; % 定义全局变量

%% Transmitter
% System parameters
c0 = 3e+8;  % 光速(m/s)
fc = 30e+9;  % 载波频率，30GHz
lambda = c0 / fc; % 波长(1cm)
M = 1024;     % 子载波数（每个subframe的）
N = 112;       % subframe长度=N*符号长度
delta_f = 120e+3; % 子载波间隔，120kHz
T = 1 / delta_f; % 符号长度，8.3us
Tcp = T / 4;    % 保护间隔取符号长度的1/4, 2.1us
Ts = T + Tcp;  % 总符号时延=符号持续时间+循环前缀
CPsize = M / 4; % 循环前缀码的长度，子载波数量的1/4
NR = 10; %接收天线数（与发射天线相等）
NT = 10; %发射天线数
bitsPerSymbol = 4; % 每符号比特数
qam = 2^(bitsPerSymbol); % QAM调制，每正交方向4个幅度取值，使用16-QAM

% Transmit data
data = randi([0 qam-1], M, N); % M行N列整数矩阵，使用M个子载波的N个OFDM符号
                                                    % randi([0 qam-1], ___) 包含从区间 [0,qam-1] 的均匀离散分布中得到的整数
TxData = qammod(data, qam, 'gray'); % bit->symbol，使用16-QAM，Gray码
%y = qammod(data, 16, 'UnitAveragePower', true, 'PlotConstellation', true);
Tx=zeros(M,N,NT);
for ii=1:NT
    Tx(1:M,1:N,ii)=TxData;
end
%% Channel
SNR = 30;
r = [10,80,20,10];  % range.(it can be a row vector)
v = [0,10,-15,20];  %velocity
theta = [-pi/4, -pi/6, pi/12,pi/6];
L=length(theta);
Ar=steeringGen(theta,NR);% steering matrix of receiver
At=steeringGen(theta,NT);% steering matrix of transmitter
Beta=coefGen(r,v); % phase shift by doppler and delay 

%% Receiver
Rx=zeros(M,N,NR);
P_symbol=0;
for ii=1:M
    P_symbol_ii=0;
    for jj=1:N
        Rx_ii_jj=Ar*diag(reshape(Beta(ii,jj,1:L),1,L))*...
                      At.'*reshape(Tx(ii,jj,1:NT),NT,1);
        Rx(ii,jj,1:NR)=Rx_ii_jj;
        P_symbol_ii=P_symbol_ii+1/NR*(Rx_ii_jj'*Rx_ii_jj);
    end
    P_symbol_ii=P_symbol_ii/N;
    P_symbol=P_symbol+P_symbol_ii;
end
P_symbol=P_symbol/M;%平均每个天线单元的接收符号功率
N0=P_symbol*10^(-SNR/10);%N0：单边噪声功率
for ii=1:M
    for jj=1:N
        for kk=1:NR
            Rx(ii,jj,kk)=Rx(ii,jj,kk)+sqrt(N0)...
                               .*(randn() + 1j*randn())/sqrt(2);
            %此处设信自干扰比SSIR为0
        end
    end
end

%% Sample Covariance Matrix
R=zeros(NR,NR);
for ii=1:M
    R0=zeros(NR,NR);
    for jj=1:N
       y=reshape(Rx(ii,jj,1:NR),10,1);
       R0=R0+y*y';
    end
    R=R+R0./N;
end
R=R./M;


%% 3. Super resolution sensing method
% 3.1 MUSIC based (a time consuming but precise method)
[P_music_theta] = MUSICforOFDMsensing(R); % 调用MUSIC算法进行感知
%[P_music_range,P_music_velo] = MUSICforOFDMsensing(R); % 调用MUSIC算法进行感知
% plot the MUSIC power spectrum
figure(1); 
% subplot(3,1,1); 
% plot(linspace(0,100,length(P_music_range)),abs(P_music_range)/max(abs(P_music_range))); % 绘制距离估计的MUSIC谱
% ylabel('Pmusic'); % 设置y轴标签
% xlabel('range(m)'); % 设置x轴标签
% ylim([10^-3,1]); % 设置y轴的显示范围
% 
% subplot(3,1,2);
% plot(linspace(0,100,M),abs(P_music_velo)/max(abs(P_music_velo))); % 绘制速度估计的MUSIC谱
% ylabel('Pmusic'); % 设置y轴标签
% xlabel('velocity(m/s)'); % 设置x轴标签
% ylim([10^-3,1]); % 设置y轴的显示范围

plot(linspace(-60,60),abs(P_music_theta)/max(abs(P_music_theta))); % 绘制角度估计的MUSIC谱
ylabel('Pmusic(\theta)'); % 设置y轴标签
xlabel('\theta(°)'); % 设置x轴标签
ylim([10^-3,1]); % 设置y轴的显示范围
title('MUSIC DoA谱'); 
savefig('fig/figure3.fig'); % 保存图像到指定文件夹

% % 3.2 ESPRIT based method
% [range,velocity] = ESPRITforOFDMsensing(CIM,k); % 调用ESPRIT算法进行感知
% fprintf('The estimation result of TLS-ESPRIT is :\n'); % 打印估计结果的标题
% fprintf('Range = %f\n',range); % 打印距离估计结果
% fprintf('Velocity = %f\n',velocity); % 打印速度估计结果

