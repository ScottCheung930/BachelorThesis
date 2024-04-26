clc; 
clear; 
addpath('fig'); % 添加路径
global c0 fc lambda M N delta_f Ts CPsize NR NT Fp; % 定义全局变量

%% Transmitter 发射机
% System parameters 系统参数
c0 = 3e+8;                      % 光速(m/s)
fc = 28e+9;                     % 载波频率30GHz
lambda = c0 / fc;              % 波长(1cm)
M = 3168;                       % 子载波数
N = 112;                        % 每帧的OFDM符号数
delta_f = 120e+3;              % 子载波间隔120kHz
T = 1 / delta_f;               % OFDM符号长度 8.3us
Tcp = T / 16;                    % 保护间隔长度（数据符号长度的1/4, 2.1us）
Ts = T + Tcp;                   % 总符号长度
CPsize = M / 4;                % 循环前缀码长度 1/4
NR = 10;                        % 发射天线数量
NT = 10;                        % 接收天线数量（=发射天线数量）
Fp = 10;                        % 2D-FFT补零因子
bitsPerSymbol = 4;             % 每符号比特数
qam = 2^(bitsPerSymbol);       % QAM调制，此处使用16-QAM

% Transmit data 发送数据
data = randi([0 qam-1], M, N); % M行N列整数矩阵，使用M个子载波的N个OFDM符号
                                  % randi([0 qam-1], ___)包含从区间 [0,qam-1]   
                                  % 的均匀离散分布中得到的整数
TxData = qammod(data, qam, 'gray'); % bit->symbol，使用16QAM，Gray码
%y = qammod(data, 16, 'UnitAveragePower', true, 'PlotConstellation', true);
Tx=zeros(M,N,NT);
for ii=1:NT
    Tx(1:M,1:N,ii)=TxData;
end

%% Channel Parameters 信道参数
SNR = 30;
r = [25,10, 20];   % 距离
v = [20, 40 ,30];  % 速度
theta = [-15 -20, 15];
L=length(theta);
Ar=steeringGen(theta,NR);% 接收天线阵列的方向矩阵（steering matrix）
At=steeringGen(theta,NT);% 发射天线阵列的方向矩阵（steering matrix）
Beta=coefGen(r,v);% 多普勒与往返时延的相移 
 
%% Receiver 接收机
Rx=zeros(M,N,NR); 
% Rx(m,n,:)为第m个子载波，第n个OFDM符号被NR个接收天线接收的数据
P_symbol=0;
for ii=1:M
    P_symbol_ii=0;
    for jj=1:N
        Rx_ii_jj=Ar*diag(reshape(Beta(ii,jj,1:L),1,L))*At.'...
                  *reshape(Tx(ii,jj,1:NT),NT,1);
        Rx(ii,jj,1:NR)=Rx_ii_jj;
        P_symbol_ii=P_symbol_ii+1/NR*(Rx_ii_jj'*Rx_ii_jj);
    end
    P_symbol_ii=P_symbol_ii/N;
    P_symbol=P_symbol+P_symbol_ii;
end
% Add noise
P_symbol=P_symbol/M;     %计算平均每个天线单元的接收符号功率
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
% receive symbols
RxData=zeros(M,N);
for ii=1:M
    for jj=1:N
            RxData(ii,jj)=sum(Rx(ii,jj,1:NR));
    end
end

%% Sample Covariance Matrix 样本协方差矩阵
R=zeros(NR,NR);
for ii=1:M
    R0=zeros(NR,NR);
    for jj=1:N
       y=reshape(Rx(ii,jj,1:NR),NR,1);
       R0=R0+y*y';
    end
    R=R+R0./N;
end
R=R./M;
%%  MUSIC DoA Estimation 波达方向估计

[P_music_theta,L_estimate]=MUSIC_OFDMsensing(R);

% plot 
figure(1); 
title('MUSIC for OFDM sensing'); 
theta_e=linspace(-30,30);
plot(theta_e,abs(P_music_theta)/max(abs(P_music_theta))); 
% 在(-30，30)之间搜索，间隔20°/33
ylabel('Pmusic'); 
xlabel('theta(°)'); 
ylim([10^-3,1]); 
 
savefig('fig/figure1.fig'); % 保存图像

Peaks=imregionalmax(abs(P_music_theta));
% index=find(Peaks==1);
% PeaksValue=P_music_theta(index);
% [~,indexx]=maxk(PeaksValue,Le);
% fprintf('\nEstimated DoA =');
% disp(theta_e(index(indexx)));
% fprintf('\n')

%% Periodogram/FFT-based range-velocity estimation 周期图 距离-速度估计
[range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData, TxData);
%load('../RxData_sensing.mat');
%load('../TxData_origin.mat');
%[range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData_sensing, TxData_origin);
Kp=2^(ceil(log2(M)));% tau->range
Mp=2^(ceil(log2(Fp*N)));% doppler->velocity
velocity_1=kron(ones(Kp,1),velocity);
range_1=kron(ones(1,Mp),range.');

f=figure(2);
f.OuterPosition=[100,50,800,720];

subplot(2,2,1) %投影速度谱
mesh(range_1,velocity_1,P_TauDoppler);
xlabel('r / m');
ylabel('v / m·s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
view(90,0)
title('velocity-P');
 
subplot(2,2,2) %投影距离谱
mesh(range_1,velocity_1,P_TauDoppler);
xlabel('r / m');
ylabel('v / m·s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
view(0,0);
title('range-P');
 
subplot(2,2,3) %距离-速度
mesh(range_1,velocity_1,P_TauDoppler);
xlabel('r / m');
ylabel('v / m¡¤s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
title('range-velocity');
view(0,90);
 
subplot(2,2,4) %3D视角
mesh(range_1,velocity_1,P_TauDoppler);
xlabel('r / m');
ylabel('v / m¡¤s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
title('2D Periodogram');
 
%% Find peak
matrix=P_TauDoppler;
        matrix(matrix<0.0001*max(max(matrix)))=0; % 小于最大值0.1%的数值丢弃
        PeaksMap = imregionalmax(matrix);
        [maxTau,maxDoppler]=find(PeaksMap==1);
        peaks=zeros(1,length(maxDoppler));
        for iii=1:length(maxTau)
            peaks(iii)=matrix(maxTau(iii),maxDoppler(iii));
        end
        %先按功率升序排列
        [peaks_sort,index] = sort(peaks);
        %再降序排列
        peaks_sort=fliplr(peaks_sort);
        index=fliplr(index);
        %找到需要的若干个波峰
        peaks_num=L_estimate;%估计的目标数
        peaks_found=peaks_sort(1:peaks_num);%降序排列的估计峰值
        fprintf("Estimated range=");
        for jj=1:peaks_num
            fprintf("%f\t",range(maxTau(index(jj))) );
        end
        fprintf("\nEstimated velocity=");
        for jj=1:peaks_num
            fprintf("%f\t",velocity(maxDoppler(index(jj))) );
        end
        fprintf("\n")
        
% P_TauDoppler(P_TauDoppler<0.0001*max(max(P_TauDoppler)))=0; 
% %小于最大值0.01%的数值丢弃
% PeaksMap = imregionalmax(P_TauDoppler);
% [maxTau,maxDoppler]=find(PeaksMap==1);
% peaks=zeros(1,length(maxDoppler));
% for ii=1:length(maxTau)
%     peaks(ii)=P_TauDoppler(maxTau(ii),maxDoppler(ii));
% end
% %升序排列
% [peaks_sort,index] = sort(peaks);
% %所有波峰降序排列
% peaks_sort=fliplr(peaks_sort);
% index=fliplr(index);
% %找到需要的若干个波峰
% peaks_num=L_estimate;
% peaks_found=peaks_sort(1:peaks_num);%降序排列的估计峰值
% fprintf("Estimated range=");
% for jj=1:peaks_num
%     fprintf("%f\t",range(maxTau(index(jj))) );
% end
% fprintf("\nEstimated velocity=");
% for jj=1:peaks_num
%     fprintf("%f\t",velocity(maxDoppler(index(jj))) );
% end
% fprintf("\n")

