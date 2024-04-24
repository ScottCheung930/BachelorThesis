clc; 
clear; 
addpath('fig'); % 添加路径

global c0 fc lambda M N delta_f Ts CPsize NR NT Fp; % 定义全局变量

%% Transmitter
% System parameters
c0 = 3e+8;  % 光速(m/s)
fc = 28e+9;  % 载波频率，28GHz
lambda = c0 / fc; % 波长(1cm)
M = 3168;     % 子载波数
N = 112;       % 符号个数
delta_f = 120e+3; % 子载波间隔，120kHz
T = 1 / delta_f; % 符号长度，8.33us
Tcp = T / 16;    % 保护间隔取符号长度的1/16, 0.58us
Ts = T + Tcp;  % 总符号时延=符号持续时间+循环前缀
CPsize = M / 16; % 循环前缀码的长度，子载波数量的1/16
NR = 10; %接收天线数（与发射天线相等）
NT = 10; %发射天线数
Fp=10;
bitsPerSymbol = 4; % 每符号比特数
qam = 2^(bitsPerSymbol); % QAM调制，每正交方向4个幅度取值

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
SNR = 0;
r = [50,10, 40];  %[25,10, 20 range.(it can be a row vector)
v = [20, 40 ,30];  %[20, 40 ,30] velocity
theta = [-pi/12, -pi/9, pi/12];
L=length(theta);
Ar=steeringGen(theta,NR);% steering matrix of receiver
At=steeringGen(theta,NT);% steering matrix of transmitter
Beta=coefGen(r,v); % phase shift by doppler and delay 

%% Receiver
Rx=zeros(M,N,NR);% Rx(m,n,:)为第m个子载波，第n个OFDM符号被NR个接收天线接收的数据向量
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

% Add noise
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

% receive symbols
RxData=zeros(M,N);
for ii=1:M
    for jj=1:N
            RxData(ii,jj)=sum(Rx(ii,jj,1:NR));
    end
end

%% Sample Covariance Matrix
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


%% 3. Super resolution sensing method
[P_music_theta, Le] = MUSIC_OFDMsensing(R); % 调用MUSIC算法进行感知
% plot the MUSIC power spectrum
figure(1); 
theta_e=linspace(-30,30);
plot(theta_e,abs(P_music_theta)/max(abs(P_music_theta))); % 绘制角度估计的MUSIC谱
ylabel('Pmusic(\theta)'); % 设置y轴标签
xlabel('\theta(°)'); % 设置x轴标签
ylim([10^-3,1]); % 设置y轴的显示范围
title('MUSIC DoA谱'); 
savefig('fig/figure1.fig'); % 保存图像到指定文件夹

Peaks=imregionalmax(abs(P_music_theta));
index=find(Peaks==1);
PeaksValue=P_music_theta(index);
[~,indexx]=maxk(PeaksValue,Le);
fprintf('\nEstimated DoA =');
disp(theta_e(index(indexx)));
fprintf('\n')

%% 4  Periodogram/FFT-based range-doppler estimation

[range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData, TxData);

f=figure(2);
f.OuterPosition=[100,50,800,720];

subplot(2,2,1)
mesh(range,velocity,P_TauDoppler);
xlabel('r / m');
ylabel('v / m·s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
view(90,0)
title('velocity-P');

subplot(2,2,2)
mesh(range,velocity,P_TauDoppler);
xlabel('r / m');
ylabel('v / m·s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
view(0,0);
title('range-P');

subplot(2,2,3)
mesh(range,velocity,P_TauDoppler);
xlabel('r / m');
ylabel('v / m·s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
title('range-velocity');
view(0,90);

subplot(2,2,4)
mesh(range,velocity,P_TauDoppler);
xlabel('r / m');
ylabel('v / m·s^-^1');
zlabel('P');
xlim([0 80]);
ylim([0 80]);
title('2D Periodogram');
%% Find peak
matrix=P_TauDoppler;
matrix(matrix<0.0001*max(max(matrix)))=0; % 小于最大值0.1%的数值丢弃
PeaksMap = imregionalmax(matrix);
surf(double(PeaksMap));
[maxRow,maxCol]=find(PeaksMap==1);
peaks=zeros(1,length(maxCol));
for ii=1:length(maxRow)
    peaks(ii)=matrix(maxRow(ii),maxCol(ii));
end
%升序排列
[peaks_sort,I] = sort(peaks);
%所有波峰降序排列
peaks_sort=fliplr(peaks_sort);
%找到需要的若干个波峰
peaks_num=Le;%想搜索的波峰数目
maxRow1=zeros(1,peaks_num);
maxCol1=zeros(1,peaks_num);
findpeak=peaks_sort(1:peaks_num);
fprintf("Estimated [range,velocity]=\n");
for jj=1:length(findpeak)
    [maxRow1(jj),maxCol1(jj)]=find(matrix==findpeak(jj));
    fprintf("[%f,%f]\n",range(maxRow1(jj),maxCol1(jj)),velocity(maxRow1(jj),maxCol1(jj)));
end

    %savefig('fig/figure2.fig'); 