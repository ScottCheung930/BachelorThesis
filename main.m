clc; 
clear; 
addpath('fig'); % 添加路径

global c0 fc lambda M N Ns delta_f Ts  NR NT Fp Delta_theta; % 定义全局变量 %CPsize
%% Targets Information
r = [25,10, 20];  %Range (it can be a row vector)
v = [20, 40 ,30];  %Velocity
theta = [-pi/12, -pi/9, pi/12]; % Angle

%% System parameters
NR = 10; %接收天线数（与发射天线相等）
NT = 10; %发射天线数
[Pattern, Delta_theta]=plotPattern();
%Delta_theta为NT个振子，半波长间隔ULA天线的-10dB主瓣宽度
%Pattern为归一化阵列因子
SNR = 10; %信噪比(dB)=平均符号功率/单边白噪声功率N0
Ndir = ceil(120/Delta_theta); %为了扫描[-60,60]的范围需要扫描Ndir个方向
K = 1120;      % 每帧的OFDM符号数
Ns = 112;     % 在信号波束内进行一次DoA估计所需的OFDM符号数    
Nf = ceil(Ns*Ndir/K); % 总帧数个数
N=Nf*K; %总OFDM符号数

c0 = 3e+8;  % 光速(m/s)
fc = 28e+9;  % 载波频率，28GHz
lambda = c0 / fc; % 波长(1cm)
delta_f = 120e+3; % 子载波间隔，120kHz
M = 3168;     % 子载波数
T = 1 / delta_f; % 符号长度，8.33us
Tcp = T / 16;    % 保护间隔取符号长度的1/16, 0.58us
Ts = T + Tcp;  % 总符号时延=符号持续时间+循环前缀
%CPsize = M / 16; % 循环前缀码的长度，子载波数量的1/16

Fp=10;
bitsPerSymbol = 4; % 每符号比特数
qam = 2^(bitsPerSymbol); % QAM调制，每正交方向bitsPerSymbo个幅度取值

%% Generate Transmission data
data = randi([0 qam-1], M, N); % M行N列整数矩阵，使用M个子载波的N个OFDM符号
                                                    % randi([0 qam-1], ___) 包含从区间 [0,qam-1] 的均匀离散分布中得到的整数
TxData = qammod(data, qam, 'gray'); % bit->symbol，使用16-QAM，Gray码
%y = qammod(data, 16, 'UnitAveragePower', true, 'PlotConstellation', true);

% Generate 
Tx=zeros(M,N,NT);
for ii=1:NT
    Tx(1:M,1:N,ii)=TxData;
end
%% Set Channel Parameters
L=length(theta);


%% Receiver
Tx=zeros(NT,1);
Rx=zeros(NR,1);
Rx_sensing=zeros(M,Ns,NR);% Rx_sensing(m,n,:)表示第n个OFDM符号中的第m个子载波的阵列接收矢量，
RxData_sensing=zeros(M,Ns); %接收符号，本仿真不考虑CSI，因而对于接收分集，选择功率最大的信道
P_symbol=0; %平均每符号功率

angle_dir=-60;%波束朝向角度初始为-60，每次扫描Delta_theta角度范围，直至完成[-60,60]的扫描
Ar=steeringGen(theta-deg2rad(angle_dir),NR);%初始steering matrix
At=steeringGen(theta-deg2rad(angle_dir),NT);
Z=zeros(1,5);

%依次接收N个OFDM符号
for ii = 1:N 
   P_symbol_ii=0;
   
   %依次计算M个子载波的接收信号
   for jj=1:M
       Beta=coefGen(ii,jj,r,v,theta-deg2rad(angle_dir)); % 阵列的幅度响应，多普勒，时延
       
       Tx=repmat(TxData(jj,ii),[NT,1]);
       Rx=Ar*diag(reshape(Beta,1,L))*...
                      At.'*Tx; 
       Rx_sensing(jj,mod(ii,Ns)+1,1:NR)=Rx;
       P_symbol_ii=P_symbol_ii+1/NR*norm(Rx);
   end
   P_symbol_ii=P_symbol_ii/M;
   P_symbol=P_symbol+P_symbol_ii;
   

   %% Process a direction
   if(mod(ii,Ns)==0)%当接收到Ns个OFDM符号时，进行一次DoA估计
       angle_start=angle_dir-0.5*Delta_theta;%搜索区间
       angle_end=angle_dir+0.5*Delta_theta;
       
       P_symbol=P_symbol/Ns;%计算平均每符号功率
       N0=P_symbol*10^(-SNR/10);%N0：单边噪声功率
       % Add Noise  (此处设信自干扰比SSIR为0)
       Rx_sensing=Rx_sensing+sqrt(N0)...
                         .*(randn(M,Ns,NR)+1j*randn(M,Ns,NR))/sqrt(2);
                     
       % Received Symbol Decision
       for iii=1:Ns
         for jjj=1:M
             %[~,k]=max(conj(Rx_sensing(jjj,iii,1:NR)).*Rx_sensing(jjj,iii,1:NR));
              RxData_sensing(jjj,iii)=sum(Rx_sensing(jjj,iii,1:NR));
         end
       end
      
       % Covariance Matrix
        R=zeros(NR,NR);
        for jjj=1:M
            R0=zeros(NR,NR);
            for iii=1:Ns
               y=reshape(Rx_sensing(jjj,iii,1:NR),NR,1);
               R0=R0+y*y';
            end
            R=R+R0./Ns;
        end
        R=R./M;
        
        %% Parameter Estimation
        % MUSIC 进行 DoA估计，得到角度谱和估计目标数
        [theta_estimate, P_music_theta, L_estimate] = MUSIC_OFDMsensing(R, Delta_theta, angle_dir); 
        Peaks=imregionalmax(abs(P_music_theta));
        index=find(Peaks==1);
        PeaksValue=P_music_theta(index);
        [~,indexx]=maxk(PeaksValue,L_estimate);
        str=["in["+num2str(angle_start)+"°,"+num2str(angle_end)+"°], \n"];
        fprintf(str);
        fprintf('%d target(s) detected\n', L_estimate);
        fprintf('Estimated DoA =');
        disp(theta_estimate(index(indexx)));
        
        %保存当前方向doa估计结果
        doa_result=vertcat(theta_estimate(index(indexx)),PeaksValue(indexx));
        doa_result=doa_result.';
        
        % Periodogram 进行 Doppler-Delay估计，得到多普勒-时延谱
        [range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData_sensing, TxData(:,ii-Ns+1:ii));
        matrix=P_TauDoppler;
        matrix(matrix<0.0001*max(max(matrix)))=0; % 小于最大值0.1%的数值丢弃
        PeaksMap = imregionalmax(matrix);
        %surf(double(PeaksMap));
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
        fprintf("\n");
        
        %保存当前方向r-v估计结果
        r_v_result=vertcat(range(maxTau(index(1:peaks_num))),velocity(maxDoppler(index(1:peaks_num))),peaks_found);
        r_v_result=r_v_result.';
        
        %保存当前结果
        result=horzcat(r_v_result, doa_result);
        if (ii==Ns)
            Z=result;
        else
            Z1=[Z; result];
            Z=Z1;
        end
        
%         fprintf("Estimated [range,velocity]=\n");
%         for jj=1:peaks_num
%             range(maxTau(index(jj)));
%             [maxRow1(jj),maxCol1(jj)]=find(matrix==peaks_found(jj));
%             fprintf("[%f,%f]  ",range(maxRow1(jj),maxCol1(jj)),velocity(maxRow1(jj),maxCol1(jj)));
%         end
        
        
        %下一个扫描方向的信道参数
        angle_dir=angle_dir+Delta_theta;
        if(angle_dir>60)
            break;
        end
        Ar=steeringGen(theta-deg2rad(angle_dir),NR);%更新steering matrix
        At=steeringGen(theta-deg2rad(angle_dir),NT);
   end
end

figure(2);
polarplot(deg2rad(Z(:,4)),Z(:,1),'*');
thetalim([-60,60]);

% Rx=zeros(M,N,NR);% Rx(m,n,:)为第m个子载波，第n个OFDM符号被NR个接收天线接收的数据向量
% P_symbol=0;
% for ii=1:M
%     P_symbol_ii=0;
%     for jj=1:N
%         Rx_ii_jj=Ar*diag(reshape(Beta(ii,jj,1:L),1,L))*...
%                       At.'*reshape(Tx(ii,jj,1:NT),NT,1);
%         Rx(ii,jj,1:NR)=Rx_ii_jj;
%         P_symbol_ii=P_symbol_ii+1/NR*(Rx_ii_jj'*Rx_ii_jj);
%     end
%     P_symbol_ii=P_symbol_ii/N;
%     P_symbol=P_symbol+P_symbol_ii;
% end
% 
% % Add noise
% P_symbol=P_symbol/M;%平均每个天线单元的接收符号功率
% N0=P_symbol*10^(-SNR/10);%N0：单边噪声功率
% for ii=1:M
%     for jj=1:N
%         for kk=1:NR
%             Rx(ii,jj,kk)=Rx(ii,jj,kk)+sqrt(N0)...
%                                .*(randn() + 1j*randn())/sqrt(2);
%             %此处设信自干扰比SSIR为0
%         end
%     end
% end

% % receive symbols
% RxData=zeros(M,N);
% for ii=1:M
%     for jj=1:N
%             RxData(ii,jj)=sum(Rx(ii,jj,1:NR));
%     end
% end



% %% 3. Super resolution sensing method
% [P_music_theta, Le] = MUSIC_OFDMsensing(R); % 调用MUSIC算法进行感知
% % plot the MUSIC power spectrum
% figure(1); 
% theta_e=linspace(-30,30);
% plot(theta_e,abs(P_music_theta)/max(abs(P_music_theta))); % 绘制角度估计的MUSIC谱
% ylabel('Pmusic(\theta)'); % 设置y轴标签
% xlabel('\theta(°)'); % 设置x轴标签
% ylim([10^-3,1]); % 设置y轴的显示范围
% title('MUSIC DoA谱'); 
% savefig('fig/figure1.fig'); % 保存图像到指定文件夹
% 
% Peaks=imregionalmax(abs(P_music_theta));
% index=find(Peaks==1);
% PeaksValue=P_music_theta(index);
% [~,indexx]=maxk(PeaksValue,Le);
% fprintf('\nEstimated DoA =');
% disp(theta_e(index(indexx)));
% fprintf('\n')
% 
% %% 4  Periodogram/FFT-based range-doppler estimation
% 
% [range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData, TxData);
% 
% f=figure(2);
% f.OuterPosition=[100,50,800,720];
% 
% subplot(2,2,1)
% mesh(range,velocity,P_TauDoppler);
% xlabel('r / m');
% ylabel('v / m·s^-^1');
% zlabel('P');
% xlim([0 80]);
% ylim([0 80]);
% view(90,0)
% title('velocity-P');

% subplot(2,2,2)
% mesh(range,velocity,P_TauDoppler);
% xlabel('r / m');
% ylabel('v / m·s^-^1');
% zlabel('P');
% xlim([0 80]);
% ylim([0 80]);
% view(0,0);
% title('range-P');
% 
% subplot(2,2,3)
% mesh(range,velocity,P_TauDoppler);
% xlabel('r / m');
% ylabel('v / m·s^-^1');
% zlabel('P');
% xlim([0 80]);
% ylim([0 80]);
% title('range-velocity');
% view(0,90);
% 
% subplot(2,2,4)
% mesh(range,velocity,P_TauDoppler);
% xlabel('r / m');
% ylabel('v / m·s^-^1');
% zlabel('P');
% xlim([0 80]);
% ylim([0 80]);
% title('2D Periodogram');


    %savefig('fig/figure2.fig'); 