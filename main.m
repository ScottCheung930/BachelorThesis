clc; 
clear; 
addpath('fig'); % 添加路径
delete(fullfile('./fig/','*.*'));
global c0 fc lambda M N Ns delta_f Ts  NR NT Fp Delta_theta figureNumber; % 定义全局变量 %CPsize
figureNumber=1;
%% Targets Information
r = [25,10, 20];  %Range (it can be a row vector)
v = [20, 40 ,30];  %Velocity
theta = [-15, -20, 15]; % Angle

%% System parameters
NR = 10; %接收天线数（与发射天线相等）
NT = 10; %发射天线数
[~, Delta_theta]=plotPattern(0,1,1);
%Delta_theta=60;
%Delta_theta为NT个振子，半波长间隔ULA天线的-10dB主瓣宽度
%Pattern为归一化阵列因子
SNR = 30; %信噪比(dB)=平均符号功率/单边白噪声功率N0
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

%% Receiver
Rx_sensing=zeros(M,Ns,NR);% Rx_sensing(m,n,:)表示第n个OFDM符号中的第m个子载波的阵列接收矢量，
RxData_sensing=zeros(M,Ns); %接收符号
P_symbol=0; %平均每符号功率

angle_dir=-60;%波束朝向角度初始为-60，每次扫描Delta_theta角度范围，直至完成[-60,60]的扫描
L=length(theta);
Ar=steeringGen(theta,NR);%steering matrix
At=steeringGen(theta,NT);
Z=zeros(1,5);
%load('./DoA/Beta.mat');
%load('./DoA/At.mat');
%load('./DoA/Ar.mat');
%% 
for kk = 1:Ndir 
   P_symbol=0;
   P_symbol_ii=0;
   [Pattern, ~]=plotPattern(angle_dir,0,0);%更新阵列因子
   %N个OFDM符号
   for ii = 1:Ns
   %M个子载波
       for jj=1:M
           %Beta=coefGen(ii,jj,r,v,theta,Pattern); % 阵列的幅度响应，多普勒，时延
           %Beta_jj_ii=reshape(Beta(jj,ii,1:L),1,L);
           %Tx=repmat(TxData(jj,ii),[NT,1]);
           %Rx=Ar*diag(Beta,1,L)*At.'*Tx; 
           Rx=Ar*diag(coefGen(ii,jj,r,v,theta,Pattern))*At.'*reshape(Tx(jj,(Ns*(kk-1)+ii),1:NT),NT,1); 
           Rx_sensing(jj,ii,1:NR)=Rx;
           P_symbol_ii=P_symbol_ii+1/NR*(Rx'*Rx);
       end
   P_symbol_ii=P_symbol_ii/M;
   P_symbol=P_symbol+P_symbol_ii;
   end
   P_symbol=P_symbol/Ns;%计算平均每符号功率
   
   %% Process a direction
    %当接收到Ns个OFDM符号时，进行一次DoA估计
       angle_start=angle_dir-0.5*Delta_theta;%搜索区间
       angle_end=angle_dir+0.5*Delta_theta;
       
       N0=P_symbol*10^(-SNR/10);%N0：单边噪声功率
       % Add Noise  (此处设信自干扰比SSIR为0)
        for iii=1:M
            for jjj=1:Ns
                for kkk=1:NR
                    Rx_sensing(iii,jjj,kkk)=Rx_sensing(iii,jjj,kkk)+sqrt(N0)...
                                  .*(randn() + 1j*randn())/sqrt(2);
                    %此处设信自干扰比SSIR为0
                end
            end
        end   
       %Rx_sensing=Rx_sensing+sqrt(N0).*(randn(M,Ns,NR)+1j*randn(M,Ns,NR))/sqrt(2);
                     
       % Received Symbol Decision
       for iii=1:Ns
         for jjj=1:M
             %[~,k]=max(conj(Rx_sensing(jjj,iii,1:NR)).*Rx_sensing(jjj,iii,1:NR));
              RxData_sensing(jjj,iii)=sum(Rx_sensing(jjj,iii,1:NR));
         end
       end
       %RxData_sensing=RxData_sensing/NR;
      
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
        %L_estimate=4;
        %% MUSIC 进行 DoA估计，得到角度谱和估计目标数
        [theta_estimate, P_music_theta, L_estimate] = MUSIC_OFDMsensing(R, Delta_theta, angle_dir); 
        [Peaks,index]=findpeaks(P_music_theta);
        if(P_music_theta(1)>0.90)
            Peaks=[Peaks,P_music_theta(1)];
            index=[index,1];
        end
        if(P_music_theta(end)>0.90)
            Peaks=[Peaks,P_music_theta(1)];
            index=[index,length(P_music_theta)];
        end
        L_estimate=min(L_estimate,length(index));
        [~,indexx]=maxk(Peaks,L_estimate);
        str=["["+num2str(angle_start)+"°,"+num2str(angle_end)+"°]"];
        fprintf(str);fprintf('\n')
        fprintf('%d target(s) detected\n', L_estimate);
        fprintf('Estimated DoA =');
        disp(theta_estimate(index(indexx)));
        
        % plot
        figure(figureNumber);
        title('MUSIC for OFDM sensing'); 
        theta_e=linspace(angle_start,angle_end);
        plot(theta_e,abs(P_music_theta)/max(abs(P_music_theta))); 
        % 在(angle_start,angle_end)之间搜索
        title(["P_M_U_S_I_C(\theta), \theta ∈ "+str])
        ylabel('Pmusic'); 
        xlabel('theta(°)'); 
        ylim([10^-3,1]);
        str=['./fig/Figure ',num2str(figureNumber),'_MUSIC.png'];
        saveas(gcf, str);
        close(gcf);
        figureNumber=figureNumber+1;
        
        %保存当前方向doa估计结果
        doa_result=vertcat(theta_estimate(index(indexx)),Peaks(indexx));
        doa_result=doa_result.';
        
        %% Periodogram 进行 Doppler-Delay估计，得到多普勒-时延谱
        %TxData_origin=TxData(:,(ii-Ns+1):ii);
        TxData_origin=TxData(:,(Ns*(kk-1)+1):(Ns*(kk-1)+ii));
        [range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData_sensing, TxData_origin);
        %[range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData, TxData);
        Kp=2^(ceil(log2(M)));% tau->range
        Mp=2^(ceil(log2(Fp*Ns)));% doppler->velocity
        velocity_1=kron(ones(Kp,1),velocity);
        range_1=kron(ones(1,Mp),range.');

        f=figure(figureNumber);
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
        
        str=['./fig/Figure ',num2str(figureNumber),'_Periodogram.png'];
        saveas(gcf, str);
        close(gcf);
        figureNumber=figureNumber+1;
        
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
        
        %保存当前方向r-v估计结果
        r_v_result=vertcat(range(maxTau(index(1:peaks_num))),velocity(maxDoppler(index(1:peaks_num))),peaks_found);
        r_v_result=r_v_result.';
        
        %保存当前结果
        result=horzcat(r_v_result, doa_result);
        if (kk==1)
            Z=result;
        else
            Z1=[Z; result];
            Z=Z1;
        end
        
        %下一个扫描方向的信道参数
        angle_dir=angle_dir+Delta_theta;
         if(angle_dir>=60)
             break;
         end

end

figure(figureNumber);
polarplot(deg2rad(Z(:,4)),Z(:,1),'*');
hold on
polarplot(deg2rad(theta),r,'o');
hold off
thetalim([-60,60]);