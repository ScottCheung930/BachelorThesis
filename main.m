clc; 
clear; 
addpath('fig'); % 添加路径

ifPlotMusicSpectrum=1;%是否画MUSIC角度谱
ifPlotPeriodogram=1;%是否画距离-速度谱

if(ifPlotPeriodogram==1)
    delete(fullfile('./fig/','*_Periodogram.png'));
end
if(ifPlotMusicSpectrum==1)
    delete(fullfile('./fig/','*_MUSIC.png'));
end
delete(fullfile('./fig/','*_Performance.png'));
global c0 fc lambda M N Ns Ndir delta_f Ts  NR NT Fp Delta_theta figureNumber; % 定义全局变量 %CPsize
figureNumber=1;
tic;

%% Targets Information
%r = [80,20, 30, 37, 25];  %Range
%v = [20, 33 ,15,10, 26];  %Velocity
%theta = [-60, -20, 15,44,-54.9]; % Angle

 r=    [ 59    61    66    45    75    49    28    35    32    53];
 v=   [ 11    31    20    15    38    39    25    16    38    10];
 theta=[    -42   -50    45    54    42   13   -34   -47    56    -2];
%  r=    [      55    65    66    66    37    69    59    72    75    43];
%   v=    [  20    37    22     5    27    15    11    31    11    10];
%   theta=    [  26   -19    59   -25   -11   -49   -15   -31   -38    19];
% r = round(55 * rand(1,10) + 25);
% v = round(40 * rand(1,10) );
% theta =round(120 * rand(1,10) - 60);
disp(r);
disp(v);
disp(theta);

%% System parameters
NR = 50; %接收天线数（与发射天线相等）
NT = 50; %发射天线数
[~, Delta_theta]=plotPattern(0,0,1);
%Delta_theta=60;
%Delta_theta为NT个振子，半波长间隔ULA天线的-10dB主瓣宽度
%Pattern为归一化阵列因子
SNR = 25; %信噪比(dB)=平均符号功率/单边白噪声功率N0
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

prepare_time=toc;
receiving_time=0;
doa_time=0;
r_v_time=0;
%% Start receiving
h = waitbar(0, 'Processing...'); % 创建进度条
flag=0;
[Pattern_comm, ~]=plotPattern(-2,0,0);
for kk = 1:Ndir 
       temp=toc;
       P_symbol=0;
       [Pattern_sense, ~]=plotPattern(angle_dir,0,0);%更新阵列因子
       Pattern=(Pattern_sense+Pattern_comm)./2;
       %N个OFDM符号
       for ii = 1:Ns
           P_symbol_ii=0;
           %M个子载波
           for jj=1:M
               Rx=Ar*diag(coefGen(ii,jj,r,v,theta,Pattern))*At.'*reshape(Tx(jj,(Ns*(kk-1)+ii),1:NT),NT,1); 
               Rx_sensing(jj,ii,1:NR)=Rx;
               P_symbol_ii=P_symbol_ii+1/NR*(Rx'*Rx);
           end
           P_symbol_ii=P_symbol_ii/M;
           P_symbol=P_symbol+P_symbol_ii;
       end
       P_symbol=P_symbol/Ns;%计算平均每符号功率
   
       %当接收到Ns个OFDM符号时，进行一次DoA估计
       angle_start=angle_dir-0.75*Delta_theta;%搜索区间
       angle_end=angle_dir+0.75*Delta_theta;
       str=["["+num2str(angle_start)+"°,"+num2str(angle_end)+"°]"];
       
       Rx_sensing=Rx_sensing/P_symbol;%归一化接收信号能量；
       N0=10^(-SNR/10);%N0：单边噪声功率
       
       % Add Noise  (此处设信自干扰比SSIR为0)
       Rx_sensing=Rx_sensing+sqrt(N0).*(randn(M,Ns,NR)+1j*randn(M,Ns,NR))/sqrt(2);
                     
       % Received Symbol Decision
       for iii=1:Ns
         for jjj=1:M
              RxData_sensing(jjj,iii)=sum(Rx_sensing(jjj,iii,1:NR));
         end
       end
       receiving_time=receiving_time+toc-temp;
%% Process a direction
       temp=toc;
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
        doa_time=doa_time+toc-temp;
    %% Parameter Estimation
        %L_estimate=4;
        
        temp=toc;
        %Periodogram 进行 Doppler-Delay估计，得到多普勒-时延谱
        TxData_origin=TxData(:,(Ns*(kk-1)+1):(Ns*(kk-1)+ii));
        [ifExist, delta_r, delta_v, range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData_sensing, TxData_origin,N0);
        r_v_time=r_v_time+toc-temp;
        
        temp=toc;
        % MUSIC 进行 DoA估计，得到角度谱和估计目标数
        if(ifExist==1)
            [theta_estimate, P_music_theta, L_estimate] = MUSIC_OFDMsensing(R, angle_start, angle_end); 
        else
            angle_dir=angle_dir+Delta_theta;
            fprintf(str);fprintf('\n')
            fprintf('0 target detected!\n\n');
            waitbar(kk/Ndir, h, sprintf('Processing... %d%%', round(kk/Ndir*100)));
            continue
        end
        
        %% Find MUSIC peaks 
%         if(L_estimate==0)
%             angle_dir=angle_dir+Delta_theta;
%             fprintf(str);fprintf('\n');
%             fprintf('L_estimate=0');
%             fprintf('\n');            
%             waitbar(kk/Ndir, h, sprintf('Processing... %d%%', round(kk/Ndir*100)));
%             continue;
%         end
        [Peaks,index]=findpeaks(P_music_theta,'MinPeakProminence',5);
        if(P_music_theta(1)>-0.1 && P_music_theta(1)-min(P_music_theta)>10)
            Peaks=[Peaks,P_music_theta(1)];
            index=[index,1];
        end
        if(P_music_theta(end)>-0.1&& P_music_theta(end)-min(P_music_theta)>10)
            Peaks=[Peaks,P_music_theta(1)];
            index=[index,length(P_music_theta)];
        end
        L_estimate=min(L_estimate,length(index));
        if(L_estimate==0)
            angle_dir=angle_dir+Delta_theta;
            fprintf(str);fprintf('\n');
            fprintf('can not find siginficant peaks!\n\n');
            waitbar(kk/Ndir, h, sprintf('Processing... %d%%', round(kk/Ndir*100)));
            continue;
        end
        [~,indexx]=maxk(Peaks,L_estimate);
        str=["["+num2str(angle_start)+"°,"+num2str(angle_end)+"°]"];
        fprintf(str);fprintf('\n')
        fprintf('%d target(s) detected\n', L_estimate);
        fprintf('Estimated DoA =');
        disp(theta_estimate(index(indexx)));
        
        if (ifPlotMusicSpectrum==1)
            plotMusicSpectrum(theta_estimate,P_music_theta,str);
        end
        
        %保存当前方向doa估计结果
        doa_result=vertcat(theta_estimate(index(indexx)),Peaks(indexx));
        doa_result=doa_result.';
        doa_time=doa_time+toc-temp;
        %%  Find Tau-Doppler peaks
        temp=toc;
        matrix=P_TauDoppler;
        [max_P,V_Ind]=max(matrix.');
        
        %max_P=conv(max_P,[0.33,0.33,0.33]);
        [pks,maxTau]=findpeaks(max_P,'MinPeakProminence',7,'MinPeakDistance',2,'MinPeakHeight',max(max_P)-100);
        meanHeight=mean(pks);
        [pks,index]=findpeaks(pks,'MinPeakProminence',7,'MinPeakHeight',75);
        maxTau=maxTau(index);
        %[pks,maxTau]=findpeaks(max_P,'MinPeakDistance',3,'MinPeakHeight',max(max_P)-100);
        %[pks,maxTau]=findpeaks(max_P,'MinPeakDistance',3,'MinPeakHeight',max(max_P)-50);
        %directionality=length(pks);
        %index=find(pks>max(pks)-22);
        %pks= pks(index);
        %maxTau=maxTau(index);
        %[pks,maxTau]=findpeaks(max_P,'MinPeakDistance',1,'MinPeakHeight',max(max_P)-22);
       
        maxDoppler=V_Ind(maxTau);
        peaks_num=length(pks);
        
        fprintf("Estimated range=");
        for jj=1:peaks_num
            fprintf("%f\t",range(maxTau(jj) ));
        end
        fprintf("\nEstimated velocity=");
        for jj=1:peaks_num
            fprintf("%f\t",velocity(maxDoppler(jj)) );
        end
        fprintf("\n\n");
        
        if (ifPlotPeriodogram==1)
            plotPeriodogram(range,velocity,P_TauDoppler,str);
        end
        
        %保存当前方向r-v估计结果
        %r_v_result=vertcat(range(maxTau(1:peaks_num)),velocity(maxDoppler(1:peaks_num)),pks,repmat(directionality,1,peaks_num));
        r_v_result=vertcat(range(maxTau(1:peaks_num)),velocity(maxDoppler(1:peaks_num)),pks,round(10*pks./meanHeight));
        r_v_result=r_v_result.';
        
        r_v_time=r_v_time+toc-temp;
        
        %% 保存当前结果
        if(peaks_num>L_estimate)
            doa_result=repmat(doa_result,ceil(peaks_num/L_estimate),1);
            doa_result=doa_result(1:peaks_num,:);
        end
        result=horzcat(r_v_result, doa_result);
        if (flag==0)
            flag=1;
            Z=result;
        else
            Z1=[Z; result];
            Z=Z1;
        end
        
        %下一个扫描方向的信道参数
        angle_dir=angle_dir+Delta_theta;
        
        waitbar(kk/Ndir, h, sprintf('Processing... %d%%', round(kk/Ndir*100)));
end

%% 去除冗余目标点
Zsort=sortrows(Z,[1,-3]);
%Zsort(:,4)=1./Zsort(:,4);
%Zsort(:,4)=floor(abs(Zsort(:,4)-mean(Zsort(:,4)))/(2*std(Zsort(:,4))));
%Zsort(:,4)=round(10*Zsort(:,4));
%Zsort=sortrows(Zsort,[1,-4,-3]);%先按距离r升序，然后按相对峰值功率和绝对峰值降序排列
Zsort=sortrows(Zsort,[1,-4,-3]);
Zprun=Zsort(1,:);
Target_Num=1;
[Num0,~]=size(Zsort);
delta_v=3*delta_v;
for i=2:Num0
    count=0;
    for j=1:Target_Num
        if ((Zprun(j,1)-delta_r<Zsort(i,1))&&(Zprun(j,1)+delta_r>Zsort(i,1))&&...
            (Zprun(j,2)-delta_v<Zsort(i,2))&&(Zprun(j,2)+delta_v>Zsort(i,2)))
            count=count+1;
            break;
        end
    end
    if(count==0)
        Target_Num=Target_Num+1;
        Zprun=[Zprun;Zsort(i,:)];
    end
end
fprintf("after pruning, TargetNum=%d\n",Target_Num);

%% 
figure(figureNumber);
subplot(1,2,1)
polarplot(deg2rad(Z(:,5)),Z(:,1),'*');
hold on
polarplot(deg2rad(theta),r,'o');
hold off
thetalim([-60,60]);
title("Before Pruning");
subplot(1,2,2)
polarplot(deg2rad(Zprun(:,5)),Zprun(:,1),'*');
hold on
polarplot(deg2rad(theta),r,'o');
hold off
thetalim([-60,60]);
title("After Pruning");
str=['./fig/Figure ',num2str(figureNumber),'_Performance.png'];
saveas(gcf, str);
%close(h);
total_time=toc;
fprintf("总耗时：%.2f秒\n",total_time);

function plotPeriodogram(range,velocity,P_TauDoppler,str)
    global figureNumber;
    Kp=length(range);
    Mp=length(velocity);
    velocity_1=kron(ones(Kp,1),velocity);
    range_1=kron(ones(1,Mp),range.');

    f=figure(figureNumber);
    f.OuterPosition=[100,50,800,720];

    subplot(2,2,1) %投影速度谱
    mesh(range_1,velocity_1,P_TauDoppler);
    xlabel('r / m');
    ylabel('v / m·s^-^1');
    zlabel('P');
    ylim([0 60]);
    xlim([10 90]);
    view(90,0)
    title(["P(v), \theta ∈ "+str]);

    subplot(2,2,2) %投影距离谱
    mesh(range_1,velocity_1,P_TauDoppler);
    xlabel('r / m');
    ylabel('v / m·s^-^1');
    zlabel('P');
    ylim([0 60]);
    xlim([10 90]);
    view(0,0);
    title(["P(r), \theta ∈ "+str]);

    subplot(2,2,3) %距离-速度
    mesh(range_1,velocity_1,P_TauDoppler);
    xlabel('r / m');
    ylabel('v / m·s^-^1');
    zlabel('P');
    ylim([0 60]);
    xlim([10 90]);
    title(["P(r,v) top view, \theta ∈ "+str]);
    view(0,90);

    subplot(2,2,4) %3D视角
    mesh(range_1,velocity_1,P_TauDoppler);
    xlabel('r / m');
    ylabel('v / m·s^-^1');
    zlabel('P');
    ylim([0 60]);
    xlim([10 90]);
    title(["P(r,v), \theta ∈ "+str]);

    str=['./fig/Figure ',num2str(figureNumber),'_Periodogram.png'];
    saveas(gcf, str);
    close(gcf);
    figureNumber=figureNumber+1;
end

function plotMusicSpectrum(theta_estimate,P_music_theta,str)
    global figureNumber;
    figure(figureNumber);
    title('MUSIC for OFDM sensing'); 
    plot(theta_estimate,P_music_theta); 
    % 在(angle_start,angle_end)之间搜索
    title(["P_M_U_S_I_C(\theta), \theta ∈ "+str])
    ylabel('P_M_U_S_I_C (dB)'); 
    xlabel('\theta (°)'); 
    ylim([-70, 0]);
    str=['./fig/Figure ',num2str(figureNumber),'_MUSIC.png'];
    saveas(gcf, str);
    close(gcf);
    figureNumber=figureNumber+1;
end