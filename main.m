clc; 
clear; 
addpath('fig'); % 添加路径

ifPlotMusicSpectrum=0;%是否画MUSIC角度谱
ifPlotPeriodogram=0;%是否画距离-速度谱

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
%% System parameters
    NR = 50; %接收天线数（与发射天线相等）
    NT = 50; %发射天线数
    [~, Delta_theta]=plotAF(0,0,1);
    %Delta_theta=60;
    %Delta_theta为NT个振子，半波长间隔ULA天线的-10dB主瓣宽度
    %Pattern为归一化阵列因子

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

    Fp=5;
    bitsPerSymbol = 4; % 每符号比特数
    qam = 2^(bitsPerSymbol); % QAM调制，每正交方向bitsPerSymbo个幅度取值


%% 蒙特卡洛实验
%SNR = 0; %信噪比(dB)=平均符号功率/单边白噪声功率N0SNR = 0; %信噪比(dB)=平均符号功率/单边白噪声功率N0
Num_SNR=6;
Index_SNR=1;
Total_MC_time=10;
Total_Ndir=Num_SNR*Ndir*Total_MC_time;
RMSE_DoA=zeros(1,Num_SNR);
RMSE_range=zeros(1,Num_SNR);
RMSE_velocity=zeros(1,Num_SNR);
RMSE_pos=zeros(1,Num_SNR);
Pr_detect=zeros(1,Num_SNR);

%开启并行池
if isempty(gcp('nocreate'))%如果并行未开启
    p = parcluster('local');  % 获取当前的本地集群配置
    p.NumWorkers = 12; 
    parpool(4);       %开启6个并行工作池
end

 h = waitbar(0, 'Processing...'); % 创建进度条
prepare_time=toc;
receiving_time=0;
doa_time=0;
r_v_time=0;
 %开始循环
 for MC_time=0:(Total_MC_time-1)
    Index_SNR=1;
    for SNR=-5:5:20
        temp=toc;
        %% Targets Information
        % r = [80,20, 30, 37, 25];  %Range
        % v = [20, 23 ,27,32, 38];  %Velocity
        % theta = [-60, -20, 15,44,-54.9]; % Angle
    %     r=    [ 59    61    66    45    75    49    28    35    32    53];
    %     v=   [ 11    31    20    15    38    39    25    16    38    10];
    %     theta=[ -42   -50    45    54    42   13   -34   -47    56    -2];
        L=10;
        r_max = 80; r_min=25;% 区间界
        k = 2 / (r_max^2-r_min^2); % 距离概率密度的常数
        % 生成[25,80]的随机距离
        r = zeros(1, L);
        i=1;
        while i <= L
            r_temp = rand() * (r_max-r_min)+r_min;
            f_r = k * r_temp;
            if rand <= f_r
                r(i) = r_temp;
                i=i+1;
            else
                continue;% 如不符合概率密度函数要求，则重新生成
            end
        end
        v = 40 * rand(1,L);
        theta =round(120 * rand(1,L) - 60);
        disp(r);
        disp(v);
        disp(theta);
        h_pathloss = c0./((fc).*(r).^2); %路径损耗

        delay = 2 .* r ./ c0; %各目标的往返时延。
        doppler = 2*v/lambda;    %各目标的多普勒频移，这里的速度是在LoS方向的分量
        h_doppler=exp(1j*2*pi*Ts.*kron((0:(Ns-1)).',doppler));
        h_delay=exp(-1j*2*pi*delta_f.*kron((0:(M-1)).',delay));
        L=length(theta);
        Ar=steeringGen(theta,NR);%steering matrix
        At=steeringGen(theta,NT);
        P_symbol=0; %平均每符号功率
        N0=10^(-SNR/10);%N0：单边噪声功率

        angle_sense=-60;%波束朝向角度初始为-60，每次扫描Delta_theta角度范围，直至完成[-60,60]的扫描
        angle_comm=theta(end);%令最后一个位置为通信用户
        [Pattern_comm, ~]=plotAF(angle_comm,0,0);
        pos=0:(NT-1);
        beamforming_comm=exp(-1j*pi*sind(angle_comm).*(pos.'));
        %% Generate Transmission data
        data = randi([0 qam-1], M, N); % M行N列整数矩阵，使用M个子载波的N个OFDM符号
                                                            % randi([0 qam-1], ___) 包含从区间 [0,qam-1] 的均匀离散分布中得到的整数
        TxData = qammod(data, qam, 'gray'); % bit->symbol，使用16-QAM，Gray码
        %y = qammod(data, 16, 'UnitAveragePower', true, 'PlotConstellation', true);

        Tx=zeros(M,N,NT);
        for ii=1:NT
            Tx(:,:,ii)=TxData;
        end
        %% Receiver
        Rx_sensing=zeros(M,Ns,NR);% Rx_sensing(m,n,:)表示第n个OFDM符号中的第m个子载波的阵列接收矢量，
        RxData_sensing=zeros(M,Ns); %接收符号

        Z=zeros(1,5);
       prepare_time=prepare_time+toc-temp;
        %% Start receiving
        flag=0;
        %扫描Ndir个方向
        for kk = 1:Ndir
                temp=toc;
               %下一个扫描方向
                if(kk~=1)
                    angle_sense=angle_sense+Delta_theta;
                    waitbar((kk+MC_time*Num_SNR*Ndir)/Total_Ndir, h,...
                        sprintf('Processing... %d%%，remain：%.2f min', ...
                        round(100*(kk+MC_time*Num_SNR*Ndir)/Total_Ndir),...
                        (Total_Ndir-(Ndir*Num_SNR*MC_time+kk))/(kk+MC_time*Num_SNR*Ndir)*(temp/60)));
                end

               %当接收到Ns个OFDM符号时，进行一次DoA估计
               angle_start = angle_sense-0.65*Delta_theta;%搜索区间
               angle_end  = angle_sense+0.65*Delta_theta;
               str=["["+num2str(angle_start)+"°,"+num2str(angle_end)+"°]"];
               fprintf("Within "); fprintf(str); fprintf('\n');


               %更新阵列因子
               [Pattern_sense, ~]=plotAF(angle_sense,0,0);
               Pattern=(Pattern_sense+Pattern_comm)./2;
               %更新波束成形矢量
               pos=0:(NT-1);
               beamforming_sense=exp(-1j*pi*sind(angle_sense).*(pos.'));
               beamforming=0.5*sqrt(2)*beamforming_sense+0.5*sqrt(2)*beamforming_comm;
               %更新阵列因子
               h_pattern = Pattern(round((theta+180)*10+1));

               Ns_temp=Ns;
               NR_temp=NR;
               NT_temp=NT;
               P_symbol=0;
               %Ns个OFDM符号
               for ii = 1:Ns
                   P_symbol_ii=0;
                   %M个子载波
                   ii_temp=Ns_temp*(kk-1)+ii;
                   h_doppler_ii=h_doppler(ii,:);
                   for jj=1:M
    %                    Rx=Ar*diag(coefGen(ii,jj,r,v,theta,Pattern))*At.'...
    %                    *(beamforming.*reshape(Tx(jj,ii_temp,:),NT_temp,1)); 
    %                    Rx=Ar*diag(h_pathloss.*h_pattern.*coefGen(ii,jj,r,v,theta,Pattern))*At.'...
    %                    *(beamforming.*reshape(Tx(jj,ii_temp,:),NT_temp,1));
                       Rx=Ar*diag(h_pathloss.*h_pattern.*h_doppler_ii.*h_delay(jj,:))*At.'...
                        *(beamforming.*reshape(Tx(jj,ii_temp,:),NT_temp,1));
                       Rx_sensing(jj,ii,:)=Rx;
                       P_symbol_ii=P_symbol_ii+(Rx'*Rx);
                   end
                   P_symbol_ii=P_symbol_ii/(M*NR);
                   P_symbol=P_symbol+P_symbol_ii;
               end
               P_symbol=P_symbol/Ns;%计算平均每符号功率
               Rx_sensing=Rx_sensing/P_symbol;%归一化接收信号能量；


               % Add Noise  (此处设信自干扰比SSIR为0)
               Rx_sensing=Rx_sensing+sqrt(N0).*(randn(M,Ns,NR)+1j*randn(M,Ns,NR))/sqrt(2);

               % Received Symbol
               for iii=1:Ns
                 for jjj=1:M
                      RxData_sensing(jjj,iii)=beamforming_sense'*reshape(Rx_sensing(jjj,iii,:),NR_temp,1);
                 end
               end
               receiving_time=receiving_time+toc-temp;
        %% Process a direction
               temp=toc;
               % Covariance Matrix
                R=zeros(NR,NR);
                M_temp=M;
                for jjj=1:Ns
                    parfor iii=1:M
                       y=reshape(Rx_sensing(iii,jjj,:),NR_temp,1);
                       R=R+y*y';
                    end
                end
                R=R./(M*Ns);
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
                    [theta_samples, P_music_theta, L_estimate] = MUSIC_OFDMsensing(R, angle_start, angle_end);
                    if(L_estimate==0)
                        continue;
                    end
                else
                    fprintf('0 target detected!\n\n');
                    continue
                end

                %% Find MUSIC peaks 

                [ifNoPeak, peak_theta , peak_value]=Find_MUSIC_Peaks(theta_samples,P_music_theta,L_estimate);
                if(ifNoPeak==1)
                    fprintf('No MUSIC peak!\n\n');
                    continue;
                end
                L_estimate=length(peak_value);
                %作图 MUSIC伪谱
                if (ifPlotMusicSpectrum==1)
                    plotMusicSpectrum(theta_samples,P_music_theta,str);
                end

                %保存当前方向doa估计结果
                doa_result=vertcat(peak_theta , peak_value);
                doa_result=doa_result.';

                doa_time=doa_time+toc-temp;

                %%  Find Tau-Doppler peaks
                temp=toc;

                [range_value,velocity_value,pks,peaks_num,meanHeight]=Find_2DFFT_Peaks(P_TauDoppler,range,velocity);

                if (ifPlotPeriodogram==1)
                    plotPeriodogram(range,velocity,P_TauDoppler,str);
                end

                %保存当前方向r-v估计结果
                r_v_result=vertcat(range_value,velocity_value,pks,10*round(pks./meanHeight));
                r_v_result=r_v_result.';

                r_v_time=r_v_time+toc-temp;

                %% 保存当前结果

                doa_result=kron(doa_result,ones(peaks_num,1));
                r_v_result=repmat(r_v_result,L_estimate,1);

                result=horzcat(r_v_result, doa_result);

                if (flag==0)
                    flag=1;
                    Z=result;
                else
                    Z1=[Z; result];
                    Z=Z1;
                end
        end

        %% 去除冗余目标点
        temp=toc;
        %Zsort=sortrows(Z,[1,2,-3,-4]);%先按第1列距离r和第二列v升序，然后按相对峰值(第4列)和绝对峰值(第3列)降序排列
        Zsort=sortrows(Z,[-3,-4]);
        Zprun=Zsort(1,:);
        Target_Num=1;
        [Num0,~]=size(Zsort);
        delta_v=3*delta_v;
%         for i=2:Num0%遍历所有行
%             if ((Zprun(Target_Num,1)-delta_r<Zsort(i,1))&&(Zprun(Target_Num,1)+delta_r>Zsort(i,1))&&...
%                 (Zprun(Target_Num,2)-delta_v<Zsort(i,2))&&(Zprun(Target_Num,2)+delta_v>Zsort(i,2)) )
%                 if((Zprun(Target_Num,3)-3<Zsort(i,3))&&(Zprun(Target_Num,3)+3>Zsort(i,3))...
%                         &&(Zsort(i,4)>Zprun(Target_Num,4)))
%                     Zprun(Target_Num,:)=Zsort(i,:);
%                 end
%             else
%                 Target_Num=Target_Num+1;
%                 Zprun=[Zprun;Zsort(i,:)];
%             end
%         end
        for i=2:Num0%遍历所有行
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
        prunning_time=toc-temp;

        %% 
        figure(figureNumber);
        set(gcf,'Position',[488.2,205,801.6,389.6])
        subplot(1,2,1)
        polarplot(deg2rad(Z(:,5)),Z(:,1),'*');
        hold on
        polarplot(deg2rad(theta),r,'o');
        %legend("Estimated Targets","True Targets");
        hold off
        thetalim([-60,60]);

        title("Before Pruning");
        subplot(1,2,2)
        polarplot(deg2rad(Zprun(:,5)),Zprun(:,1),'*');
        hold on
        polarplot(deg2rad(theta),r,'o');
        hold off
        legend("Estimated Targets","True Targets");
        thetalim([-60,60]);

        title("After Pruning");
        str=['./fig/Figure ',num2str(figureNumber),'_MC',num2str(MC_time),'_SNR',num2str(SNR),'dB_Performance.png'];
        figureNumber=figureNumber+1;
        saveas(gcf, str);
        
        %% 计算误差
        Position=[r;v;theta];
        Position=Position';
        Position=sortrows(Position,1);
        e_DoA=zeros(Target_Num,length(theta));
        e_range=zeros(Target_Num,length(theta));
        e_velocity=zeros(Target_Num,length(theta));
        e_pos=zeros(Target_Num,length(theta));
        RMSE_pos=zeros(1,6);
        for ii=1:Target_Num
            for jj=1:length(theta)
                 e_pos(ii,jj)=(Position(jj,1)*cosd(Position(jj,3))-Zprun(ii,1)*cosd(Zprun(ii,5)))^2+(Position(jj,1)*sind(Position(jj,3))-Zprun(ii,1)*sind(Zprun(ii,5)))^2;
                 e_DoA(ii,jj)=(Position(jj,3)-Zprun(ii,5))^2;
                 e_range(ii,jj)=(Position(jj,1)-Zprun(ii,1))^2;
                 e_velocity(ii,jj)=(Position(jj,2)-Zprun(ii,2))^2;
            end
        end
        %检测概率
        success_detected=0;
        while (min(min(e_pos))<9)%当位置误差超过3m，视为检测失败
            min_value=min(min(e_pos));
            [row,col]=find(e_pos==min_value);

            RMSE_pos(Index_SNR)=RMSE_pos(Index_SNR)+min_value;%位置RMSE
            RMSE_DoA(Index_SNR)=RMSE_DoA(Index_SNR)+e_DoA(row,col);%DoA估计
            RMSE_range(Index_SNR)=RMSE_range(Index_SNR)+e_range(row,col);  %距离估计
            RMSE_velocity(Index_SNR)=RMSE_velocity(Index_SNR)+e_velocity(row,col);%速度估计

            e_pos(row,:)=[];
            e_pos(:,col)=[];
            success_detected= success_detected+1;
        end
            
        Pr_detect(Index_SNR)=Pr_detect(Index_SNR)+success_detected;
        Index_SNR=Index_SNR+1;
    end
 end
 %%
RMSE_pos=sqrt(RMSE_pos./Pr_detect);%位置RMSE
RMSE_DoA=sqrt(RMSE_DoA./Pr_detect);%DoA估计
RMSE_range=sqrt(RMSE_range./Pr_detect);  %距离估计
RMSE_velocity=sqrt(RMSE_velocity./Pr_detect);%速度估计
%%
Pr_detect= Pr_detect./(Total_MC_time*length(theta));%检测概率

 
total_time=toc;
fprintf("总耗时：%.2f s\n",total_time);
fprintf("准备时间占%.2f %%\t 发送接收耗时占 %.2f%%\t DoA估计耗时占%.2f%%\t 距离速度估计耗时占%2f%%\t"...
       ,100*prepare_time/total_time, 100*receiving_time/total_time, 100*doa_time/total_time,100*r_v_time/total_time);
%% 画图分析
SNR=-5:5:20;
figure(figureNumber);
subplot(1,3,1);
plot(SNR,RMSE_range,'-*');
title("(a)")
xlabel("SNR / dB");
ylabel("距离估计RMSE/m");
subplot(1,3,2);
plot(SNR,RMSE_velocity,'-*');
title("(b)")
xlabel("SNR / dB");
ylabel("速度估计RMSE/m·s^-^1");
subplot(1,3,3);
plot(SNR,RMSE_DoA,'-*');
title("(c)")
xlabel("SNR / dB");
ylabel("到达角估计RMSE/ °");
str=['./fig/Figure ',num2str(figureNumber),'_参数估计RMSE.png'];
figureNumber=figureNumber+1;
saveas(gcf, str);

figure(figureNumber);
subplot(1,2,1);
plot(SNR,RMSE_pos,'-*');
title("(d)")
xlabel("SNR / dB");
ylabel("位置估计RMSE/m");
subplot(1,2,2);
plot(SNR,Pr_detect,'-*');
title("(e)")
xlabel("SNR / dB");
ylabel("检测概率");
ylim([0 1]);
str=['./fig/Figure ',num2str(figureNumber),'_位置RMSE和检测概率.png'];
figureNumber=figureNumber+1;
saveas(gcf, str);



