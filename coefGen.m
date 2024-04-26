%第n个OFDM符号，第m个子载波
function Beta = coefGen(n,m,range,velocity,theta)
    global c0 lambda Ts delta_f fc NT
    L=length(range);
    Beta=zeros(1,L); %M行N列，M个子载波N个OFDM符号
    delay = 2 * range / c0; %各目标的往返时延。
   % h_phase = exp(1j*pi*randn(size(delay)));    % 各目标随机的信道相位响应
    h_pathloss = c0*1./((4*pi)^(1.5)*(fc)*(range).^2); %路径损耗
    h_pattern = abs((sin(NT*0.5*pi*sin(theta)))./(sin(0.5*pi*sin(theta))));
    h_pattern = h_pattern/max(h_pattern);%归一化阵列因子
    doppler = 2*velocity/lambda;    %各目标的多普勒频移，这里的速度是在LoS方向的分量
            %Beta(ii,jj,1:L)=h_pathloss.*h_phase.*exp(1j*2*pi*(jj-1)*Ts.*doppler).*exp(-1j*2*pi*(ii-1)*delta_f.*delay);
            Beta=h_pathloss.*h_pattern.*exp(1j*2*pi*(n-1)*Ts.*doppler).*exp(-1j*2*pi*(m-1)*delta_f.*delay);
end