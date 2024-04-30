%第n个OFDM符号，第m个子载波
%theta为角度制
function Beta = coefGen(n,m,range,velocity,theta,Pattern)
    global c0 lambda Ts delta_f fc
    %L=length(range);
    delay = 2 * range / c0; %各目标的往返时延。
   % h_phase = exp(1j*pi*randn(size(delay)));    % 各目标随机的信道相位响应
    h_pathloss = c0./((fc).*(range).^2); %路径损耗
    h_pattern = Pattern(round((theta+180)*10+1));
    doppler = 2*velocity/lambda;    %各目标的多普勒频移，这里的速度是在LoS方向的分量
    %Beta=h_pathloss.*exp(1j*2*pi*(n-1)*Ts.*doppler).*exp(-1j*2*pi*(m-1)*delta_f.*delay);
    Beta=h_pathloss.*h_pattern.*exp(1j*2*pi*(n-1)*Ts.*doppler).*exp(-1j*2*pi*(m-1)*delta_f.*delay);
end