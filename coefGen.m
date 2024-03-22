function Beta = coefGen(range,velocity)
    global c0 lambda M N Ts delta_f fc
    L=length(range);
    Beta=zeros(M,N,L); %M行N列，M个子载波N个OFDM符号
    delay = 2 * range / c0; %各目标的往返时延。
    h_gain_phase = exp(1j*2*pi*rand(size(delay)));    % 各目标随机的信道相位响应
    h_gain_amp = c0*1./((4*pi)^(1.5)*(fc)*(range).^2); %路径损耗
    doppler = 2*velocity/lambda;    %各目标的多普勒频移，这里的速度是在LoS方向的分量
    for ii=1:M
        for jj=1:N
            Beta(ii,jj,1:L)=h_gain_amp.*h_gain_phase.*exp(1j*2*pi*(jj-1)*Ts.*doppler).*exp(-1j*2*pi*(ii-1)*delta_f.*delay);
        end
    end
end