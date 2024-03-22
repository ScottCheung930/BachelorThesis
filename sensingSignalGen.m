%% Sensing Channel Generation
% Description: Generating a simulated time shifted and doppler modulaed
% signal.
% TxSignal_cp: transmit signal
%
% range: target range. When it is a row vector, it indicated that there are
% more than one target;
%
% velocity:target relative velocity. Its size should be the same as "range"
%
% SNR: Signal-to-noise Ratio
%
% Author: Yunbo HU (SIMIT, UCAS)
% GitHub: https://github.com/edenhu1111
%% Code
function RxSignal = sensingSignalGen(TxSignal_cp,range,velocity,SNR)
    global c0 lambda M delta_f
    delay = 2 * range / c0; %各目标的往返时延。
    h_gain = exp(1j*2*pi*rand(size(delay)));    %各目标随机的信道相位响应
    doppler = 2*velocity/lambda;    %各目标的多普勒频移，这里的速度是在LoS方向的分量
%     max_delay = max(delay,[],'all');
%     RxSignal = zeros(size(TxSignal_cp,1) + max_delay,1);
    RxSignal = zeros(size(TxSignal_cp,1)    ,1);
    d = zeros(size(TxSignal_cp)); %多普勒项，列向量
    for p = 1:length(delay)% 分别计算每个目标产生的回波
        ii = 0:length(d)-1;% ii ，接收信号下标，行向量
        % 多普勒引起的相移：
        % 对于第m个OFDM符号，exp(j*2*pi*m*Ts*doppler)
        d = exp(1j*2*pi*doppler(p)*ii'/(delta_f*M)); %d
        % 往返时延引起的相移：
        % 对于第k个子载波，exp(-j*2*pi*k*delta_f*delay)
        tau = exp(-1j*2*pi*delay(p)*delta_f*M/length(d)*ii');
%         RxSignal = RxSignal + h_gain(p) * ...
%             [zeros(delay(p),1);...
%             TxSignal_cp .* d...
%             ;zeros(max_delay-delay(p),1)];
        RxSignal = RxSignal + h_gain(p).*ifft(fft(TxSignal_cp.*d).*tau);
        %结果相当于h_gain*TxSignal_cp*(cycconv(d,tau))
    end
    RxSignal = RxSignal + 10^(-SNR/10)*(randn(size(RxSignal)) + 1j*randn(size(RxSignal)))/sqrt(2);
end

