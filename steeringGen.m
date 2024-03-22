%考虑ULA阵列天线之间半波长间隔
%生成阵列流形
function A=steeringGen(theta,N)
    L=length(theta);
    A=zeros(N,L);
    for ii=1:length(theta)
        theta_ii=theta(1,ii);
        for jj=1:N
            A(jj,ii)=exp((jj-1)*1j*pi*sin(theta_ii));
        end
    end
end