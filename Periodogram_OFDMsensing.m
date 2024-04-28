function [ifExist, delta_r,delta_v,range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData_sensing, TxData_origin,N0)
global M Ns c0 fc Ts delta_f Fp;

G=RxData_sensing./TxData_origin;
Kp=2^(ceil(log2(M)));% tau->range
Mp=2^(ceil(log2(Fp*Ns)));% doppler->velocity
delta_r=c0/(2*delta_f*Kp);
delta_v=c0/(2*fc*Ts*Mp);

W=kron(ones(1,Mp),barthannwin(M));
P_TauDoppler=ifft(fft(G,Mp,2).*W, Kp , 1);%P_TauDoppler(tau,doppler);

tau1=round(10*(2*delta_f*Kp)/c0);
tau2=round(90*(2*delta_f*Kp)/c0);
dopp1=1;
dopp2=round(60*(2*fc*Ts*Mp)/c0);
P_TauDoppler=P_TauDoppler(tau1:tau2,dopp1:dopp2);
P_TauDoppler=conj(P_TauDoppler).*P_TauDoppler;
P_TauDoppler = 10*log10(P_TauDoppler/N0);
fprintf("peak SNR=%.2f\n",max(max(P_TauDoppler)));
if(max(max(P_TauDoppler))<135)%80m,主瓣边缘=148.772
    ifExist=0;
else
    ifExist=1;
end
%velocity=1:Mp;
velocity=dopp1:dopp2;
velocity=velocity.*delta_v;
%velocity=kron(ones(Kp,1),velocity);
%range=1:Kp;
range=tau1:tau2;
range=range.*delta_r;
%range=kron(ones(1,Mp),range.');

end