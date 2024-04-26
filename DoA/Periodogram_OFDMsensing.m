function [range, velocity, P_TauDoppler]=Periodogram_OFDMsensing(RxData, TxData)
global M N c0 fc Ts delta_f Fp;
 
G=RxData./TxData;
Kp=2^(ceil(log2(M)));% tau->range
Mp=2^(ceil(log2(Fp*N)));% doppler->velocity
P_TauDoppler=ifft( fft(G,Mp,2) , Kp , 1);%P_TauDoppler(tau,doppler);
P_TauDoppler=conj(P_TauDoppler).*P_TauDoppler;
velocity=1:Mp;
velocity=velocity.*c0./(2*fc*Ts*Mp);
%velocity=kron(ones(Kp,1),velocity);
range=1:Kp;
range=range.*c0./(2*delta_f*Kp);
%range=kron(ones(1,Mp),range.');
 
end
