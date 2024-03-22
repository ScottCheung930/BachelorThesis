%% MUSIC for OFDM sensing
% CIM: Input channel information matrix(pre-processed by known transmitted symbols)
% k:target number
% Author: Yunbo HU(SIMIT, UCAS)
% GitHub: https://github.com/edenhu1111
function [P_music_theta] = MUSICforOFDMsensing(R)
global  NR
%% DoA estimation
[V,D]=eig(R);
[d,ind_D] = sort(diag(D),'descend');
L=MDL(d);
U_n=V(:,ind_D(L+1:end));
theta=linspace(-pi/3,pi/3);
A=steeringGen(theta, NR);
P_music_theta=zeros(1,length(theta));
for ii=1:length(theta)
    P_music_theta(ii)=1/(A(:,ii)'*(U_n*U_n')*A(:,ii));
end



%% range estimation
% [V,D]=eig(R);
% 
% [~,ind_D] = sort(diag(D),'descend');
% U_n = V(:,ind_D(k+1:end));
% 
% delay = linspace(0,2*100/c0*delta_f,M);
% ii = 0:M-1;
% A = exp(-1j*2*pi*kron(ii',delay));
% P_music_range = zeros(size(A,2),1);
% for jj = 1:size(A,2)
%     P_music_range(jj) = 1/(A(:,jj)'*(U_n*U_n')*A(:,jj));
% end


%% Velocity Estimation
% R_dop = CIM.'*conj(CIM)/M;
% [V,D]=eig(R_dop);
% [~,ind_D] = sort(diag(D),'descend');
% U_n = V(:,ind_D(k+1:end));
% 
% doppler = linspace(0,2*100/lambda*Ts,M);
% ii = 0:N-1;
% A = exp(1j*2*pi*kron(ii',doppler));
% P_music_velo = zeros(size(A,2),1);
% for jj = 1:size(A,2)
%     P_music_velo(jj) = 1/(A(:,jj)'*(U_n*U_n')*A(:,jj));
% end


end

%% find the number of targets using MDL criterion
function L=MDL(S) %S=sort(diag(D),'decend')
global NR M N
    s=0;
    A=S.^(1/NR-s);
    mdl_min=-((NR-s)*M*N)*log(prod(A(s+1:NR))/(1/(NR-s)*sum(S(s+1:NR))))+0.5*s*(2*NR-s)*log(M*N);
    L=0;
    for s=1:length(S)-1
        A=S.^(1/NR-s);
        mdl=-((NR-s)*M*N)*log(prod(A(s+1:NR))/(1/(NR-s)*sum(S(s+1:NR))))+0.5*s*(2*NR-s)*log(M*N);
        if(mdl<mdl_min)
            mdl_min=mdl;
            L=s;
        end
    end
end