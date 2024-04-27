%% MUSIC_OFDMsensing
% CIM: Input channel information matrix(pre-processed by known transmitted symbols)
% k:target number
% Author: Yunbo HU(SIMIT, UCAS)
% GitHub: https://github.com/edenhu1111
function [theta, P_music_theta, L] = MUSIC_OFDMsensing(R, angle_start,angle_end,N0)
    global  NR
    %% DoA estimation
    Iv=eye(NR,NR);
    Iv=fliplr(Iv);
    R=R+Iv*conj(R)*Iv;%Toeplitz化
    [V,D]=eig(R);
    [d,ind_D] = sort(diag(D),'descend');
    L=MDL(d);
    if(L==0)
         fprintf("According to MDL criterion, L=0\n");
        return;
    end
    U_n=V(:,ind_D(L+1:end));
    theta=linspace(angle_start,angle_end,min(100,round((angle_end-angle_start)/0.1)));
    A=steeringGen(theta, NR);
    P_music_theta=zeros(1,length(theta));
    for ii=1:length(theta)
        P_music_theta(ii)=1/(A(:,ii)'*(U_n*U_n')*A(:,ii));
    end
    P_music_theta=abs(P_music_theta);
    P_music_theta=P_music_theta/max(P_music_theta);
    P_music_theta=10*log10(P_music_theta);
    if(range(P_music_theta)<0.05)
        fprintf("No MUSIC peak!\n");
        L=0;
    end
end

%% find the number of targets using MDL criterion
function L=MDL(S) %S=sort(diag(D),'decend')
global NR M Ns
    s=0;
    A=S.^(1/NR-s);
    mdl_min=-((NR-s)*M*Ns)*log(prod(A(s+1:NR))/(1/(NR-s)*sum(S(s+1:NR))))+0.5*s*(2*NR-s)*log(M*Ns);
    L=0;
    for s=1:length(S)-1
        A=S.^(1/NR-s);
        mdl=-((NR-s)*M*Ns)*(log(prod(A(s+1:NR)))-log((1/(NR-s)*sum(S(s+1:NR)))))+0.5*s*(2*NR-s)*log(log(M*Ns));
        if(mdl<mdl_min)
            mdl_min=mdl;
            L=s;
        end
    end
end