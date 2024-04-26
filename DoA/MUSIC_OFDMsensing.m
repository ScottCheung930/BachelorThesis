%% MUSIC for OFDM sensing
function [P_music_theta,L]= MUSIC_OFDMsensing(R)
    global  NR
    %% DoA estimation
    [V,D]=eig(R);
    [d,ind_D] = sort(diag(D),'descend');
    L=MDL(d);
    U_n=V(:,ind_D(L+1:end));
    theta=linspace(-30,30);
    A=steeringGen(theta, NR);
    P_music_theta=zeros(1,length(theta));
    for ii=1:length(theta)
        P_music_theta(ii)=1./(A(:,ii)'*(U_n*U_n')*A(:,ii));
    end
    
end

%% find the number of targets using MDL criterion
function L=MDL(S) %S=sort(diag(D),'decend')
global NR M N
    s=0;
    A=S.^(1/NR-s);
    mdl_min=-((NR-s)*M*N)...
             *log(prod(A(s+1:NR))/(1/(NR-s)*sum(S(s+1:NR))))...
             +0.5*s*(2*NR-s)*log(M*N);
    L=0;
    for s=1:length(S)-1
        A=S.^(1/NR-s);
        mdl=-((NR-s)*M*N)...
             *log(prod(A(s+1:NR))/(1/(NR-s)*sum(S(s+1:NR))))...
             +0.5*s*(2*NR-s)*log(M*N);
        if(mdl<mdl_min)
            mdl_min=mdl;
            L=s;
        end
    end
end