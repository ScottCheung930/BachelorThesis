%计算并画出阵列天线方向图，计算-10dB主瓣宽度亦即单次波束扫描宽度(平均)
%%
%function [Pattern, Delta_theta]=plotPattern(phi,ifPlot,ifBeamWidth)%phi为天线最大辐射方向(角度制)
    NT=20;phi_1=50;phi_2=-10;
    theta=-180:0.1:179.9;
    theta_radians=deg2rad(theta);
    phi_1_radians=deg2rad(phi_1);
    Pattern_1=abs((sin(NT*0.5*pi*(sin(theta_radians)-sin(phi_1_radians))))./(sin(0.5*pi*(sin(theta_radians)-sin(phi_1_radians)))));
    for theta=-180:0.1:179.9%处理阵列因子中的最大值点
        if ((theta-phi_1)==0 || theta+phi_1 == 180 || theta+phi_1 == -180)
            Pattern_1((phi_1+180)*10+1)=NT;
        end
    end
    
    phi_2_radians=deg2rad(phi_2);
    Pattern_2=abs((sin(NT*0.5*pi*(sin(theta_radians)-sin(phi_2_radians))))./(sin(0.5*pi*(sin(theta_radians)-sin(phi_2_radians)))));
    for theta=-180:0.1:179.9%处理阵列因子中的最大值点
        if ((theta-phi_2)==0 || theta+phi_2 == 180 || theta+phi_2 == -180)
            Pattern_2((phi_2+180)*10+1)=NT;
        end
    end
     Pattern_1= Pattern_1/NT; %归一化
     Pattern_2= Pattern_2/NT; %归一化
        figure(1);
        polarplot(theta_radians, Pattern_1);
        thetalim([-90,90]);
        hold on
        polarplot(theta_radians, Pattern_2);
        thetalim([-90,90]);
        hold off;
        str=['./ArrayFactor__.png'];
        saveas(gcf,str);
        close(gcf);
    
    %Delta_theta=Delta_theta*1.5;
    %在扫描过程中，波束宽度会逐渐增大
    %增大的比例近似为1/cos(方向角)，当方向角为60°时，增大至2倍
%end