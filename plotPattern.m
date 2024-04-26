%计算并画出阵列天线方向图，计算-10dB主瓣宽度亦即单次波束扫描宽度
%%
function [Pattern, Delta_theta]=plotPattern()
    global NT;
    %NT=10;
    theta=-180:0.1:179.9;
    theta_radians=deg2rad(theta);
    Pattern=abs((sin(NT*0.5*pi*sin(theta_radians)))./(sin(0.5*pi*sin(theta_radians))));
    Pattern(1801)=NT;
    %画方向图
    figure(1);
    polarplot(theta_radians, Pattern);
    text=["ULA Pattern with NT=NR="+num2str(NT)];
    title(text);
    
    %找-10dB波束宽度（边射）
    maximum=Pattern(1801);
    for i=1:1800
        if(20*log10(maximum/Pattern(1800+i))>=10)
            Delta_theta=2*i*0.1;
            break;
        end
    end
    
    Pattern=Pattern./maximum;
    Delta_theta=Delta_theta*1.5;
    %在扫描过程中，波束宽度会逐渐增大
    %增大的比例为1/cos(方向角)，当方向角为60°时，增大至2倍
end