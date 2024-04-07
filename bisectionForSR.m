%Find the RPY32, which can satisfy beta (under this condition)
%Using: Bisection
s_MIN = 0.01;
s_MAX = 0.20;
if iterPP > 1
    s_MAX = s_rec(iterPP - 1);
end
biFlag=0; %Judge if found
R1 = 1; % Just DSR, no pyrite formation thus no influence
while s_MAX-s_MIN>1e-6
    s_mid = (s_MIN+s_MAX)/2;
    s = s_mid;
    CH0 = PP_SET(iterPP) * 6 / 30 / s * 0.4 / 0.6 * 1000;
    SO4_H2S_Steady_Full
    S32a=S32; S34a=S34; H32a=H32; H34a=H34; CHa=CH; Fea=Fe; % Profiles
    consumedOM=1-ResOrg;
    DPYa=CDeltaPY;       
    Pa=Percent;
    Rfluxa=Rflux;
    if(abs(consumedOM-0.6)<1e-4) %Find the RPY32, which satisfies beta
        biFlag=1;
        break
    end
    if(consumedOM > 0.6)
        s_MIN=s_mid;
    else
        s_MAX=s_mid;
    end
end

