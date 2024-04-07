%Find the RPY32, which can satisfy beta (under this condition)
%Using: Bisection

R1_MIN = 1e-6;
R1_MAX = 1e+3;
if iterBeta > 1
    R1_MAX = Rpy_rec(iterPP, iterBeta - 1);
end

biFlag=0; %Judge if found
while R1_MAX-R1_MIN>1e-9
    R1_mid = (R1_MIN+R1_MAX)/2;
    R1 = R1_mid;
    %[trying,iterS,iterS0,iterEPP,iterBeta,R1_mid]
	SO4_H2S_Steady_Full
    S32a=S32; S34a=S34; H32a=H32; H34a=H34; CHa=CH; Fea=Fe; % Profiles
    consumedOM=1-ResOrg;
    DPYa=CDeltaPY;       
    Pa=Percent;
    Rfluxa=Rflux;
    if(abs(Rfluxa-beta)<1e-5) %Find the RPY32, which satisfies beta
        biFlag=1;
        break
    end
    if(Rfluxa>beta)
        R1_MIN=R1_mid;
    else
        R1_MAX=R1_mid;
    end
end

