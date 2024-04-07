close all; clc; 
clear all; format long

DS = 4.24e-10*365*24*3600*1000;  %(Diffusivity of [SO4] m^2/kyear, same for S32 and S34)
DH2S = 6.14e-10*365*24*3600*1000;  %(Diffusivity of [H2S] m^2/kyear, same for S32 and S34)
% 0.0033
R0 = 0.0123; %(RBSR32, 1/kyear/(M/m^3) tune for depth, then fix, 0.0123, balanced depletion)

DeltaS_SET = 46; %RBSR34 vs. RBSR32
d34Sso4 = 35; %d34Sso4_seawater
Fe0 = 300; %Boundary condition for Fe, mM/L
maxFe = Fe0; %Max Fe0

%Numerical control
zL=100;  %Height below which [H2S] approach 0, must be large enough!
N=zL*20; %Fine grid to resolve boundary layer, 200 per meter is better, usually 2000 for saving time(also my computer)(T_T)
ni=1000; %number of iterations for SO4 equations, ni=1 means [CH2O] unlimited
nni=3000; %number of iterations for H2S equations, nni=1 means [Fe] unlimited

s = 0.05; %m/kyr, sedimentation rate
S0 = 2; % mMol/L, seawater sulfate concentration
PP_SET = [0.4 0.5 0.6 0.7 0.8 0.9 1];
beta_SET = [0.1 0.25 0.5 0.75 0.9];
d34Spy_X = zeros(length(PP_SET),length(beta_SET));
PYwt_Y = zeros(length(PP_SET),length(beta_SET));
Rec = zeros(length(PP_SET)*length(beta_SET), 8);
cnt = 0;
s_rec = zeros(length(PP_SET), 1);
Rpy_rec = 1e-6 + zeros(length(PP_SET), length(beta_SET));
for iterPP = 1:length(PP_SET) %Enumerate EPP, which decide CH0
    bisectionForSR
    s_rec(iterPP) = s;
    if biFlag == 0
        continue;
    end
    for iterBeta = 1:length(beta_SET) %Enumerate beta
        beta = beta_SET(iterBeta);
        bisectionForRPY32
        Rpy_rec(iterPP, iterBeta) = (R1_MIN+R1_MAX)/2;
        d34Spy_X(iterPP, iterBeta) = DPYa;
        PYwt_Y(iterPP, iterBeta) = Pa * 100;
        cnt = cnt + 1;
        Rec(cnt, :) = [S0, PP_SET(iterPP), Fe0, s, (R1_MIN+R1_MAX)/2, DPYa, Pa * 100, Rfluxa];
        [iterPP iterBeta]
    end
end

xlswrite('Result.xlsx', [["SO4", "PP", "Fe0", "s", "Rpy", "d34Spy", "[pyrite]", "beta"]; Rec]);
