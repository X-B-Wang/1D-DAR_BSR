zL=100;  %Height below which [H2S] approach 0, must be large enough!
N=zL*20; %Fine grid to resolve boundary layer, 200 per meter is better, usually 2000 for saving time(also my computer)(T_T)
ni=1000; %number of iterations for SO4 equations, ni=1 means [CH2O] unlimited
nni=3000; %number of iterations for H2S equations, nni=1 means [Fe] unlimited
dz=zL/N; %grid size; 
zz=linspace(0,zL,N+1);zz=zz';


welling = 0;
CH0 = 10; %uM
oxygen0 = 300; % uM
Rox = 0.03;
sink = 1*1000; % m/kyear

DO2 = 1.46e-9*365*24*3600*1000; % m^2/kyear


oxygen = oxygen0 + 0 * zz;
CH = CH0 + 0 * zz;

for iter = 1:1000
    oxygenp = oxygen;
    
    Gu=[1:1:N-1 2:1:N 1:1:N-1 N];
    Gv=[2:1:N 1:1:N-1 1:1:N-1 N];
    GS=[ones(1,N-1)*(DO2/dz^2-welling/2/dz) ones(1,N-1)*(DO2/dz^2+welling/2/dz) -2*DO2/dz^2-Rox*CH(2:N)' -DO2/dz^2-welling/2/dz-Rox*CH(N+1)];
    G=sparse(Gu,Gv,GS);
    
    RHS=0*zz(1:N);
    RHS(1)=-oxygen0*(DO2/dz^2+welling/2/dz);

    f = RHS;
    n = length(f);
    v = zeros(n,1);   
    y = v;
    w = G(1,1);
    y(1) = f(1)/w;
    for j=2:n
        v(j-1) = G(j-1,j)/w;
        w = G(j,j) - G(j,j-1)*v(j-1);
        y(j) = ( f(j) - G(j,j-1)*y(j-1) )/w;
    end
    for j=n-1:-1:1
        y(j) = y(j) - v(j)*y(j+1);
    end

    oxygen = [oxygen0; y];
    
    Ioxygen = zeros(1, N+1);
    for i = 2:N+1
        Ioxygen(i) = ((oxygen(1) + oxygen(i)) / 2 + sum(oxygen(2:i-1))) * dz;
    end
    
    CHp = CH;
    for i = 1:N+1
        CH(i) = CH0 * exp(-Rox / (sink + welling) * Ioxygen(i));
    end
    if max(abs(CH - CHp)) < 1e-3 && max(abs(oxygenp - oxygen)) < 1e-3
        break;
    end
end

subplot(1, 2, 1);
plot(oxygen, zz); set(gca, 'yDir', 'reverse');

subplot(1, 2, 2);
plot(CH, zz); set(gca, 'yDir', 'reverse');

xlswrite('Results.xlsx', [zz(1:20:end) oxygen(1:20:end) CH(1:20:end)]);