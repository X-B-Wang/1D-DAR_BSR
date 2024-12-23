%%%This is the portion of program that actually computes the coupled DAR system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Do not further adjust code below here unless changing dynamics%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set up grids
dz=zL/N; %grid size; 
zz=linspace(0,zL,N+1);zz=zz'; %z grid, range 0 -- zL

%Tied physical properties
r0=R0*(1-DeltaS_SET/1000); %(RBSR34, tied to RBSR32)
r1=R1; %(RPY34, tied to RPY32)
% (S34/S32)std=0.0450045
% ((S34/S32)i-(S34/S32)std)/(S34/S32)std*1000=d34Si
% (S34/S32)i/(S34/S32)std-1=d34Si/1000
% (S34/S32)i=(d34Si/1000+1)*(S34/S32)std
% S34 : S32 = 0.0450045(d34Si/1000+1) : 1 = 0.0450045(d34Si+1000) : 1000 =
% 0.0450045*d34Si+45.0045 : 1000
S32_SET=1000;
S34_SET=0.0450045*d34Sso4+45.0045;
S032=S0*S32_SET/(S32_SET+S34_SET);
S034=S0*S34_SET/(S32_SET+S34_SET);
%S032=S0*95/100.101; S034=S0*5.101/100.101; %(Scale into actual [SO4]32 and [SO4]34 B.C.) 

%Initialize CH2O and Fe profiles for iterative process, if not iterating, can be interpretted as unlimited
CH=CH0+0*zz;
Fe=Fe0+0*zz;

%CH2O-SO4 reactions
for ii=1:ni
    CHO=CH; %intermediate variable check convergence
  %Set up differentiation operator matrix for S32
    Gu=[1:1:N-1 2:1:N 1:1:N-1 N];
    Gv=[2:1:N 1:1:N-1 1:1:N-1 N];
    GS=[ones(1,N-1)*(DS/dz^2-s/2/dz) ones(1,N-1)*(DS/dz^2+s/2/dz) -2*DS/dz^2-R0*CH(2:N)' -DS/dz^2-s/2/dz-R0*CH(N+1)];
    G=sparse(Gu,Gv,GS);
  %Set up forcing for S32
    RHS=0*zz(1:N);
    RHS(1)=-S032*(DS/dz^2+s/2/dz);
  %[SO4]32 profile based on [CH2O] guess
    
  %Use Thomas method to solve for tridiagonal system
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
        ys32=y;
        S32=[S032;y;];
    

   
  %Set up differentiation operator matrix for S34
    Gu=[1:1:N-1 2:1:N 1:1:N-1 N];
    Gv=[2:1:N 1:1:N-1 1:1:N-1 N];
    GS=[ones(1,N-1)*(DS/dz^2-s/2/dz) ones(1,N-1)*(DS/dz^2+s/2/dz) -2*DS/dz^2-r0*CH(2:N)' -DS/dz^2-s/2/dz-r0*CH(N+1)];
    G=sparse(Gu,Gv,GS);
  %Set up forcing for S34
    RHS=0*zz(1:N)';
    RHS(1)=-S034*(DS/dz^2+s/2/dz);
  %[SO4]34 profile based on [CH2O] guess
  
    %Use Thomas method to solve for tridiagonal system
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
    ys34=y;
    S34=[S034;y;];
   
  %Integrate [SO4]32 and [SO4]34 with 2nd order accuracy
    S32I=zeros(1,N+1); S34I=zeros(1,N+1);
    for i=2:N+1
        S32I(i)=((S32(1)+S32(i))/2+sum(S32(2:i-1)))*dz;
        S34I(i)=((S34(1)+S34(i))/2+sum(S34(2:i-1)))*dz;
    end
  
    if ni~=1 %ni~=1 means iterate CH2O for depletion
        for i=1:N+1
            CH(i)=CH0*exp(-2*(R0/s*S32I(i)+r0/s*S34I(i)));%Update CH2O using analytic expression
        end
    end
  %Check CH2O convergence
    if (mean((CHO(:)-CH(:)).^2)<1e-16)
        break
    end
end
if (ii==ni && ni~=1)
    ii
end


SSA=R0*S32(:).*CH(:); %[CH2O]*[SO4]32 as available for [H2S]32 reaction
ssa=r0*S34(:).*CH(:); %[CH2O]*[SO4]34 as available for [H2S]34 reaction

%Fe=Fe0*CH/CH0;
for ii=1:nni
    FeO=Fe;%Intermediate variable check convergence
  %Set up differentiation operator matrix for S32
    Hu=[1:1:N-1 2:1:N 1:1:N-1 N];
    Hv=[2:1:N 1:1:N-1 1:1:N-1 N];
    HS=[ones(1,N-1)*(DH2S/dz^2-s/2/dz) ones(1,N-1)*(DH2S/dz^2+s/2/dz) -2*DH2S/dz^2-R1*Fe(2:N)' -1*DH2S/dz^2-s/2/dz-R1*Fe(N+1)];
    H=sparse(Hu,Hv,HS);
  %Forcing [H2S]32
    RHSH=-SSA(2:end);
    
    %Use Thomas method to solve for tridiagonal system
    f = RHSH;
    n = length(f);
    v = zeros(n,1);   
    y = v;
    w = H(1,1);
    y(1) = f(1)/w;
    for j=2:n
        v(j-1) = H(j-1,j)/w;
        w = H(j,j) - H(j,j-1)*v(j-1);
        y(j) = ( f(j) - H(j,j-1)*y(j-1) )/w;
    end
    for j=n-1:-1:1
        y(j) = y(j) - v(j)*y(j+1);
    end
    yh32=y;
    H32=[0;y;];    
    
  %Set up differentiation operator matrix for S34
    Hu=[1:1:N-1 2:1:N 1:1:N-1 N];
    Hv=[2:1:N 1:1:N-1 1:1:N-1 N];
    HS=[ones(1,N-1)*(DH2S/dz^2-s/2/dz) ones(1,N-1)*(DH2S/dz^2+s/2/dz) -2*DH2S/dz^2-r1*Fe(2:N)' -1*DH2S/dz^2-s/2/dz-r1*Fe(N+1)];
    H=sparse(Hu,Hv,HS);
  %Forcing [H2S]34
    RHSH=-ssa(2:end);
    
    %Use Thomas method to solve for tridiagonal system
    f = RHSH;
    n = length(f);
    v = zeros(n,1);   
    y = v;
    w = H(1,1);
    y(1) = f(1)/w;
    for j=2:n
        v(j-1) = H(j-1,j)/w;
        w = H(j,j) - H(j,j-1)*v(j-1);
        y(j) = ( f(j) - H(j,j-1)*y(j-1) )/w;
    end
    for j=n-1:-1:1
        y(j) = y(j) - v(j)*y(j+1);
    end
    H34=[0;y;];    
    yh34=y;
    
  %Integrate [H2S]32 and [H2S]34 with 2nd order accuracy
    H32I=zeros(1,N+1);H34I=zeros(1,N+1);
    for i=2:N+1
        H32I(i)=((H32(1)+H32(i))/2+sum(H32(2:i-1)))*dz;
        H34I(i)=((H34(1)+H34(i))/2+sum(H34(2:i-1)))*dz;
    end   
  %Update [Fe] profile
  
    if nni==1 %nni~=1 means iterate Fe for depletion
    else
        for i=1:N+1
            Fe(i)=Fe0*exp(-0.5*R1/s*H32I(i)-0.5*r1/s*H34I(i));%Update Fe
        end
    end
  %Check Fe convergence
    if (max(abs(FeO(:)./Fe(:)-1))<1e-16)
        break
    end
end
if (ii==nni && nni~=1)
    ii
end

% N layers, so N+1 grids, using Simpson for integration
%Compute ratio of diffusive flux vs reaction
gr=-(H32(1)-H32(2)+H34(1)-H34(2))*DH2S/dz;
%dr=-s/2*(yh32(1)-2*yh32(2000));
PY=((H32(1)+H34(1))*Fe(1) ...
    +4*sum((H32(2:2:end-1)+H34(2:2:end-1)).*Fe(2:2:end-1)) ...
    +2*sum((H32(3:2:end-2)+H34(3:2:end-2)).*Fe(3:2:end-2)) ...
    +(H32(end)+H34(end))*Fe(end))*R1*dz/6;
COST=((S32(1)*R0+S34(1)*r0)*CH(1) ...
    +4*sum((S32(2:2:end-1)*R0+S34(2:2:end-1)*r0).*CH(2:2:end-1)) ...
    +2*sum((S32(3:2:end-2)*R0+S34(3:2:end-2)*r0).*CH(3:2:end-2)) ...
    +(S32(end)*R0+S34(end)*r0)/2*dz*CH(end))*dz/6;
Rflux=gr/COST; %Flux at the boundaries
Rdep=PY/COST;
%Rs=dr/COST %Check if H2S lower boundary is proper. =1 should be ok

%weighted average delta H2S
delta=log(H34(:)./H32(:)/0.0450045)*1000;
ind=H32==0|H34==0;delta(ind)=0;
ind=find(zz==zL);
PY32=(H32(1)*Fe(1) ...
    +4*sum(H32(2:2:end-1).*Fe(2:2:end-1)) ...
    +2*sum(H32(3:2:end-2).*Fe(3:2:end-2)) ...
    +H32(end)*Fe(end))*R1*dz/6;
PY34=(H34(1)*Fe(1) ...
    +4*sum(H34(2:2:end-1).*Fe(2:2:end-1)) ...
    +2*sum(H34(3:2:end-2).*Fe(3:2:end-2)) ...
    +H34(end)*Fe(end))*r1*dz/6;
CDeltaPY=log((PY34/PY32)/0.0450045)*1000;

%mass ratio
ZPY=120;
WPY=ZPY*0.5*((H32(1)*R1+H34(1)*r1).*Fe(1) ...
    +4*sum((H32(2:2:end-1)*R1+H34(2:2:end-1)*r1).*Fe(2:2:end-1)) ...
    +2*sum((H32(3:2:end-2)*R1+H34(3:2:end-2)*r1).*Fe(3:2:end-2)) ...
    +(H32(end)*R1+H34(end)*r1).*Fe(end))*dz/6;
Percent=WPY/800000*.6/s;

%Residual Organic Matter ratio
ResOrg=CH(end)/CH0;
TOCa=CH(end)*2.5/1000*12*0.6/0.4/2000;
% mmol/L * 1 * 0.001mol/mmol * 12g/mol = g/L =kg/m3



