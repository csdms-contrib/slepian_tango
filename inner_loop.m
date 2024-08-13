function [OutValZ, OutValG, OutValH, OutValK] = inner_loop(i,lcount,Z,G,H,K,R,C,B,Bp,Cp,Visc,Grav,gc,gb,load,Te,SurfBC,deg,gr)


%defval('deg',2:20)
%defval('Bp',1:-.01:.8)
%defval('Visc','isoviscous')
%defval('Te',0)
%defval('planet','Earth')
% defval('rho_m',3200) % Mantle density
% defval('drho_cm',6000) % Core-mantle density contrast
% defval('mu0',10^19)
% defval('SurfBC',0)
% defval('poisson',.25) % Poisson's ratio, dimensionless
% defval('E',10^11) % Young's modulus, units of Pa
rho_m=3200;
drho_cm=6000;
mu0=10^19;
poisson=0.25;
E=10^11;


for l=deg
    lcount=lcount+1;
    load_star=load/(R*rho_m);   
    
    [PBR,PCR]=ReturnPropagator(Bp(i),Cp,l,Visc);
    
    % Construct the system of equations
    
    % Nondimensionalization:
    % V{star} = V * mu_0 / (g_r*rho_m*R)
    % dr{star} = dr / R
    % Tau{star} = Tau / (g_r*rho_m*R)
    % load{star} = load / (rho_m*R)
    IX = 4*pi*Grav*R/(gr*(2*l+1));
    k1 = (-l^3*(l+1)^3+4*l^2*(l+1)^2)/(-l*(l+1)+1-poisson);
    k2 = (-l*(l+1)+2)/(-l*(l+1)+1-poisson);
    Elastic = E*Te^3/(12*(1-poisson^2)*rho_m*gr*(R-Te/2)^4) * k1 ...
        + E*Te/(rho_m*gr*(R-Te/2)^2) * k2;
    A=zeros(4,4);
    b=zeros(4,1);
    
    for ii=1:4
        di2=0; di3=0; di4=0;
        if ii==2; di2=1; elseif ii==3; di3=1; elseif ii==4; di4=1; end
        A(ii,1) = PCR(ii,2);
        A(ii,2) = PCR(ii,3)*drho_cm/rho_m*gc/gr*C/R ...
            - di3*IX*drho_cm*(C/R)^(l+2);
        A(ii,3) = di3*(1 + Elastic - IX*rho_m);
        if SurfBC==0
            A(ii,4) = di4;
        else
            A(ii,4) = di2;
        end
        b(ii) = (-PBR(ii,3)*B(i)/R*gb/gr ...
            + di3*IX*rho_m*(B(i)/R)^(l+2)) * load_star;
    end
    
    % Solution vector components: vcTH (lateral velocity at the core), 
    % perturbation at C, perturbation at R, Tau_theta at the surface

    % Improve the condition number of the system (
    ScaleFactor=NaN(4,1);
    for ii=1:4
        ScaleFactor(ii)=norm(A(:,ii))*10^15;
        A(:,ii)=A(:,ii)/ScaleFactor(ii);
    end
    
     Xvec=A\b;
%    Xvec=pinv(A)*b;
    
    Xvec=Xvec./ScaleFactor;
    dC_star=Xvec(2);
    dH_star=Xvec(3);

    % Re-dimensionalize:
    dH = dH_star * R;
    dC = dC_star * R;
    % Core velocities and surface velocities / shear stresses, if you  
    % desire them:
    vcTH = Xvec(1) * gr*rho_m*R^2 / mu0;
    if SurfBC==0
        tauTH = Xvec(4) * gr*rho_m*R;
    else
        vrTH = Xvec(4) * gr*rho_m*R^2 / mu0;
    end
    load = load_star * R*rho_m;

    UR = 4*pi*Grav/(2*l+1)*(R*rho_m*dH + C*(C/R)^(l+1)*drho_cm*dC + ...
        B(i)*(B(i)/R)^(l+1)*load);
    N=UR/gr;
    dUB = 4*pi*Grav*B(i)/(2*l+1)*(B(i)/R)^(l+1)*load;
    %dB = 4*pi*Grav*B(i)/(2*l+1)*load;
    dB = 4*pi*Grav*R/(2*l+1)*load;
    
    %if l==2
    %    [B(i) dB(1) dUB(1)]
    %end
    
    Z(i,lcount)={UR/(gr*dH)};  
    G(i,lcount)={UR/dB};
    H(i,lcount)={dH*rho_m/load};
    K(i,lcount)={UR/dUB};
end
OutValZ = Z;
OutValG = G;
OutValH = H;
OutValK = K;