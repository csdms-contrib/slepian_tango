function [Z,G,H,K]=dyn_kernels_MCMC(deg,Bp,Visc,Te,planet,SurfBC)
%
% Calculates the flow in a viscous, self-gravitating sphere due to a 
% harmonic load, and returns a set of response kernels related to the
% resulting dynamic topography and geoid.  A good derivation can be found
% in Hager and Clayton 1989, and some terminology is borrowed from James et
% al. 2012. The effects of bending and membrane stresses are included using
% the formalism of Turcotte et al. 1981.
% Uses the function 'fralmanac.m' to retrieve planetary parameter values
% for Earth and 'ReturnPropagator.m' to calculate propagator matrices.
%
% INPUT:
% 
% deg      Desired spherical harmonic degree (can be a vector of degrees).
% Bp       Radius of the loading interface normalized by the planetary
%          radius (should be less than 1 and greater than the normalized
%          CMB radius).  Can also be a vector of loading radii, in which 
%          case the radii should be sorted in decreasing order.
% Visc     A two-columned matrix in which the first column contains
%          normalized interface radii and the second column contains
%          viscosities normalized by the reference viscosity. OR: 'PREM',
%          'isoviscous' [default]
%          (Note that none of the returned kernels depend on the absolute 
%          value of the reference viscosity, only relative viscosity
%          contrasts)
% Te       Elastic thickness, in units of km [default: 0]
% planet   String containing the name of the desired planet, from which 
%          physical constants are assumed, e.g. 'Earth', 'Mercury', 'Mars',
%          or 'Venus' [default: 'Earth']
% SurfBC   0 No-slip surface boundary condition [default]
%          1 Free-slip surface boundary condition
% 
% OUTPUT:
%
% Z        A non-dimensional (ndepths)x(ndegs) matrix giving the
%          degree-dependent ratios of geoid to topography for each inputted
%          loading depth (the "admittance kernel").
% G        A non-dimensional (ndepths)x(ndegs) matrix giving the "geoid
%          kernel".  This represents the ratio of geoid height for a given
%          mass load, relative to the static geoid anomaly produced by a
%          similar mass anomaly at the surface.
% H        A non-dimensional (ndepths)x(ndegs) matrix giving the normalized
%          surface displacements.
% K        A non-dimensional (ndepths)x(ndegs) "potential kernel" used by
%          some authors (related to the Geoid kernel g by an
%          upward-continuation factor.
%
% dyn_kernels('demo1') - Reproduces figures 9.21(a,d) of Hager and Clayton
%                        1989
% dyn_kernels('demo2') - Reproduces figures C1(b,c) of Katzman et al. 1998
% dyn_kernels('demo3') - Reproduces figures 6(a,b) of Herrick and Phillips
%                        1992
% 
% Created by Peter James (pjames@mit.edu) 6/10/12

% defval('deg',2:20)
% defval('Bp',1:-.01:.8)
% defval('Visc','isoviscous')
% defval('Te',0)
% defval('planet','Earth')
% defval('rho_m',3200) % Mantle density
% defval('drho_cm',6000) % Core-mantle density contrast
% defval('mu0',10^19)
% defval('SurfBC',0)
% defval('poisson',.25) % Poisson's ratio, dimensionless
% defval('E',10^11) % Young's modulus, units of Pa

if isempty(strfind(deg(:)','demo'))

if ischar(planet)
    if strcmp(planet,'Venus')
        R = 6051.880596*10^3;
        C = 3000*10^3;
        gr = 8.87;
    elseif strcmp(planet,'Mercury')
        R = 2440*10^3;
        C = 2040*10^3;
        gr = 3.7;
    elseif strcmp(planet,'Mars')
        R = 3396*10^3;
        C = 1500*10^3;
        gr = 3.71;
    elseif strcmp(planet,'Earth')
        R=6371000;
        C=3476000;
        gr=9.8100;
    else
        warning('Planet name not recognized')
    end
else
    warning('Please enter a string for the input parameter ''planet''');
end
% disp(['Radius of ',planet,' (km): ',num2str(R/1000,6)])
% disp(['CMB of ',planet,' (km): ',num2str(C/1000,6)])
% disp(['Gravity on ',planet,' (m/s^2): ',num2str(gr,3)])

% Normalized core radius
Cp=C/R;

% Convert to meters
Te=Te*1000;

if ischar(Visc)
    if strcmp(Visc,'PREM')
        Visc = ...
         [(R-670*10^3)/R, 10;...
         (R-400*10^3)/R, 1;...
         (R-60*10^3)/R, 1/30;...
         R/R, 100];
    elseif strcmp(Visc,'isoviscous')
        Visc = [R/R, 1];
    end
end

B=Bp*R;

Grav=6.6743*10^(-11); % Gravitational constant, units of m^3/kg/s^2
% Assume the gravitational acceleration throughout the interior is 
% approximately equal to the surface acceleration.
gc=gr;
gb=gr;
% Unitary driving load (kg/m^2)
load=1;

ndegs=length(deg);
ndepths=length(Bp);

Z={[]};
G={[]};
H={[]};
K={[]};

% L=1:lmax;
parfor (i=1:length(B),12)
    lcount=0;
    [OutValZ{i}, OutValG{i}, OutValH{i}, OutValK{i}]  = inner_loop(i,lcount,Z,G,H,K,R,C,B,Bp,Cp,Visc,Grav,gc,gb,load,Te,SurfBC,deg,gr);
end

OutValZ = OutValZ';
OutValG = OutValG';
OutValH = OutValH';
OutValK = OutValK';

Z=[];
G=[];
H=[];
K=[];
for jj = 1: ndepths 
Z(jj,:) = cell2mat(OutValZ{jj}(jj,:));
G(jj,:) = cell2mat(OutValG{jj}(jj,:));
H(jj,:) = cell2mat(OutValH{jj}(jj,:));
K(jj,:) = cell2mat(OutValK{jj}(jj,:));
end

elseif strcmp(deg,'demo1')
    R=6371;
    CMB=3476;
    depths=0:1:R-CMB;
    Bp = (R*ones(size(depths))-depths)/R;
    [~,G,H,~]=dyn_kernels([2,4,8],Bp,'isoviscous',0,'Earth',1);
    figure
    plot(depths,-H(:,1))
    hold on
    plot(depths,-H(:,2),':')
    plot(depths,-H(:,3),'--')

    view(90,90)
    hleg1 = legend('L = 2','L = 4','L = 8');
    title('L = 2,4,8 displacement kernels with depth','FontWeight','bold',...
        'FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Surface displacement kernels (dimensionless)','FontSize',10)
    ylim([0 1])
    xlim([0 R-CMB])
        hleg1 = legend('Te = 0','Te = 40','Te = 130');

    figure
    plot(depths,G(:,1))
    hold on
    plot(depths,G(:,2),':')
    plot(depths,G(:,3),'--')
    
    view(90,90)
    hleg2 = legend('L = 2','L = 4','L = 8');
    title('L = 2,4,8 geoid kernels with depth','FontWeight','bold',...
        'FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Geoid kernels (dimensionless)','FontSize',10)
    ylim([-.5 .5])
    xlim([0 R-CMB])

elseif strcmp(deg,'demo2')
    depths=0:10:1500;
    R=6371;
    Bp = (R*ones(size(depths))-depths)/R;
    Visc = ...
         [(R-670)/R, 10;...
         (R-400)/R, 1;...
         (R-60)/R, 1/30;...
         R/R, 100];
    [~,G1,H1,~]=dyn_kernels(26,Bp,Visc);
    [~,G2,H2,~]=dyn_kernels(26,Bp,Visc,40);
    [~,G3,H3,~]=dyn_kernels(26,Bp,Visc,130);
    figure
    plot(depths,G1)
    hold on
    plot(depths,G2,'.')
    plot(depths,G3,'--')

    view(90,90); legend;
    hleg1 = legend('Te = 0','Te = 40','Te = 130');
    title('L = 26 geoid kernels with depth','FontWeight','bold',...
        'FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Geoid kernel (dimensionless)','FontSize',10)
    ylim([-.3 .3])
    xlim([0 1500])

    figure
    plot(depths,H1)
    hold on
    plot(depths,H2,'.')
    plot(depths,H3,'--')
    
    view(90,90);
    hleg2 = legend('Te = 0','Te = 40','Te = 130');
    title('L = 26 displacement kernels with depth',...
        'FontWeight','bold','FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Surface displacement kernel (dimensionless)','FontSize',10)
    ylim([-1 1])
    xlim([0 1500])
    
elseif strcmp(deg,'demo3')
    deg=2:18;
    R=6051;
    Bp=(R-200)/R;
    [Z1,~,~,K1]=dyn_kernels(deg,Bp,'isoviscous',0,'Venus');
    Visc2 = [(R-800)/R, 1;...
             R/R, .1];
    [Z2,~,~,K2]=dyn_kernels(deg,Bp,Visc2,0,'Venus');
    Visc3 = [(R-800)/R, 1;...
             R/R, .01];
    [Z3,~,~,K3]=dyn_kernels(deg,Bp,Visc3,0,'Venus');
    Visc4 = [(R-200)/R, 1;...
             R/R, .1];
    [Z4,~,~,K4]=dyn_kernels(deg,Bp,Visc4,0,'Venus');
    Visc5 = [(R-200)/R, 1;...
             R/R, .01];
    [Z5,~,~,K5]=dyn_kernels(deg,Bp,Visc5,0,'Venus');

    figure
    plot(deg,K1)
    hold on
    plot(deg,K2,'linestyle','--','Color','r')
    plot(deg,K3,'linestyle','--','Color','g')
    plot(deg,K4,'linestyle','--','Color','c')
    plot(deg,K5,'linestyle','--','Color','m')
    hleg1 = legend('isoviscous','Model B','Model C','Model D','Model E');
    title('Potential kernels at 200 km depth',...
        'FontWeight','bold','FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('K','FontSize',10)
    ylim([-1 1])
    xlim([2 18])

    figure
    plot(deg,Z1)
    hold on
    plot(deg,Z2,'linestyle','--','Color','r')
    plot(deg,Z3,'linestyle','--','Color','g')
    plot(deg,Z4,'linestyle','--','Color','c')
    plot(deg,Z5,'linestyle','--','Color','m')
    hleg2 = legend('isoviscous','Model B','Model C','Model D','Model E');
    title('Admittance kernels at 200 km depth',...
        'FontWeight','bold','FontSize',12)
    xlabel('Loading depth (km)','FontSize',10)
    ylabel('Z','FontSize',10)
    ylim([-.1 .1])
    xlim([2 18])

end

end