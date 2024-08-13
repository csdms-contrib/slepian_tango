function varargout=geoidkernel_MCMC(L,noden,Visc,region,seistomo,scaling,step)
% [geoid,slepgeoid]=GEOIDKERNEL_MCMC(L,noden,Visc,region,seistomo,scaling,step)
%
% Returns the geoid, localized geoid, surface deflection, and localized
% surface deflection for a specific REGION. The SEISTOMO and SCALE determine
% the seismic model input and the relevant scaling factor as a function of 
% depth. VISC determines the radial viscosity profile. L gives the maximum
% bandwidth and NODEN determines the minimum depth for which the seismic 
% model is used. This method uses geoid kernels combined with Slepian 
% functions to return the localized geoid and dynamic topography for a REGION.
%
% INPUT:
%
% L        Bandwidth (maximum angular degree) [default: 20]
% noden    Minumum depth [km] for which seismic tomography model is used.
%          Above this depth, geoid and topography kernels are not computed.
%          OR Vector of depths whose contributions to the geoid and dynamic
%          topography are to be zeroed out. [default: 300]
% Visc     A two-columned matrix in which the first column contains
%          normalized interface radii and the second column contains
%          viscosities normalized by the reference viscosity. OR 
%          'PREM' [default] OR
%          'isoviscous'
% region   Angular extent of a spherical cap, in degrees OR
%          'england', 'eurasia',  'namerica', 'australia', 'greenland', 
%          'africa', 'samerica', 'amazon', 'orinoco', 'antarctica', 
%          'contshelves', 'alloceans' OR
%          [lon lat] an ordered list defining a closed curve [degrees] OR
%          {'region' buf} where buf is the distance in degrees that 
%          the region outline will be enlarged by BUFFERM
% seistomo String name of seismic tomography model you wish to use for the 
%          density structure of the mantle. Must be one of {'S40RTS',
%          'SEMUCB_WM1','S362WMANIM','SEISGLOB2','SAW642ANb','STB00'}. 
%          [default: S40RTS]
% scaling  Either a single scalar for uniform shear velocity to density
%          conversion OR a two column matrix with an arbitrary number of 
%          radii in the first column and the corresponding scaling 
%          factors in the second column. Can also input 's20rtsgiven' to 
%          use the density conversion from Simmons et al. (2007). 
%          [default: 0.3]
% step     Factor by which to step through depths in the mantle [km]. Must 
%          be a multiple of 10. If using model 'STB00', step is ignored.
%          [default: 10]
%
% OUTPUT:
%
% geoid     cosi mastrix for angular degree 0 to L with coefficients of the
%           full geoid.
% slepgeoid cosi mastrix for angular degree 0 to L with coefficients of the
%           geoid localized over REGION.
%
% SEE ALSO:
%
% DYN_KERNELS_CH, RETURNPROPAGATOR, VS2DENSITY, LOADTOMOSH, ADDMON, GLMALPHA
%
% Last modified by kengourley-at-arizona.edu, 8-13-2024
warning('off','all')

% defval('L',20);
% defval('noden',300);
% defval('toplot',true);
% defval('region','none');
% defval('seistomo','S40RTS');
% defval('scaling',0.3);
% defval('step',10);

if isequal(scaling,'s20rtsgiven')
    %Load scaling values as a function of depth 
    scaling=load('Simmonsetal2007scaling.mat','myline');
    scaling=scaling.myline;
end

%Get the dynamic kernel for every km of depth in the mantle
R=6371; %radius of Earth's surface [km]
CMB=3476; %radius of CMB [km]
Grav=6.6743*10^(-11); %gravitational constant

%Compute kernels for non-zero depths only
alllayers=0:step:R-CMB-5;
if isscalar(noden)
    lowerlayers=alllayers>(noden-1);
    depths=alllayers(lowerlayers);
elseif isvector(noden)
    nonzerolayers=~ismember(alllayers,noden);
    depths=alllayers(nonzerolayers);
else 
    error('noden must be a scalar or a vector of depths.')
end

%Normalize depths
Bp=(R*ones(size(depths))-depths)/R;
%Normalize viscosities
% defval('Visc',[(R-670)/R, 10;...
%                (R-400)/R, 1;...
%                (R-60)/R, 1/30;...
%                R/R, 100]);

%Compute geoid kernels for every depth for every sph har 0 to L
[~,K,~,~]=dyn_kernels_MCMC(0:L,Bp,Visc,0,'Earth',1);

%Zero out degree 0 and degree 1 contribution
K = [zeros(size(K(:,1:2))) K(:,3:end)];

%Compute sph har for bandwidth L
degord=((L+1)^2+L+1)/2;
[~,dels,~,lmcosi,~,~,~,bigl,~,ronm]=addmon(L);

%Compute the density structure at all depths
cosirho=vs2density_MCMC(depths,seistomo,L,scaling).*step;

%Compute the product of kernels and density structure
Ko=K(:,dels+1); Ko(:,:,2)=K(:,dels+1);
cosi3=cosirho.*Ko;

%Integrate and normalize the global geoid
R=6371000;
cosinorm=Grav*R*cosi3./repmat(2*dels'+1,length(depths),1,2);
geoid=squeeze(sum(cosinorm,1));

slepgeoid=zeros(degord,2);

if ~isequal(region,'none')
    %Determine eigens for region
    [G,V,~,~,N,~,~,~]=glmalpha(region,L);
    [~,k]=sort(V,'descend');
    G=G(:,k);

    %Repeat the above for the localized geoid
    R=6371000;
    Knorm=Grav*R*K.*repmat((1./(2*[0:L]+1)),size(K,1),1);
    GN=G(:,1:round(N));

    %Preallocate matrices used in for-loop
    lmcosirho=zeros(size(cosirho,2),size(cosirho,3));
    KGi=zeros(size(G,1),size(G,2));
    RHOalpha=zeros(round(N),1);
    slepgeoidU=zeros((L+1)^2,1);
    cosi=zeros(size(lmcosi,2),2);

    %Matrix multiplication must be performed inside loop
    parfor i=1:length(depths)
        %Pick out density and kernels at depth i
        lmcosirho=squeeze(cosirho(i,:,:));
        KGi = G.*repmat(Knorm(i,bigl+1)',1,(L+1)^2);

        %Localize product of density and geoid kernels
        RHOalpha=KGi(:,1:round(N))'*lmcosirho(ronm);
        slepgeoidU=GN*RHOalpha;
        cosi=lmcosi(:,3:4);
        cosi(ronm)=slepgeoidU(:);
        slepgeoid=slepgeoid+cosi;

    end
end

% if toplot
%     %Whole geoid from kernel method
%     figure
%     plotplm([lmcosi(:,1:2) geoid],[],[],4,1)
%     title('Global Geoid')
%     h1 = colorbar;
%     set(get(h1,'title'),'string','[meters]');
%     %Whole dynamic topography
%     figure
%     plotplm([lmcosi(:,1:2) dyntopo],[],[],4,1)
%     title('Global Dynamic Topography')
%     h2 = colorbar;
%     set(get(h2,'title'),'string','[meters]');
%     %Local geoid from kernel method
%     figure
%     plotplm([lmcosi(:,1:2) slepgeoid],[],[],4,1)
%     title('Localized Geoid')
%     h3 = colorbar;
%     set(get(h3,'title'),'string','[meters]');
%     %Local dynamic topography
%     figure
%     plotplm([lmcosi(:,1:2) sleptopo],[],[],4,1)
%     title('Localized Dynamic Topography')
%     h4 = colorbar;
%     set(get(h4,'title'),'string','[meters]');
% end

varns={geoid(:),slepgeoid(:)};
varargout=varns(1:nargout);
