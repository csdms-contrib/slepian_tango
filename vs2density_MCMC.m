function varargout=vs2density_MCMC(depth,model,L,scaling)
% cosirho=VS2DENSITY(depth,model,L,scaling)
%
% Compute the spherical harmonic coefficients of the absolute density 
% distribution of a layer at a given DEPTH using the MODEL of choice with a
% SCALING factor of choice up to a bandwidth L.
%
% INPUT:
%
% depth      Depth of MODEL you desire in km. Must be greater than or equal
%            to 80 km (except for SEMUCB_WM1, which has a minimum depth of 
%            90 km) and less than or equal to 2890 km. Must be a multiple 
%            of 10 km. OR
%            If using STB00, must be one of the pre-defined depths given in 
%            the directory (e.g. 60, 124, 189, etc.). [default: 400]
% model      String name of model you wish to load. Must be one of
%            {'S40RTS','SEMUCB_WM1','S362WMANIM','SEISGLOB2','SAW642ANb',
%            'STB00'}. [default: S40RTS]
% L          Bandwidth of spherical harmonic expansion. [default: 20]
% scaling    Either a single scalar for uniform shear velocity to density
%            conversion OR a two column matrix with an arbitrary number of 
%            radii in the first column and the corresponding scaling 
%            factors in the second column. [default: 0.3]
% 
% OUTPUT:
%
% cosirho  The desired matrix of spherical harmonic coefficients. Does not 
%          include degrees and orders (only two columns).
%
% See also LOADTOMOSH, INTERPPREM, GEOIDKERNEL
%
% Last modified by kengourley-at-arizona.edu 4/12/2018

% defval('depth',400);
% defval('model','S40RTS');
% defval('L',20);
% defval('scaling',0.3);

if L > 40 || (L > 31 && isequal(model,'STB00'))
    error('Bandwidth L must be less than 41, unless your model is STB00, in which case L must be less than 32.')
end

if isscalar(depth)
    %load model
    lmcosi=loadtomoSH_MCMC(model,depth);
    %filter to bandwidth L
    lmcosi=lmcosi(1:((L+1)^2+L+1)/2,:);
    %load PREM density
    [~,~,~,Den]=interpprem(depth*1000);

    if ismatrix(scaling) && isequal(size(scaling,2),2)
        scalev=interp1(scaling(:,1),scaling(:,2),depth,'linear');
    elseif isscalar(scaling)
        scalev=scaling;
    else
        error('Scaling factor must either be a single number or two column matrix of values.')
    end

    cosirho=lmcosi(:,3:4)*scalev*Den*10;
elseif isvector(depth)    
    %load model
    lmcosi=loadtomoSH_MCMC(model,depth);
    %filter to bandwidth L
    lmcosi=lmcosi(:,1:((L+1)^2+L+1)/2,3:4);
    [m,n,~]=size(lmcosi);
    Den=zeros(m,1);
    for i=1:length(depth)
        [~,~,~,rho]=interpprem(depth(i)*1000);
        Den(i)=rho;
    end
    
    if ismatrix(scaling) && isequal(size(scaling,2),2)
        scalev=zeros(m,1);
        for i=1:length(depth)
            scalev(i)=interp1(scaling(:,1),scaling(:,2),depth(i),'linear');
        end
    elseif isscalar(scaling)
        scalev=repmat(scaling,length(depth),1);
    else
        error('Scaling factor must either be a single number or two column matrix of values.')
    end
    cosirho=lmcosi.*repmat(scalev,1,n,2).*repmat(Den,1,n,2).*10;
end

% Provide output
varns={cosirho};
varargout=varns(1:nargout);