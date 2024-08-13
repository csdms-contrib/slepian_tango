function varargout=loadtomoSH_MCMC(model,depth)
% lmcosi=LOADTOMOSH(model,depth)
%
% Load the precomputed spherical harmonic coefficients for a specified
% MODEL and DEPTH. These files should be in IFILES/TOMOGRAPHY/{MODEL NAME}.
%
% INPUT:
%
% model      String name of model you wish to load. Must be one of
%            {'S40RTS','SEMUCB_WM1','S362WMANIM','SEISGLOB2','SAW642ANb',
%            'STB00'}. [default: S40RTS]
% depth      A scalar for the depth of MODEL you desire in km. Must be a 
%            multiple of 10 km. OR
%            A vector of depths you wish to load; must be a subset of
%            25:1:2890 for S40RTS
%            30:10:2890 for SAW642ANb;
%            80:10:2890 for SEISGLOB2;
%            90:10:2890 for SEMUCB_WM1 and S362WMANIM;
%            30:1:2890 for STB00. [default: 400]
%
% OUTPUT:
%
% lmcosi     The desired matrix of spherical harmonic coefficients.
%            If input DEPTH is a scalar, lmcosi has shape 861 x 4.
%            If input DEPTH is a vector, lmcosi has shape d x 861 x 4,
%            where d is the number of depths in the input vector.
%
% See also VS2DENSITY, GEOIDKERNEL
%
% Last modified by kengourley-at-arizona.edu 9/1/2018

% defval('depth',400);
% defval('model','S40RTS');

% Check if we're importing slab model depths
% if isequal(model,'STB00')
%     numb_depth=round(linspace(2890,60,45));
%     if ~ismember(depth, numb_depth)
%         error('Depth does not exist.');
%     end
% Else try to import one of the four tomography models
% else
if ~ismember(model,{'S40RTS','SEMUCB_WM1','S362WMANIM','SEISGLOB2','SAW642ANb','STB00'})
    error('Model does not exist. Please check your spelling.');
%     elseif (ischar(depth) && (~strcmp(depth,'all') || ~strcmp(depth,'all_1km'))) || (isscalar(depth) && ~ismember(depth,80:10:2890))
%         error('Depth does not exist.');
elseif (isvector(depth) && ~all(depth<2891)) || (isscalar(depth) && depth>2890)
    error('Requested depth(s) fall outside permitted range.')
elseif isequal(model,'S40RTS') && ((isvector(depth) && ~all(depth>24)) || (isscalar(depth) && depth<25))
    error('Requested depth(s) fall outside permitted range.')
elseif isequal(model,'SEISGLOB2') && ((isvector(depth) && ~all(depth>79)) || (isscalar(depth) && depth<80))
    error('Requested depth(s) fall outside permitted range.')
elseif (isequal(model,'SEMUCB_WM1') || isequal(model,'S362WMANIM')) && ((isvector(depth) && ~all(depth>89)) || (isscalar(depth) && depth<90))
    error('Requested depth(s) fall outside permitted range.')
elseif (isequal(model,'SAW642ANb') || isequal(model,'STB00')) && ((isvector(depth) && ~all(depth>29)) || (isscalar(depth) && depth<30))
    error('Requested depth(s) fall outside permitted range.')
end
%end

% If scalar depth, load single depth slice
if isscalar(depth)
    loadfile=sprintf('%s_%dkm.dat',model,depth);
    lmcosi=load(fullfile(getenv('IFILES'),'TOMOGRAPHY',model,loadfile));
% Else, load vector of depths
elseif isvector(depth)
    if isequal(model,'S40RTS')
        if min(depth) >= 80 && all(mod(depth,10) == 0)
            numb_depth=80:10:2890;
            loadfile=sprintf('%s_all.dat',model);
        else
            numb_depth=25:1:2890;
            loadfile=sprintf('%s_all_1km.dat',model);
        end
    elseif isequal(model,'SEISGLOB2')
        numb_depth=80:10:2890;
        loadfile=sprintf('%s_all.dat',model);
    elseif isequal(model,'SEMUCB_WM1') || isequal(model,'S362WMANIM')
        numb_depth=90:10:2890;
        loadfile=sprintf('%s_all.dat',model);
    elseif isequal(model,'SAW642ANb')
        numb_depth=30:10:2890;
        loadfile=sprintf('%s_all.dat',model);
    elseif isequal(model,'STB00')
        if min(depth) >= 30 && all(mod(depth,10) == 0)
            numb_depth=30:10:2890;
            loadfile=sprintf('%s_all.dat',model);
        else
            numb_depth=30:1:2890;
            loadfile=sprintf('%s_all_1km.dat',model);
        end
    end
    
    llmcosi=load(fullfile(getenv('IFILES'),'TOMOGRAPHY',model,loadfile));
    
    % Filter all depths to the ones requested
    
    ndepth=1:length(numb_depth);
    nd=ismember(numb_depth,depth);
    n=ndepth(nd);
    llmcosi=llmcosi(n',:,:);
    [~,llength]=size(llmcosi);
    % Reconstruct lmcosi matrix so it is depths x 861 x 4
    lmcosi=llmcosi(:,1:(llength/4));
    lmcosi(:,:,2)=llmcosi(:,(llength/4 + 1):(llength/2));
    lmcosi(:,:,3)=llmcosi(:,(llength/2 + 1):(3*llength/4));
    lmcosi(:,:,4)=llmcosi(:,(3*llength/4 + 1):llength);
else
    error('Depth must either be a single number or a vector of depths.');
end

% Provide output
varns={lmcosi};
varargout=varns(1:nargout);
