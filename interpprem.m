function varargout=interpprem(depth)
% [Pdepth,Vp,Vs,Den]=interpprem(depth)
%
%
% INPUT:
%
% depth     The depth you want.
%
% OUTPUT:
% 
% Depth     Your depth back to you. 
% Vp        P wave speed (m/s)
% Vs        S wave speed (m/s)
% Den       Density (kg/m^3)
%
% Last modified by charig @ princeton.edu, 05/28/2013

% Load the PREM model
mystruct = load(fullfile(getenv('IFILES'),'EARTHMODELS','MATFILES','premiso'));
radius = mystruct.radius;
psd1=mystruct.psd;

% Make a combined matrix with depth instead of radius
myprem = flipud([(-radius + radius(end)) psd1]);

% Find what level our depth is in.
indeks = find(myprem(:,1) >= depth,1);
% This code ensures that if we get to a discontinuity, then we use the value
% above that discontinuity as the value there.

% And return that value
if depth == myprem(indeks,1)
    % Use the upper value of the discontinuity
    Pdepth = myprem(indeks,1);
    Vp = myprem(indeks,2);
    Vs = myprem(indeks,3);
    Den = myprem(indeks,4);
else
    % Make an interpolation
    Pdepth = interp1(myprem([indeks-1 indeks],1),myprem([indeks-1 indeks],1),depth);
    Vp = interp1(myprem([indeks-1 indeks],1),myprem([indeks-1 indeks],2),depth);
    Vs = interp1(myprem([indeks-1 indeks],1),myprem([indeks-1 indeks],3),depth);
    Den = interp1(myprem([indeks-1 indeks],1),myprem([indeks-1 indeks],4),depth);
end




%save(fnpl,'K','r','th','N','V','MTAP')
  
% Prepare output
varns={Pdepth,Vp,Vs,Den};
varargout=varns(1:nargout);