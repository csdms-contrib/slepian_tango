function [PBR,PCR]=ReturnPropagator(B,C,l,MantleProperties)
%
% Returns the propogator matrix used to calculate geoid kernels.
%
% INPUT:
%
% B                  Normalized Depth. If you submit a vector of depths
%                     then PBR and PCR are returned for each depth with 
%                     size 4 x length(B)*4 
% C                  Normalized core radius (Core mantle boundary)
% l                  The spherical harmonic degree (single)
% MantleProperties   A 3xn matrix with non-dimensional radii in the first
%                    column, viscosities in the second column, and 
%                    densities (if desired) in the third column.  
%                    Radii in the first column should increase
%
% OUTPUT:
%
% PBR
% PCR
%
% SEE ALSO: DYN_KERNELS, GKERNEL
%
% Created by Peter James (pjames@mit.edu) 6/10/12
% Last modified by charig-at-princeton.edu on 06/01/2015
%

% Preallocate
PCR=repmat(eye(4),1,length(B)); 
PBR=repmat(eye(4),1,length(B));

% Do this for each depth
for j=1:length(B)

  % Check the input arguments
  if max(MantleProperties(:,1)) > 1
     warning('Interface radii should be normalized');
  elseif norm(MantleProperties(:,1)-sort(MantleProperties(:,1))) > 1E-8
     warning('Interface radii in MantleProperties are not correctly sorted');
  elseif MantleProperties(1,1) < C
     warning('Mantle Properties are specified deeper than the lower BC');
     disp(MantleProperties(1,1))
     disp(C)
  elseif MantleProperties(size(MantleProperties,1)) ~= 1
     warning('Properties should be specified all the way to the surface');
  elseif B(j) <= C
     warning('Mantle load should be above the lower BC')
  end


  if size(MantleProperties,2)>2
    error('I haven''t gotten to chemical stratification yet')         
  else
    for i=1:size(MantleProperties,1)
      rlayer=MantleProperties(i,1);
      if i==1
        rlayer2=C;
      else
        rlayer2=MantleProperties(i-1,1);
      end
      mustar=MantleProperties(i,2);
      A = [-2, l*(l+1), 0, 0;...
          -1, 1, 0, 1/mustar;...
          12*mustar, -6*l*(l+1)*mustar, 1, l*(l+1);...
          -6*mustar, 2*(2*l*(l+1)-1)*mustar, -1, -2];

      % Make PCR
      PCR(:,(4*j-3):(4*j))=expm((log(rlayer)-log(rlayer2))*A)*PCR(:,(4*j-3):(4*j));
      % Make PBR
      if B(j)<rlayer
        if B(j)>rlayer2
          PBR(:,(4*j-3):(4*j))=expm((log(rlayer)-log(B(j)))*A)*PBR(:,(4*j-3):(4*j));
        else
          PBR(:,(4*j-3):(4*j))=expm((log(rlayer)-log(rlayer2))*A)*PBR(:,(4*j-3):(4*j));
        end
      end
    end
  end
end