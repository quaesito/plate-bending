function [YieldIndex,OvYieldMom,ElastMom,distBis] = ...
    CheckYield(Mppx,Mnpx,Mppy,Mnpy,Mx,My,Mxy)
%----------------------------------------------------------
%  Purpose:
%     Check numerically if the Yield Criterion is exceeded in any location
%
%  Synopsis:
%     [YieldIndex,OvYieldMom,ElastMom] = CheckYield(Mppx,Mnpx,Mppy,Mnpy,Mx,My,Mxy)
%
%  Variable Description:
%     YieldIndex - 0 Not Yielding / 1 Yielding
%     OvYieldMom - Over Yielded Moments ~0 / Non-Yielded Moments = 0
%     ElastMom  - Over Yielded Moments = 0 / Non-Yielded Moments ~= 0
%     MomP - Amount on moment points in the physical domain
%     distBis - Distances of each moment point to the bisectrice
%     Mppx,Mnpx,Mppy,Mnpy  - Resisting Moments
%     Mx,My,Mxy - Acting Moments
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: Nov 2017
%--------------------------------------------------------------------------

% Define Operative Matrices 
% Distances of each moment point to the bisectrice (negative if on the side
% of Cone2)
distBis = zeros(length(Mx),1);
% Cone1 is the one on the left side of the bisectrice, Cone2 is on the
% right one
Cone1 = distBis; Cone2 = distBis; 
% Yield will determine whether yielding is occurring or not
OvYieldMom = zeros(3,length(Mx)); 
ElastMom = OvYieldMom;
YieldIndex = zeros(1,length(Mx));

% Start Loop over every point of loading Mx,My,Mxy
for MomP = 1:1:length(Mx)
% Calculate distances between point and the bisectrice in the Mx-My-plane
% to determine the location of the point
distBis(MomP) = ((Mnpy - Mppy)*Mx(MomP) - (Mppx - Mnpx)*My(MomP) ...
    + Mppx*Mppy - Mnpy*Mnpx)/sqrt((Mnpy - Mppy).^2 + (Mppx - Mnpx).^2);

% Calculate the Mxy of the surface for each couple of points (Mx,My)
[Mnpxy1, Mp, Mnpxy2, Mppxy2] = ComputeMrxyatMxMylocations(Mnpx,Mppx,Mnpy,Mppy,Mx,My);

% Locate each point (Mx,My) respectively in the left(1) of right(2) side of 
% the cone mid-line
   if distBis(MomP) > 0
     if Mx(MomP) < Mnpx  || Mx(MomP) > Mppx || My(MomP) < Mnpx || My(MomP) > Mppy
     YieldIndex(MomP) = 1; 
     OvYieldMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';

     elseif Mxy(MomP) > 0
         if Mxy(MomP) > Mp(MomP)
         YieldIndex(MomP) = 1; 
         OvYieldMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
         else
         ElastMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
         end
         
      else % if Mxy(MomP) < 0
         if Mxy(MomP) < Mnpxy1(MomP)
         YieldIndex(MomP) = 1; 
         OvYieldMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
         else
         ElastMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
         end   
     end     
 
   elseif distBis(MomP) < 0
     if Mx(MomP) < Mnpx  || Mx(MomP) > Mppx || My(MomP) < Mnpx || My(MomP) > Mppy
     YieldIndex(MomP) = 1; % disp('Yielding is occurring')
     OvYieldMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
     
     elseif Mxy(MomP) > 0
         if Mxy(MomP) > Mppxy2(MomP)
         YieldIndex(MomP) = 1; % disp('Yielding is occurring')
         OvYieldMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
         else
         ElastMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';         
         end
      else % if Mxy(MomP) < 0
         if Mxy(MomP) < Mnpxy2(MomP)
         YieldIndex(MomP) = 1; % disp('Yielding is occurring')
         OvYieldMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';
         else
         ElastMom(:,MomP) = [Mx(MomP) My(MomP) Mxy(MomP)]';         
         end
     end     
   end
end
end