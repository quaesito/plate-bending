function [ps]=PlateShear(nnel,dshapedx,dshapedy,shape)

%------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic matrix for shear
%
%  Synopsis:
%     [ps]=PlateShear(nnel,dshapedx,dshapedy,shape) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dshapedx - derivatives of shape functions with respect to x   
%     dshapedy - derivatives of shape functions with respect to y
%     shape - shape function
%------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*3+1;  
 i2=i1+1;
 i3=i2+1;
 ps(1,i1) = dshapedx(i);
 ps(2,i1) = dshapedy(i);
 ps(1,i2) = -shape(i); 
 ps(2,i3) = -shape(i);
 end
