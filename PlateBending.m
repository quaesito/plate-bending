function [pb]=PlateBending(nnel,dshapedx,dshapedy)

%--------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic matrix for bending 
%
%  Synopsis:
%     [pb]=PlateBending(nnel,dshapedx,dshapedy) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dshapedx - derivatives of shape functions with respect to x   
%     dshapedy - derivatives of shape functions with respect to y
%--------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*3+1;  
 i2=i1+1;
 i3=i2+1;
 pb(1,i2)=-dshapedx(i);
 pb(2,i3)=-dshapedy(i);
 pb(3,i2)=-dshapedy(i);
 pb(3,i3)=-dshapedx(i);
 pb(3,i1)=0;
 end
