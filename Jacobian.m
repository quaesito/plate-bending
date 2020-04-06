function [detjacobian,invjacobian]=Jacobian(nnel,dshapedxi,dshapedeta,xcoord,ycoord)

%------------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian for two-dimensional mapping
%
%  Synopsis:
%     [detjacobian,invjacobian]=Jacobian(nnel,dshapedxi,dshapedeta,xcoord,ycoord) 
%
%  Variable Description:
%     jacobian - Jacobian 
%     nnel - number of nodes per element   
%     dshapedxi - derivative of shape functions w.r.t. natural coordinate xi
%     dshapedeta - derivative of shape functions w.r.t. natural coordinate eta
%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
%------------------------------------------------------------------------

 jacobian = zeros(2,2);

 for i=1:nnel
 jacobian(1,1) = jacobian(1,1) + dshapedxi(i)*xcoord(i);
 jacobian(1,2) = jacobian(1,2) + dshapedxi(i)*ycoord(i);
 jacobian(2,1) = jacobian(2,1) + dshapedeta(i)*xcoord(i);
 jacobian(2,2) = jacobian(2,2) + dshapedeta(i)*ycoord(i);
 end

 detjacobian = det(jacobian) ;  % Determinant of Jacobian matrix
 invjacobian = inv(jacobian) ;  % Inverse of Jacobian matrix