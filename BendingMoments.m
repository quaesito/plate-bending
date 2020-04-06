function [BendStrains,BendMoments]=BendingMoments...
    (displacement,ElemNode,B_pb,D_pb,ndof,sysdof)

%----------------------------------------------------------
%  Purpose:
%     Evaluate Bending Strains and Moments in each element
%
%  Synopsis:
%     [BendStrains,BendMoments]=BendingMoments(displacement,ElemNode,B_pb,D_pb,ndof,sysdof)
%
%  Variable Description:
%     BendStrains - Bending Strains 
%     BendMoments - Bending Moments
%     B_pb  - Matrix for Kinematic Equation of Bending
%     D_pb  - matrix for material property for bending
%     ndof - number of d.o.f. per node
%     sysdof - number of d.o.f. in the whole system
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: June 2017
%--------------------------------------------------------------------------

nelem = length(ElemNode) ;          % Number of elements
nnelem = size(ElemNode,2);          % Number of nodes per element
elemdof = ndof*nnelem;              % Number of dof per element

% Extract displacements ordered by node
w = displacement(1:3:sysdof) ;      % Vertical Displacement
thetax = displacement(2:3:sysdof) ; % Rotation about x-axis
thetay = displacement(3:3:sysdof) ; % Rotation about y-axis

% Create operative matrices
ElemW = zeros(nelem,nnelem);
ElemThetax = zeros(nelem,nnelem);
ElemThetay = zeros(nelem,nnelem);
ElemDispl = zeros(elemdof,nelem);
BendStrains = zeros(ndof,nelem);
BendMoments = zeros(ndof,nelem);

for ielem = 1:nelem
    for innelem = 1:nnelem
    % Create Topology Matrices for Displacements
    ElemW(ielem, innelem) = w(ElemNode(ielem, innelem));
    ElemThetax(ielem, innelem) = thetax(ElemNode(ielem, innelem));
    ElemThetay(ielem, innelem) = thetay(ElemNode(ielem, innelem));
    
    % Create a Matrix with indices for assembling the above-displacements
    IndexMatrix = zeros(ndof,nnelem);
         for indof = 1:ndof
         ipos = 1:ndof:elemdof;
         IndexMatrix(indof,:) = (indof-1) + ipos;
         end
         
    % Assemble displacements in one matrix
    ElemDispl(IndexMatrix(1,innelem)',ielem) = ElemW(ielem,innelem);
    ElemDispl(IndexMatrix(2,innelem)',ielem) = ElemThetax(ielem,innelem);
    ElemDispl(IndexMatrix(3,innelem)',ielem) = ElemThetay(ielem,innelem);
    
    % Calculate BendStrains and Bending Moments as function of z.
    % Each column corresponds to an element (1:nelem), and each row
    % corresponds to a value of Bending Moment or Curvature, respectively,
    % 1:3 as x,y,xy.
    BendStrains = B_pb*ElemDispl;
    
    % Extract Bending Moments per each node
    BendMoments = D_pb*BendStrains;
    end  
end