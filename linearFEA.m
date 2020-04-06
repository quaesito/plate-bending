function [displacement,posrow,poscol,B_pbnode,D_pb] = ...
    linearFEA(lambda,LoadStep,UnitP,sysdof,elemdof,...
    GaussQuadratureBendOrder,E,nu,ndof,nnelem,nelem,ElemNode,NodeCoor,...
    t,glb,gls,load_elem,loadType,bcdof)
%% RUN LINEAR FEA %%
P = lambda(LoadStep)*UnitP;                      % Load Intensity [kN/m^2]

%% Initialization of matrices and vectors %%
force = zeros(sysdof,1) ;             % system Force Vector
stiffness = zeros(sysdof,sysdof);     % system stiffness matrix
index = zeros(elemdof,1);             % index vector

%%  Computation of element matrices and vectors and their assembly %% 
%%  For bending stiffness
% Sampling points & weights
[pointb,weightb] = GaussQuadrature(GaussQuadratureBendOrder);     
D_pb = E/(1-nu^2)*[1  nu 0; nu  1  0; 0  0  (1-nu)/2];  % Constitutive Matrix

%%  For shear stiffness
[points,weights] = GaussQuadrature('first');    % sampling points & weights
G = 0.5*E/(1+nu);                               % shear modulus
shcof = 5/6;                                    % shear correction factor
D_ps = G*shcof*[1 0; 0 1];                      % material matrix for shear 

% Initialization matrix of strain interpolation matrix at each node
B_pbnode = zeros(ndof*nnelem,nelem*ndof*nnelem);
% Vectors useful for storing B_pb at each node
posrow = 1:ndof:size(B_pbnode,1);              % Position in row - multiplier of ndof
poscol = 1:ndof*nnelem:size(B_pbnode,2);       % Position in col - multiplier of ndof*nnelem

for ielem = 1:nelem                        % loop for the total number of elements
for innelem = 1:nnelem                     % loop for total number of nodes for single element   
node(innelem) = ElemNode(ielem,innelem);   % extract connected node for (iel)^th element
xx(innelem) = NodeCoor(node(innelem),1);   % extract x-value of the node
yy(innelem) = NodeCoor(node(innelem),2);   % extract y-value of the node
end

B_pb = zeros(3,elemdof);               % initialization of kinematic matrix for bending
B_ps = zeros(2,elemdof);               % initialization of kinematic matrix for shear
ke = zeros(elemdof,elemdof);           % initialization of element stiffness matrix 
kb = zeros(elemdof,elemdof);           % initialization of bending stiffness matrix 
ks = zeros(elemdof,elemdof);           % initialization of shear stiffness matrix 
felem = zeros(elemdof,1) ;             % initialization of force vector                   

%%  Numerical integration for bending term %%
for int = 1:glb                        % nglb is sampling points for bending
% Local coordinates xi and eta
xi = pointb(int,1);                    % sampling point in x-axis here value of xi is given
eta = pointb(int,2);                   % sampling point in y-axis. here value of eta is given
wt = weightb(int,1);                   % weights for sampling points

% Compute shape functions and derivatives w.r.t. local coordinates at sampling point
[shape,dshapedxi,dshapedeta] = Shapefunctions(xi,eta);    
% compute Jacobian
[detjacobian,invjacobian] = Jacobian(nnelem,dshapedxi,dshapedeta,xx,yy);  
% Convert shape function derivatives into physical coordinates x and y.
[dshapedx,dshapedy] = ShapefunctionDerivatives(nnelem,dshapedxi,...
    dshapedeta,invjacobian);
                                  % derivatives w.r.t. physical coordinates  
                                     
% The DOF for single element are ordered as u=[w thetax thetay], where 
% thetax is rotation about y axis
B_pb = PlateBending(nnelem,dshapedx,dshapedy);    % bending kinematic matrix

% Store strain interpolation matrix for each node in the looping element ielem 
B_pbnode(posrow(int):posrow(int)+(ndof-1),...
    poscol(ielem):poscol(ielem)+(ndof*nnelem-1)) = B_pb;

% Compute Bending element matrix
kb = kb + t^3/12*B_pb'*D_pb*B_pb*wt*detjacobian;
end   % end of numerical integration loop for bending term

%%  Numerical Integration for Shear term %%
for int = 1:gls
% Local coordinates xi and eta
xi = points(int,1);                % sampling point in x-axis
eta = points(int,2);               % sampling point in y-axis
wt = weights(int,1);               % weights for sampling points
% Compute shape functions and derivatives at sampling point
[shape,dshapedxi,dshapedeta] = Shapefunctions(xi,eta);   
% Compute Jacobian
[detjacobian,invjacobian] = Jacobian(nnelem,dshapedxi,dshapedeta,xx,yy); 
% Convert shape function derivatives into physical coordinates x and y.
[dshapedx,dshapedy] = ShapefunctionDerivatives(nnelem,dshapedxi,...
    dshapedeta,invjacobian); 
                                  % derivatives w.r.t. physical coordinates
% Compute Shear kinematic matrix 
B_ps = PlateShear(nnelem,dshapedx,dshapedy,shape);  

%  Compute shear element matrix                
ks = ks + t*B_ps'*D_ps*B_ps*wt*detjacobian;

% Converting uniform pressure into equivalent nodal forces. For this
% integration one point Gauss Quadrature is used.

%% Force vector in each element
fe = Force(nnelem,shape,P,ielem,load_elem,loadType) ;     % Force vector
felem = felem + fe*wt*detjacobian ;
end                      % end of numerical integration loop for shear term

%%  Compute Element Matrix %%
ke = kb + ks;
% Extract system dofs associated with element
index = elementdof(node,nnelem,ndof);  
% Assemble element Stiffness and Force matrices
[stiffness,force] = assemble(stiffness,force,ke,felem,index);  
end

% Modify stiffness and force accordingly.
[stiffness,force] = constraints(stiffness,force,bcdof);

%% Solution of Linear FE Problem %%
displacement = stiffness\force ;
