%--------------------------------------------------------------------------
%% Problem : Linear-Elastic Plate bending
%--------------------------------------------------------------------------
% Problem : Find the maximum load a RC slab can carry, and its eventual
% yield-line patterns
% Two Boundary conditions are used, simply supported and clamped.
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: Dec 2018
%--------------------------------------------------------------------------

clc
clear all
close all
% profile on
% LinearMain
% profile viewer
% parpool % Enable parallel looping
tic
disp('Please wait Programme is running')

%%  Input data  %%
% Geometrical and material properties of plate
a = 1;                           % Length of the plate (along X-axes) [m]
b = 1;                           % Length of the plate (along Y-axes) [m]
E = 30e6;                         % Elastic Modulus [GPa --> kN/m^2]
nu = 0.3;                         % Poisson's ratio
t = 0.1 ;                         % Plate thickness [m]
I = t^3/12 ;                      % Moment of Inertia [m^4]

%% Number of mesh element in x and y direction %%
nelem_x = 20;
nelem_y = 20;
disp(['Selected mesh is: ' num2str(nelem_x)...
    ' by ' num2str(nelem_y)]) ;

if mod(nelem_x,2) == 1 || mod(nelem_y,2) == 1
    error('Amount of elements along each axis has to be even')
end

%% Nodal Connectivity %%
[ElemNode] = TopoQ4(nelem_y,nelem_x) ; 
% for node connectivity counting starts from 1 towards +y axis and new row 
% again start at y = 0
nodes = createNodes(a, b, nelem_x, nelem_y); % Node number + Node Coordinates
NodeCoor = nodes(:,2:3);                     % Node Coordinates

nelem = length(ElemNode) ;            % number of elements
nnelem = 4;                           % number of nodes per element
ndof = 3;                             % number of dofs per node
nnode = length(NodeCoor) ;            % total number of nodes in system
sysdof = nnode*ndof;                  % total system dofs  
elemdof = nnelem*ndof;                % degrees of freedom per element

% x- and y-coordinates in each element
ElemX = zeros(nelem,nnelem) ;     % Row no. i in ElemX contains x-nodal coordinates for element number i
ElemY = zeros(nelem,nnelem) ;     % Row no. i in ElemY contains y-nodal coordinates for element number i
        
for nodeno = 1:nnelem 
	ElemX(:,nodeno) = NodeCoor(ElemNode(:,nodeno),1) ;
	ElemY(:,nodeno) = NodeCoor(ElemNode(:,nodeno),2) ;
end
                             
% Display generated mesh
ElemDraw2D(ElemX,ElemY,'-',0,1:nelem,nodes,2)

%% Input Cross-section and Material data %%
fcu = 30;        % Concrete compressive strength [MPa] 
fyd = 450;       % Steel yield strength [MPa]  
B = b/nelem_x;   % Slab width [m]
H = 180;         % Slab height [mm]
C = 30;          % Distance of rebars from bottom/top of the slab [mm]

reb_lay_elem = zeros(nelem_x*nelem_y,1);

%% Define amount of additive rebar layers %%
reb_lay_am = 1; 
disp(['Amount of additive rebar layers: ' num2str(reb_lay_am)]);

% Layers overwrite one on the top of the other. Take that into account

%% Define master layer [0] %%
% Longitudinal Bars along x-direction
phisx0 = 8;           % Rebars diameter on the lower edge [mm]
phix0 = 8;            % Rebars diameter on the upper edge [mm]
sp_nsx0 = 0.15;        % Spacing of rebars on the lower edge [mm]
sp_nx0 = 0.15;         % Spacing of rebars on the upper edge [mm]
nsx0 = 1/sp_nsx0;  % Amount of rebars on the lower edge per m [-]
nx0 = 1/sp_nx0;    % Amount of rebars on the upper edge per m [-]
% Longitudinal Bars along y-direction
phisy0 = phisx0;      % Rebars diameter on the lower edge [mm]
phiy0 = phix0;        % Rebars diameter on the upper edge [mm]
nsy0 = nsx0;          % Amount of rebars on the lower edge per m [-]
ny0 = nx0;            % Amount of rebars on the upper edge per m [-]
sp_nsy0 = 1/nsy0;     % Spacing of rebars on the lower edge [mm]
sp_ny0 = 1/ny0;       % Spacing of rebars on the upper edge [mm]

% Plot Rebar Layer area
rl = figure;
set(rl,'name','Slab Design','numbertitle','off');
layer0 = fill([0,a,a,0],[0,0,b,b],'cyan','facealpha',.2);
legend(layer0,'Slab');
daspect([1 1 1])
set(gca,'fontsize', 20);
axis on ;
xlabel('x (m)'); ylabel('y (m)');
xlim([0 a]);
ylim([0 b]);
hold on

if reb_lay_am ~= 0
    if reb_lay_am >= 1
        
%% Define Rebar Layer [1] %%
% Longitudinal Bars along x-direction
phisx1 = 8;     % Rebars diameter on the lower edge [mm]
phix1 = 8;      % Rebars diameter on the upper edge [mm]
nsx1 = 6;        % Amount of rebars on the lower edge per m [-]
nx1 = 6;         % Amount of rebars on the upper edge per m [-]
sp_nsx1 = 1/nsx1; % Spacing of rebars on the lower edge [mm]
sp_nx1 = 1/nx1;   % Spacing of rebars on the upper edge [mm]
% Longitudinal Bars along y-direction
phisy1 = phisx1;     % Rebars diameter on the lower edge [mm]
phiy1 = phix1;      % Rebars diameter on the upper edge [mm]
nsy1 = nsx1;        % Amount of rebars on the lower edge per m [-]
ny1 = nx1;         % Amount of rebars on the upper edge per m [-]
sp_nsy1 = 1/nsy1; % Spacing of rebars on the lower edge [mm]
sp_ny1 = 1/ny1;   % Spacing of rebars on the upper edge [mm]

% Define Rebar Layer Area according to inserted coordinates of nodes
% (COUNTERCLOCKWISE!)
reb1nodes = [0.2 0.2; 0.4 0.2; 0.4 0.4; 0.2 0.4];
[reb_lay_elem] = lay_area(reb1nodes,... % Coordinates of rebar nodes [m]
    ElemNode,reb_lay_elem,NodeCoor,nelem_x,nelem_y,1,'blue'); 
legend('Master Rebar Layer 0','Rebar Layer 1');

%% Define Rebar Layer [2] %%
if reb_lay_am >= 2
% Longitudinal Bars along x-direction
phisx2 = 10;     % Rebars diameter on the lower edge [mm]
phix2 = 10;      % Rebars diameter on the upper edge [mm]
nsx2 = 4;        % Amount of rebars on the lower edge per m [-]
nx2 = 2;         % Amount of rebars on the upper edge per m [-]
sp_nsx2 = 1/nsx2; % Spacing of rebars on the lower edge [mm]
sp_nx2 = 1/nx2;   % Spacing of rebars on the upper edge [mm]
% Longitudinal Bars along y-direction
phisy2 = 10;     % Rebars diameter on the lower edge [mm]
phiy2 = 10;      % Rebars diameter on the upper edge [mm]
nsy2 = 4;        % Amount of rebars on the lower edge per m [-]
ny2 = 2;         % Amount of rebars on the upper edge per m [-]
sp_nsy2 = 1/nsy2; % Spacing of rebars on the lower edge [mm]
sp_ny2 = 1/ny2;   % Spacing of rebars on the upper edge [mm]

% Define Rebar Layer Area according to inserted coordinates of nodes
% (COUNTERCLOCKWISE!)
reb2nodes = [2 2; 4 2; 4 4; 2 4];
[reb_lay_elem] = lay_area(reb2nodes,... % Coordinates of rebar nodes [m]
    ElemNode,reb_lay_elem,NodeCoor,nelem_x,nelem_y,1,'magenta');
legend('Master Rebar Layer 0','Rebar Layer 1','Rebar Layer 2');
end
    end
end

%% Compute Plastic Moments %%
% Master Layer [0]
disp('Computing Maximum Resisting Moments for master rebar layer 0:')
% [Mppx0,Mnpx0] = PlasticMoments(fcu,B,H,C,phisx0,nsx0,phix0,nx0); 
% % Sagging Moments positive - Hogging Moments Negative
% disp('along x-direction:')
% disp(['Bottom Rebars: phi' num2str(phisx0) ' @ ' num2str(sp_nsx0*1000) ' mm']);
% disp(['Top Rebars: phi' num2str(phix0) ' @ ' num2str(sp_nx0*1000) ' mm']);
% Plastic Moments are chosen arbitrarily
Mppx0 = 1; Mnpx0 = -1;
disp(['Maximum Resisting Sagging Moment: ' num2str(Mppx0) ' kNm']) ;
disp(['Maximum Resisting Hogging Moment: ' num2str(Mnpx0) ' kNm']) ;

% [Mppy0,Mnpy0] = PlasticMoments(fcu,B,H,C,phisy0,nsy0,phiy0,ny0); 
% disp('along y-direction:')
% disp(['Bottom Rebars: phi' num2str(phisy0) ' @ ' num2str(sp_nsy0*1000) ' mm']);
% disp(['Top Rebars: phi' num2str(phiy0) ' @ ' num2str(sp_ny0*1000) ' mm']);
% Plastic Moments are chosen arbitrarily
Mppy0 = 1; Mnpy0 = -1;
disp(['Maximum Resisting Sagging Moment: ' num2str(Mppy0) ' kNm']) ;
disp(['Maximum Resisting Hogging Moment: ' num2str(Mnpy0) ' kNm']) ;

if reb_lay_am >= 1
    % Rebar Layer 1
    [Mppx1,Mnpx1] = PlasticMoments(fcu,fyd,B,H,C,phisx1,nsx1,phix1,nx1); 
    % Sagging Moments positive - Hogging Moments Negative
    disp('Computing Maximum Resisting Moments for rebar layer 1:')
    disp('along x-direction:')
    disp(['Bottom Rebars: phi' num2str(phisx1) ' @ ' num2str(sp_nsx1*1000) ' mm']);
    disp(['Top Rebars: phi' num2str(phix1) ' @ ' num2str(sp_nx1*1000) ' mm']);
    disp(['Maximum Resisting Sagging Moment: ' num2str(Mppx1) ' kNm']) ;
    disp(['Maximum Resisting Hogging Moment: ' num2str(Mnpx1) ' kNm']) ;

    [Mppy1,Mnpy1] = PlasticMoments(fcu,fyd,B,H,C,phisy1,nsy1,phiy1,ny1); 
    disp('along y-direction:')
    disp(['Bottom Rebars: phi' num2str(phisy1) ' @ ' num2str(sp_nsy1*1000) ' mm']);
    disp(['Top Rebars: phi' num2str(phiy1) ' @ ' num2str(sp_ny1*1000) ' mm']);
    disp(['Maximum Resisting Sagging Moment: ' num2str(Mppy1) ' kNm']) ;
    disp(['Maximum Resisting Hogging Moment: ' num2str(Mnpy1) ' kNm']) ;
end

if reb_lay_am >= 2
    % Rebar Layer 2
    [Mppx2,Mnpx2] = PlasticMoments(fcu,fyd,B,H,C,phisx2,nsx2,phix2,nx2); 
    % Sagging Moments positive - Hogging Moments Negative
    disp('Computing Maximum Resisting Moments for rebar layer 2:')
    disp('along x-direction:')
    disp(['Bottom Rebars: phi' num2str(phisx2) ' @ ' num2str(sp_nsx2*1000) ' mm']);
    disp(['Top Rebars: phi' num2str(phix2) ' @ ' num2str(sp_nx2*1000) ' mm']);
    disp(['Maximum Resisting Sagging Moment: ' num2str(Mppx2) ' kNm']) ;
    disp(['Maximum Resisting Hogging Moment: ' num2str(Mnpx2) ' kNm']) ;

    [Mppy2,Mnpy2] = PlasticMoments(fcu,fyd,B,H,C,phisy2,nsy2,phiy2,ny2); 
    disp('along y-direction:')
    disp(['Bottom Rebars: phi' num2str(phisy2) ' @ ' num2str(sp_nsy2*1000) ' mm']);
    disp(['Top Rebars: phi' num2str(phiy2) ' @ ' num2str(sp_ny2*1000) ' mm']);
    disp(['Maximum Resisting Sagging Moment: ' num2str(Mppy2) ' kNm']) ;
    disp(['Maximum Resisting Hogging Moment: ' num2str(Mnpy2) ' kNm']) ;
end

%% Choose Order of Nielsen Criterion %%
 YieldCriterionOrder = 1;         % Linearized Yield Criterion
% YieldCriterionOrder = 2;          % Quadratic Yield Criterion

if YieldCriterionOrder < 1 || YieldCriterionOrder > 2 
    error('Admissible yield criterion orders are 1 and 2')
end

%% Order of Gauss Quadrature %%
GaussQuadratureBendOrder = 'second'; % Select 'first' or 'second'
GaussQuadratureShearOrder = 'first'; %   //
if strcmp(GaussQuadratureBendOrder,'first')
    glb = 1; % 1x1 Gauss-Legendre quadrature for bending 
else 
    glb = 4; % 2x2 Gauss-Legendre quadrature for bending 
end
if strcmp(GaussQuadratureShearOrder,'first')
    gls = 1; % 1x1 Gauss-Legendre quadrature for shear 
else
    glb = 4; % 2x2 Gauss-Legendre quadrature for shear
end

%% Define type of BCs %%
%bc = 0; % Only external BCs
% bc = 1; % External + 1 internal BCs
 bc = 2; % External + 2 internal BCs

r2 = figure;
set(r2,'name','Slab BCs','numbertitle','off');
layer0(1) = fill([0,a,a,0],[0,0,b,b],'cyan','facealpha',.2);
legend(layer0(1),'Slab');
daspect([1 1 1]);
set(gca,'fontsize', 20);
axis on ;
xlabel('x (m)'); ylabel('y (m)');
xlim([0 a]);
ylim([0 b]);
hold on

%% External Boundary conditions (at the edges) %%
    %typeextBC = 's-s-s-s' ;        % Simply Supported
    %typeextBC = 'c-c-c-c'   ;      % Fully Clamped
    %typeextBC = 's-s-o-s'   ;      % Simply Supported at 3 out of 4 edges
    %typeextBC = 'o-s-o-s'   ;      % Simply Supported at 2 out of 4 edges
    %typeextBC = 's-o-s-o'   ;      % Simply Supported at 2 out of 4 edges
    %typeextBC = 's-o-o-s'   ;      % Simply Supported at 2 out of 4 edges 
    typeextBC = 'o-s-o-o'   ;       % Simply Supported at right edge    

% Establish which dofs to apply the BCs
extbcdof = BoundaryCondition(typeextBC,NodeCoor) ;

if bc == 1
%% Internal Boundary conditions 1 %%
%   typeintBC = 's-s-s-s' ;        % Simply Supported
%   typeintBC = 'c-c-c-c'   ;      % Fully Clamped
%   typeintBC = 'points';          % Point support
  typeintBC = 'patchs, s-s-s-s'; % Patch support
%   typeintBC = 'patchs, c-c-c-c'; % Patch clamp

    
    %dir = 'x'; % BC along x-direction
    %dir = 'y'; % BC along y-direction
    %dir = 'point'; % BC in one point support
    dir = 'patch'; % BC in a patch with simple support at edges
    
    coor_x1 = 0; % Define x-coordinate of internal BC
    coor_y1 = 0.967; % Define y-coordinate of internal BC
    coor_x2 = 0.1; % Define x-coordinate of internal BC
    coor_y2 = 1; % Define y-coordinate of internal BC
    
intbcdof = IntBoundaryCondition(typeintBC,NodeCoor,dir,...
    coor_x1,coor_y1,coor_x2,coor_y2,nelem_x,'b') ;
legend('Slab','Patch BC1');
end

if bc == 2
%% Internal Boundary conditions 2 %%

% 1st Internal Boundary condition
%   typeintBC1 = 's-s-s-s' ;        % Simply Supported
%   typeintBC1 = 'c-c-c-c'   ;      % Fully Clamped
%   typeintBC1 = 'points';          % Point support
   typeintBC1 = 'patchs, s-s-s-s';  % Patch support
%   typeintBC1 = 'patchs, c-c-c-c'; % Patch clamp

    
    %dir1 = 'x'; % BC along x-direction
    %dir1 = 'y'; % BC along y-direction
    %dir1 = 'point'; % BC in one point support
    dir1 = 'patch'; % BC in a patch with simple support at edges
    
    % Bottom left vertex
    % ROUNDED TO THREE SIGNIFICANT DIGITS!
    coor1_x1 = 0; % Define x-coordinate of internal BC
    coor1_y1 = 0; % Define y-coordinate of internal BC
    % Top right vertex
    coor1_x2 = 0.1; % Define x-coordinate of internal BC
    coor1_y2 = 0.1; % Define y-coordinate of internal BC
    
intbcdof1 = IntBoundaryCondition(typeintBC1,NodeCoor,dir1,...
    coor1_x1,coor1_y1,coor1_x2,coor1_y2,nelem_x,'b') ;

% 2nd Internal Boundary condition
%   typeintBC2 = 's-s-s-s' ;        % Simply Supported
%   typeintBC2 = 'c-c-c-c'   ;      % Fully Clamped
%   typeintBC2 = 'points';          % Point support
    typeintBC2 = 'patchs, s-s-s-s'; % Patch support
%    typeintBC2 = 'patchs, c-c-c-c'; % Patch clamp

    
    %dir2 = 'x'; % BC along x-direction
    %dir2 = 'y'; % BC along y-direction
    %dir2 = 'point'; % BC in one point support
    dir2 = 'patch'; % BC in a patch with simple support at edges
    
    % Bottom left vertex
    % ROUNDED TO THREE SIGNIFICANT DIGITS!
    coor2_x1 = 0; % Define x-coordinate of internal BC
    coor2_y1 = 0.9; % Define y-coordinate of internal BC
    % Top right vertex
    coor2_x2 = 0.1; % Define x-coordinate of internal BC
    coor2_y2 = 1; % Define y-coordinate of internal BC
    
intbcdof2 = IntBoundaryCondition(typeintBC2,NodeCoor,dir2,...
    coor2_x1,coor2_y1,coor2_x2,coor2_y2,nelem_x,'k') ;
legend('Slab','Internal BC1','Internal BC2');
% Merge the two Internal Boundary Conditions together
intbcdof = [intbcdof1, intbcdof2];
intbcdof = sort(intbcdof);
end

if bc == 0
    bcdof = extbcdof;
    for i = 1:length(bcdof)
    bcdof(any(bcdof == 0)) = []; % Remove all 0 entries
    end
elseif bc >= 1
    bcdof = [extbcdof,intbcdof];
    bcdof = sort(bcdof);
    for i = 1:length(bcdof)
    bcdof(any(bcdof == 0)) = []; % Remove all 0 entries
    end
else
    error('Please choose appropriate bcdof type')
end

%% Pressure on plate %%
MaxLoadSteps = 2;                                                      % Define Maximum Amount of Load Steps    
UnitP = -1 ;                                                           % Define Unit Load [kN/m^2]
Pmax = 1300;                                                             % Define Maximum Load [kN/m^2]
Pstart = 0;                                                            % Starting Load [kN/m^2]
lambda = Pstart:abs((Pmax - Pstart)/(MaxLoadSteps - 1)):abs(Pmax);     % Define Incrementing Load Factor [-]
lambda(1) = [];                                                        % Discard Lambda = 0

%% Define load type
%loadType = 0; % Uniform transverse pressure
loadType = 1; % One patch load applied
%loadType = 2; % Two patch loads applied
%loadType = 3; % Three patch loads applied
%loadType = 4; % Four patch loads applied

r3 = figure;
set(r3,'name','Slab Loads','numbertitle','off');
layer0(1) = fill([0,a,a,0],[0,0,b,b],'cyan','facealpha',.2);
legend(layer0(1),'Slab');
daspect([1 1 1]);
set(gca,'fontsize', 20);
axis on ;
xlabel('x (m)'); ylabel('y (m)');
xlim([0 a]);
ylim([0 b]);
hold on

load_elem = zeros(nelem,1);
if loadType == 0
    disp('Load Condition: Uniform Transverse Load')
elseif loadType == 1 % one patch load applied on specified area
    disp('Load Condition: One patch load')
    % The concept of lay_area is now recycled for determining modelling
    % patch loads
    % INSERT COORDINATES ROUNDED WITH 3 SIGNIFICANT FIGURES
    load1nodes = [0.4 0.4; 0.5 0.4; 0.5 0.5; 0.4 0.5];
    [load_elem] = lay_area(load1nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,1,'black');
legend('Slab','Patch Load');
elseif loadType == 2 % one patch load applied on specified area
    disp('Load Condition: Two patch loads')
    % Patch load 1
    load1nodes = [0.66 0.249; 0.69 0.249; 0.69 0.311; 0.66 0.311];
    [load_elem1] = lay_area(load1nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,1,'black');
    % Patch load 2
    load2nodes = [0.81 0.249; 0.84 0.249; 0.84 0.311; 0.81 0.311];
    [load_elem2] = lay_area(load2nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,2,'red');
legend('Slab','Patch Load 1','Patch Load 2');
load_elem = load_elem1 + load_elem2;
elseif loadType == 3 % one patch load applied on specified area
    disp('Load Condition: Three patch loads')
    % Patch load 1
    load1nodes = [0.66 0.249; 0.69 0.249; 0.69 0.311; 0.66 0.311];
    [load_elem1] = lay_area(load1nodes,... % Coordinates of fourth node [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,1,'black');
    % Patch load 2
    load2nodes = [0.81 0.249; 0.84 0.249; 0.84 0.311; 0.81 0.311];
    [load_elem2] = lay_area(load2nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,2,'red');
    % Patch load 3
    load3nodes = [0.81 0.249; 0.84 0.249; 0.84 0.311; 0.81 0.311];
    [load_elem3] = lay_area(load3nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,3,'green');
legend('Slab','Patch Load 1','Patch Load 2','Patch Load 3');
load_elem = load_elem1 + load_elem2 + load_elem3;
elseif loadType == 4 % one patch load applied on specified area
    disp('Load Condition: Four patch loads')
    % Patch load 1
    load1nodes = [4 3; 4.5 3; 4.5 3.5; 4 3.5];    
    [load_elem1] = lay_area(load1nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,1,'black');
    % Patch load 2
    load2nodes = [7.5 3; 8 3; 8 3.5; 7.5 3.5];
    [load_elem2] = lay_area(load2nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,2,'red');
    % Patch load 3
    load3nodes = [7.5 8.5; 8 8.5; 8 9; 7.5 9];
    [load_elem3] = lay_area(load3nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,3,'green');
    % Patch load 4
    load4nodes = [4 8.5; 4.5 8.5; 4.5 9; 4 9];
    [load_elem4] = lay_area(load4nodes,... % Coordinates of load nodes [m]
        ElemNode,load_elem,NodeCoor,nelem_x,nelem_y,4,'blue');
legend('Slab','Patch Load 1','Patch Load 2','Patch Load 3','Patch Load 4');
load_elem = load_elem1 + load_elem2 + load_elem3 + load_elem4;
else
    error('Four load types are supported so far');
end

%% ----------------- START ANALYSIS------------------------------------%%

for LoadStep = 1:MaxLoadSteps-1                    % Loop over each load increment starting here                       
P = lambda(LoadStep)*UnitP;                      % Load Intensity [kN/m^2]
disp(['Load Step ' num2str(LoadStep) ': lambda = ' num2str(lambda(LoadStep))]);
disp('Elastic FE Analysis Running');

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
% integration one point Guass Quadrature is used.

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

%% Compute Bending Strains & Bending Moments
% Extract displacements ordered by node
w = displacement(1:3:sysdof) ;      % Vertical Displacement
thetax = displacement(2:3:sysdof) ; % Rotation about x-axis
thetay = displacement(3:3:sysdof) ; % Rotation about y-axis

%% Maximum transverse displacement %%
format long 
minw = min(w) ;
% store(LoadStep+1,1) = P;
% store(LoadStep+1,2) = minw;
% P = P-1;  % This is increment of load for next step.

% % Table with values of w, thetax, thetay
% mytable(nnode,displacement,sysdof) ;

%% Plots %%
% Deformed Shape
x = NodeCoor(:,1) ;
y = NodeCoor(:,2) ;
f3 = figure ;
set(f3,'name','Postprocessing','numbertitle','off') ;
plot3(x,y,w,'.') ;
title('plate deformation') ;

%% Contour Plots %%
PlotFieldonDefoMesh(NodeCoor,ElemNode,w,w)
title('Profile of w on deformed Mesh') ;    

PlotFieldonMesh(NodeCoor,ElemNode,w)
title('Profile of w on plate')

% Create operative matrices
ElemW = zeros(nelem,nnelem);
ElemThetax = zeros(nelem,nnelem);
ElemThetay = zeros(nelem,nnelem);
ElemDispl = zeros(elemdof,nelem);

for ielem = 1:nelem
    for innelem = 1:nnelem
    % Create Topology Matrices for Displacements
    ElemW(ielem, innelem) = w(ElemNode(ielem, innelem));
    ElemThetax(ielem, innelem) = thetax(ElemNode(ielem, innelem));
    ElemThetay(ielem, innelem) = thetay(ElemNode(ielem, innelem));
    
    % Create a Matrix with indices for assembling the above-displacements
    IndexMatrix = zeros(ndof,nnelem);
         for indof = 1:ndof
         posrow = 1:ndof:elemdof;
         IndexMatrix(indof,:) = (indof-1) + posrow;
         end
         
    % Assemble displacements in one matrix
    ElemDispl(IndexMatrix(1,innelem)',ielem) = ElemW(ielem,innelem);
    ElemDispl(IndexMatrix(2,innelem)',ielem) = ElemThetax(ielem,innelem);
    ElemDispl(IndexMatrix(3,innelem)',ielem) = ElemThetay(ielem,innelem);   
    end
end

    % Calculate Strains and Moments as function of z.
    % Each column corresponds to an element (1:nelem), and rows
    % corresponds to a value of, respectively, Curvature or Moment,
    % 1:3 as x,y,xy for node 1
    % 4:6 as x,y,xy for node 2 and so on

    %------------------------Create function BendMoments-----------------
StrainsElemNode = zeros(size(ElemDispl,1),size(ElemDispl,2));
Moments = StrainsElemNode;
StrainsNodes = zeros(ndof,nelem);
% The strain interpolation matrix in each single node, is applied on the
% displacement vector containing values of displacements in all nodes in
% column
for i = 1:nelem
    for j = 1:nnelem
        StrainsElemNode(posrow(j):posrow(j)+ndof-1,i) = ...
            B_pbnode(posrow(j):posrow(j)+ndof-1,poscol(i):...
            poscol(i)+ndof*nnelem-1)*ElemDispl(:,i);
    end
end

% The strains within each element at the same node do not give back the
% same values, hence we average the strains within each element to produce
% reliable and symmetric results
Strains = zeros(ndof,nelem);
posarray = [posrow; posrow+1; posrow+2];
for i = 1:ndof
    for j = 1:nelem
        Strains(i,j) = mean(StrainsElemNode(posarray(i,:),j));
    end
end

% Obtain Moments
Moments = -D_pb*Strains*10^(-4);

%% Compute Principal Bending Moments %%
% 1) Solve Eigenvalue Problem - Output Principal Moments and Directions
PrincipalM = zeros(2,nelem); 
PrincipalD = zeros(4,nelem);
for ielem = 1:nelem
BendTensor = [Moments(1,ielem) Moments(3,ielem); ...
    Moments(3,ielem) Moments(2,ielem)];
[D,Mp] = eig(BendTensor); 
PrincipalM(:,ielem) = diag(Mp);
PrincipalD(1:length(D),ielem) = D(1,:);
PrincipalD(1+length(D):length(D)+length(D),ielem) = D(2,:);
end

%% Plot Moment Field on Mesh
% Invert Moments
momx = Moments(1,:);
Moments(1,:) = Moments(2,:);
Moments(2,:) =  momx;
clear momx

 PlotMomentFieldonMesh(NodeCoor,ElemNode,Moments,lambda,...
    LoadStep,max(Mppx0,Mppy0),min(Mnpx0,Mnpy0))          
  PlotPrincipalMomentFieldonMesh(NodeCoor,ElemNode,PrincipalM,lambda,...
     LoadStep,max(Mppx0,Mppy0),min(Mnpx0,Mnpy0)) 

% Mx, My, Mxy
% Extract values of Bending Moments in Vectors
Mx = Moments(1,:);
My = Moments(2,:);
Mxy = Moments(3,:);              

%% Compute Resisting Mxy (referring to the linearized criterion) %%
nmomel_x = 11;            % Choose amount of elements in moment field

[NegSol0, PosSol0, MomNodeCoor0] = ComputeMxyNielsenCriterion1...
    (nmomel_x,Mppx0,Mnpx0,Mppy0,Mnpy0);
if reb_lay_am ~= 0
if reb_lay_am >= 1
    [NegSol1, PosSol1, MomNodeCoor1] = ComputeMxyNielsenCriterion1...
    (nmomel_x,Mppx1,Mnpx1,Mppy1,Mnpy1);
end
if reb_lay_am >= 2
    [NegSol2, PosSol2, MomNodeCoor2] = ComputeMxyNielsenCriterion1...
    (nmomel_x,Mppx2,Mnpx2,Mppy2,Mnpy2);
end
end

%% Check numerically if the Yield Criterion is exceeded in any location %%
% Check whether the newly Reduced Distributed Moment field is yielding
% somewhere

% Check yield is run such as the whole slab would be constituted by
% different rebar layers
if reb_lay_am == 0
   [YieldIndex0,OvYieldMom0,ElastMom0,distBis0] = CheckYield(...
      Mppx0,Mnpx0,Mppy0,Mnpy0,Mx,My,Mxy); % Yield indices for rebar layer 0
elseif reb_lay_am == 1
   [YieldIndex0,OvYieldMom0,ElastMom0,distBis0] = CheckYield(...
      Mppx0,Mnpx0,Mppy0,Mnpy0,Mx,My,Mxy); % Yield indices for rebar layer 0
   [YieldIndex1,OvYieldMom1,ElastMom1,distBis1] = CheckYield(...
       Mppx1,Mnpx1,Mppy1,Mnpy1,Mx,My,Mxy); % Yield indices for rebar layer 1
elseif reb_lay_am == 2
   [YieldIndex0,OvYieldMom0,ElastMom0,distBis0] = CheckYield(...
      Mppx0,Mnpx0,Mppy0,Mnpy0,Mx,My,Mxy); % Yield indices for rebar layer 0
   [YieldIndex1,OvYieldMom1,ElastMom1,distBis1] = CheckYield(...
       Mppx1,Mnpx1,Mppy1,Mnpy1,Mx,My,Mxy); % Yield indices for rebar layer 1
   [YieldIndex2,OvYieldMom2,ElastMom2,distBis2] = CheckYield(...
       Mppx2,Mnpx2,Mppy2,Mnpy2,Mx,My,Mxy); % Yield indices for rebar layer 2
end

YieldIndex = zeros(1,nelem);
OvYieldMom = zeros(ndof,nelem);
ElastMom = zeros(ndof,nelem);
distBis = zeros(nelem,1);

% Pick the right values according to rebar layer areas
for i = 1:nelem
    if reb_lay_elem(i) == 0
        YieldIndex(i) = YieldIndex0(i);
        OvYieldMom(:,i) = OvYieldMom0(:,i);
        ElastMom(:,i) = ElastMom0(:,i);
        distBis(i) = distBis0(i);
    elseif reb_lay_elem(i) == 1
        YieldIndex(i) = YieldIndex1(i);
        OvYieldMom(:,i) = OvYieldMom1(:,i);
        ElastMom(:,i) = ElastMom1(:,i);
        distBis(i) = distBis1(i);
    elseif reb_lay_elem(i) == 2
        YieldIndex(i) = YieldIndex2(i);
        OvYieldMom(:,i) = OvYieldMom2(:,i);
        ElastMom(:,i) = ElastMom2(:,i);
        distBis(i) = distBis2(i);
    end
end

% Clear unuseful variables
clear YieldIndex0 OvYieldMom0 ElastMom0 distBis0
if reb_lay_am >= 1
    clear YieldIndex1 OvYieldMom1 ElastMom1 distBis1
end
if reb_lay_am >= 2
    clear YieldIndex2 OvYieldMom2 ElastMom2 distBis2
end

%% Plot Yield Criterion with (Elastic) Over-Yielded Solution %%
if YieldCriterionOrder == 1 
    if reb_lay_am == 2
    NielsenCriterion1_Elastic(nmomel_x,NegSol0,PosSol0,MomNodeCoor0,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,0) 
    NielsenCriterion1_Elastic(nmomel_x,NegSol1,PosSol1,MomNodeCoor1,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,1) 
    NielsenCriterion1_Elastic(nmomel_x,NegSol2,PosSol2,MomNodeCoor2,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,2)
    elseif reb_lay_am == 1
    NielsenCriterion1_Elastic(nmomel_x,NegSol0,PosSol0,MomNodeCoor0,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,0) 
    NielsenCriterion1_Elastic(nmomel_x,NegSol1,PosSol1,MomNodeCoor1,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,1) 
    else
    NielsenCriterion1_Elastic(nmomel_x,NegSol0,PosSol0,MomNodeCoor0,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,0)
    end
elseif YieldCriterionOrder == 2
    if reb_lay_am == 2
    NielsenCriterion2_Elastic(nmomel_x,Mppx0,Mnpx0,Mppy0,Mnpy0,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,0);
    NielsenCriterion2_Elastic(nmomel_x,Mppx1,Mnpx1,Mppy1,Mnpy1,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,1);
    NielsenCriterion2_Elastic(nmomel_x,Mppx2,Mnpx2,Mppy2,Mnpy2,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,2);
    elseif reb_lay_am == 1
    NielsenCriterion2_Elastic(nmomel_x,Mppx0,Mnpx0,Mppy0,Mnpy0,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,0);
    NielsenCriterion2_Elastic(nmomel_x,Mppx1,Mnpx1,Mppy1,Mnpy1,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,1);
    else 
    NielsenCriterion2_Elastic(nmomel_x,Mppx0,Mnpx0,Mppy0,Mnpy0,OvYieldMom,...
        ElastMom,lambda,LoadStep,length(Mx),YieldIndex,reb_lay_elem,0);
    end
end
end
toc
