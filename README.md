# Bending of multilayer thin/thick plate under composite transverse pressures modeled using Finite Element Method

A multi-cross section thin plate under composite transverse pressure is considered with line/area boundary conditions. 

In main.m:
	
## 1) Input Data
'''bash
% Geometrical and material properties of plate
a = 1;                           % Length of the plate (along X-axes) [m]
b = 1;                           % Length of the plate (along Y-axes) [m]
E = 30e6;                         % Elastic Modulus [GPa --> kN/m^2]
nu = 0.3;                         % Poisson's ratio
t = 0.1 ;                         % Plate thickness [m]
I = t^3/12 ;                      % Moment of Inertia [m^4]
'''

## 2) Number of mesh element in x and y direction
'''bash
nelem_x = 30;
nelem_y = 30;
'''

## 3) Input Cross-section and Material data 
'''bash
fcu = 30;        % Concrete compressive strength [MPa] 
B = b/nelem_x;   % Slab width [m]
H = 180;         % Slab height [mm]
C = 30;          % Distance of rebars from bottom/top of the slab [mm]
'''

## 4) Define amount of additive rebar layers
% In sec 3, the master layer of the plate has been defined. Now, here specify whether you want any area of the plate to have a different cross-section than the master one. In such case, the new cross-section will overwrite the existing one/ones in the area where it was specified.

'''bash
reb_lay_am = 2; 
'''

## 5) Define master layer [0] 
'''bash
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
'''

## 6) Define Rebar Layer [1]
'''bash
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
reb1nodes = [1 1; 5 1; 5 5; 1 5];
'''

## 7) Define Rebar Layer [2]
'''bash
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
reb1nodes = [1 1; 5 1; 5 5; 1 5];
'''

## 8) Choose Order of Nielsen Criterion
'''bash
 YieldCriterionOrder = 1;         % Linearized Yield Criterion
% YieldCriterionOrder = 2;          % Quadratic Yield Criterion
'''

## 9) Define type of BCs 
'''bash
%bc = 0; % Only external BCs
% bc = 1; % External + 1 internal BCs
 bc = 2; % External + 2 internal BCs

%% External Boundary conditions (at the edges) %%
    %typeextBC = 's-s-s-s' ;        % Simply Supported
    %typeextBC = 'c-c-c-c'   ;      % Fully Clamped
    %typeextBC = 's-s-o-s'   ;      % Simply Supported at 3 out of 4 edges
    %typeextBC = 'o-s-o-s'   ;      % Simply Supported at 2 out of 4 edges
    %typeextBC = 's-o-s-o'   ;      % Simply Supported at 2 out of 4 edges
    %typeextBC = 's-o-o-s'   ;      % Simply Supported at 2 out of 4 edges 
    typeextBC = 'o-s-o-o'    ;      % Simply Supported at right edge   
'''



## 10) Internal Boundary conditions 1 
'''bash
if bc == 1
%   typeintBC = 's-s-s-s' ;        % Simply Supported
%   typeintBC = 'c-c-c-c'   ;      % Fully Clamped
%   typeintBC = 'points';          % Point support
  typeintBC = 'patchs, s-s-s-s';  % Patch support
%   typeintBC = 'patchs, c-c-c-c'; % Patch clamp
    
    %dir = 'x'; % BC along x-direction
    %dir = 'y'; % BC along y-direction
    %dir = 'point'; % BC in one point support
    dir = 'patch'; % BC in a patch with simple support at edges
    
    coor_x1 = 0; % Define x-coordinate of internal BC
    coor_y1 = 0.967; % Define y-coordinate of internal BC
    coor_x2 = 0.1; % Define x-coordinate of internal BC
    coor_y2 = 1; % Define y-coordinate of internal BC
'''


## 11) Internal Boundary conditions 2 %%
'''bash
if bc == 2
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
    coor1_x2 = 0.033; % Define x-coordinate of internal BC
    coor1_y2 = 0.033; % Define y-coordinate of internal BC
    
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
    coor2_y1 = 0.967; % Define y-coordinate of internal BC
    % Top right vertex
    coor2_x2 = 0.033; % Define x-coordinate of internal BC
    coor2_y2 = 1; % Define y-coordinate of internal BC
    
intbcdof2 = IntBoundaryCondition(typeintBC2,NodeCoor,dir2,...
    coor2_x1,coor2_y1,coor2_x2,coor2_y2,nelem_x,'k') ;
legend('Slab','Internal BC1','Internal BC2');
% Merge the two Internal Boundary Conditions together
intbcdof = [intbcdof1, intbcdof2];
intbcdof = sort(intbcdof);
'''

## 12) Define Load Steps and Intensity of Load
'''bash
MaxLoadSteps = 2;                                                      % Define Maximum Amount of Load Steps    
UnitP = -1 ;                                                           % Define Unit Load [kN/m^2]
Pmax = 1300;                                                           % Define Maximum Load [kN/m^2]
Pstart = 0;                                                            % Starting Load [kN/m^2]
lambda = Pstart:abs((Pmax - Pstart)/(MaxLoadSteps - 1)):abs(Pmax);     % Define Incrementing Load Factor [-]
lambda(1) = [];              
'''

## 13) Define load type
'''bash
%loadType = 0; % Uniform transverse pressure
loadType = 1; % One patch load applied
%loadType = 2; % Two patch loads applied
%loadType = 3; % Three patch loads applied
%loadType = 4; % Four patch loads applied
'''

Specify the load location of load 1 at load1nodes, of load 2 at load2nodes etc.

'''bash
%% ----------------- START ANALYSIS------------------------------------%%
'''
