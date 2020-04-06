function bcdof = BoundaryCondition(typeBC,coordinates)

%--------------------------------------------------------------------------
%   Purpose:
%           To determine the boundary conditions degree of freedom
%   Synopsis:
%           bcdof = BoundaryCondition(typeBC,coordinates,loadstep)
%
%   Variable Description:
%           bcdof - boundary condition degree of freedom
%                   dofs are (w, thetax, thetay)
%           typeBC - string which gives the state of boundary condition
%           coordinates - geometric coordinates of nodes
%           
%--------------------------------------------------------------------------

L1 = find(coordinates(:,2) == min(coordinates(:,2))) ; % at y = 0 (along X-axes)
L2 = find(coordinates(:,1) == max(coordinates(:,1))) ; % at x = a (along Y-axes)
L3 = find(coordinates(:,2) == max(coordinates(:,2))) ; % at y = b (along X-axes)
L4 = find(coordinates(:,1) == min(coordinates(:,1))) ; % at x = 0 (along Y-axes)
m = length(L1) ;
n = length(L2) ;
 
if strcmp(typeBC,'s-s-s-s')
         disp('Boundary Conditions: The slab is simply supported at all the edges')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL1(i1) = 3*L1(i)-1;    % constraining thetax along x axis
            dofL1(i2) = 3*L1(i)-2 ;   % constraining w along x axis
            dofL3(i1) = 3*L3(i)-1 ;    
            dofL3(i2) = 3*L3(i)-2 ;   
        end   
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL2(i1) = 3*L2(i) ;    % constraining thetay along y axis
            dofL2(i2) = 3*L2(i)-2 ;   % constraining w along y axis
            dofL4(i1) = 3*L4(i) ;
            dofL4(i2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
end
if  strcmp(typeBC,'c-c-c-c')
        disp('Boundary Conditions: The slab is fully clamped')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL1(i1) = 3*L1(i);    % constraining thetay along x axis
            dofL1(i2) = 3*L1(i)-2 ; % constraining w 
            dofL3(i1) = 3*L3(i) ;
            dofL3(i2) = 3*L3(i)-2 ;
        end   
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL2(i1) = 3*L2(i)-1 ;  % constraining thetax along y axis
            dofL2(i2) = 3*L2(i)-2 ;  % constraining w 
            dofL4(i1) = 3*L4(i)-1 ;
            dofL4(i2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
end
if strcmp(typeBC,'s-s-o-s')
         disp('Boundary Conditions: The slab is simply supported at all edges except top one')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL1(i1) = 3*L1(i)-1;    % constraining thetax along x axis
            dofL1(i2) = 3*L1(i)-2 ;   % constraining w 
        end   
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL2(i1) = 3*L2(i) ;    % constraining thetay along y axis
            dofL2(i2) = 3*L2(i)-2 ;   % constraining w 
            dofL4(i1) = 3*L4(i) ;
            dofL4(i2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
end
if strcmp(typeBC,'o-s-o-s')
         disp('Boundary Conditions: The slab is simply supported at left and right edges')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ; 
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL2(i1) = 3*L2(i) ;    % constraining thetay along y axis
            dofL2(i2) = 3*L2(i)-2 ;   % constraining w along y axis
            dofL4(i1) = 3*L4(i) ;
            dofL4(i2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
end
if strcmp(typeBC,'s-o-s-o')
         disp('Boundary Conditions: The slab is simply supported at bottom and top edges')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL1(i1) = 3*L1(i)-1;    % constraining thetax along x axis
            dofL1(i2) = 3*L1(i)-2 ;   % constraining w along x axis
            dofL3(i1) = 3*L3(i)-1 ;    
            dofL3(i2) = 3*L3(i)-2 ;   
        end   
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
bcdof(any(bcdof == 0)) = []; % Remove any 0 entry
end

if strcmp(typeBC,'s-o-o-s')
         disp('Boundary Conditions: The slab is simply supported at bottom and left edges')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL1(i1) = 3*L1(i)-1;    % constraining thetax along x axis
            dofL1(i2) = 3*L1(i)-2 ;   % constraining w along x axis
        end   
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL4(i1) = 3*L4(i) ;
            dofL4(i2) = 3*L4(i)-2 ;
        end
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
bcdof(any(bcdof == 0)) = []; % Remove any 0 entry
end

if strcmp(typeBC,'o-s-o-o')
         disp('Boundary Conditions: The slab is simply supported at the right edge')
        dofL1 = zeros(1,2*m) ;
        dofL2 = zeros(1,2*n) ;
        dofL3 = zeros(1,2*m) ;
        dofL4 = zeros(1,2*n) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL2(i1) = 3*L2(i) ;    % constraining thetay along y axis
            dofL2(i2) = 3*L2(i)-2 ;   % constraining w along x axis
        end   
        L1UL3 = union(dofL1,dofL3) ;
        L2UL4 = union(dofL2,dofL4) ;
        bcdof = union(L1UL3,L2UL4) ;
bcdof(any(bcdof == 0)) = []; % Remove any 0 entry
end