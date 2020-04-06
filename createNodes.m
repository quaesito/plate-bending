function nodes = createNodes(a, b, nelem_x, nelem_y)
% Defines nodes coordinates
%
% INPUT
%   a Length of the plate along x-axis [m]
%   b Length of the plate along y-axis [m]
%	nelem_y Number of elements in the vertical direction
%   nelem_x Number of elements in x-direction of the beam
%
% OUTPUT
%   ElemNode    Matrix containing the node numbers associated
%                 with its coordinate, e.g. 
%                 ElemNode = [... 
%                             2 4 0
%                             ...           ]
%                 The nodes are be numbered within increase in x-direction
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: June 2017
%--------------------------------------------------------------------------

stepx = a/nelem_x;  
stepy = b/nelem_y; 
nodes = [];
for y = 0 : stepy : b
    for x = 0 : stepx : a
        nodes = [nodes; size(nodes, 1) + 1 x y];
    end
end
end
