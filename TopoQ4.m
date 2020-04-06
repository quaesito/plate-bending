function [ElemNode] = TopoQ4(nelem_y,nelem_x)
% Defines the topology matrix of a plate built up by isoparametric
% four noded plane elements.
%
% INPUT
%	nelem_y Number of elements in the vertical direction
%   nelem_x Number of elements in x-direction of the beam
%
% OUTPUT
%   ElemNode    Matrix containing the node numbers associated
%                 with each element, e.g. 
%                 ElemNode = [... 
%                             8 12 3 6
%                             ...     ]
%                 The nodes are be numbered counterclockwise. Dim = nelem x 4,
%				  where nelem is the number of elements.
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: June 2017
%--------------------------------------------------------------------------

% Calculation of topology matrix ElemNode 

NodeTopo = zeros(nelem_x+1,nelem_y+1) ; 

for colnr = 1:nelem_y+1
    NodeTopo(:,colnr) = [(colnr-1)*(nelem_x+1)+1 : colnr*(nelem_x+1)]' ;
end

nelem = nelem_x*nelem_y ; % Total number of elements

ElemNode = zeros(nelem,4) ;

elemnr = 1 ;
for colnr = 1:nelem_y
    for rownr = 1:nelem_x		
        % Counter-clockwise numbering of nodes
		ElemNode(elemnr,2) = NodeTopo(rownr+1,colnr)   ; % Lower right node
		ElemNode(elemnr,3) = NodeTopo(rownr+1,colnr+1) ; % Upper right node
	    ElemNode(elemnr,4) = NodeTopo(rownr,colnr+1)   ; % Upper left node
		ElemNode(elemnr,1) = NodeTopo(rownr,colnr)     ; % Lower left node
		
		elemnr = elemnr + 1 ;
    end
end
%---------------------------------------------------------