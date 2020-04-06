function [nodeL1,nodeL2,nodeL3,nodeL4,ElemL] =...
    LNodeCoorElem(NodeCoor,nelem,nnelem,ElemNode)
%----------------------------------------------------------
%  Purpose:
%     Find indices of nodes and elements along boundaries
%
%  Synopsis:
%     [nodeL1,nodeL2,nodeL3,nodeL4,ElemL1,ElemL2,ElemL3,ElemL4] =...
%                             LNodeCoorElem(NodeCoor,nelem,nnelem,ElemNode)
%
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: May 2018
%--------------------------------------------------------------------------

%% Find Indices of NodeCoor along the boundaries %%
nodeL1 = find(NodeCoor(:,2) == min(NodeCoor(:,2))) ; % at y = 0 (along X-axes)
nodeL2 = find(NodeCoor(:,1) == max(NodeCoor(:,1))) ; % at x = a (along Y-axes)
nodeL3 = find(NodeCoor(:,2) == max(NodeCoor(:,2))) ; % at y = b (along X-axes)
nodeL4 = find(NodeCoor(:,1) == min(NodeCoor(:,1))) ; % at x = 0 (along Y-axes)

%% Find Indices of elements along boundaries %%% 
% Find Elements along L1
ElemL1 = zeros(nelem,1);
for inode = 1:length(nodeL1)
    for ielem = 1:nelem
        for innelem = 1:nnelem
    if nodeL1(inode) == ElemNode(ielem,innelem)
        ElemL1(ielem) = ielem;
    end
        end
    end
end
ElemL1(ElemL1 == 0) = [];

% Find Elements along L2
ElemL2 = zeros(nelem,1);
for inode = 1:length(nodeL2)
    for ielem = 1:nelem
        for innelem = 1:nnelem
    if nodeL2(inode) == ElemNode(ielem,innelem)
        ElemL2(ielem) = ielem;
    end
        end
    end
end
ElemL2(ElemL2 == 0) = [];

% Find Elements along L3
ElemL3 = zeros(nelem,1);
for inode = 1:length(nodeL3)
    for ielem = 1:nelem
        for innelem = 1:nnelem
    if nodeL3(inode) == ElemNode(ielem,innelem)
        ElemL3(ielem) = ielem;
    end
        end
    end
end
ElemL3(ElemL3 == 0) = [];

% Find Elements along L4
ElemL4 = zeros(nelem,1);
for inode = 1:length(nodeL4)
    for ielem = 1:nelem
        for innelem = 1:nnelem
    if nodeL4(inode) == ElemNode(ielem,innelem)
        ElemL4(ielem) = ielem;
    end
        end
    end
end
ElemL4(ElemL4 == 0) = [];

% if length(ElemL1) == length(ElemL2) &&...
%    length(ElemL1) == length(ElemL3) &&...
%    length(ElemL1) == length(ElemL4)
% ElemL = [ElemL1'
%          ElemL2'
%          ElemL3'
%          ElemL4'];

ElemL = zeros(4,max([length(ElemL1),length(ElemL2),length(ElemL3),length(ElemL4)]));
ElemL(1,1:length(ElemL1)) = ElemL1';
ElemL(2,1:length(ElemL2)) = ElemL2';
ElemL(3,1:length(ElemL3)) = ElemL3';
ElemL(4,1:length(ElemL4)) = ElemL4';

% else 
%     error('Please use same amount of elements along x and y directions')
end