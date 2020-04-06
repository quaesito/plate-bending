function [lay_elem, firstnodex, firstnodey, secondnodex, secondnodey,...
                thirdnodex, thirdnodey, fourthnodex, fourthnodey] = ...
                lay_area(load_nodes,...
    ElemNode,lay_elem,NodeCoor,nelem_x,nelem_y,lay_id,color)
%-------------------------------------------------------------
% Define element at the perimeter of layer area and 
% mark elements who are part of such layer with a Lay_id
%
% INPUT  
%	firstnodex     x-coordinate of Node 1	
%	firstnodey     y-coordinate of Node 1
%   ...
%   fourthnodex    x-coordinate of Node 4
%   fourthnodey    y-coordinate of Node 4
%   ElemNode       Connectivity Matrix containing nodes of each element
%   lay_elem       Vector marking which layer each element belongs to 
%   NodeCoor       Coordinates of each node
%   nelem_x        Amount of elements in x-direction
%   nelem_y        Amount of elements in y-direction
%   lay_id         Layer ID
%   color          color for fill plots
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The University of Auckland
% Latest revision: October 2018
%--------------------------------------------------------------------------

% Extract node number according to inserted coordinates
lay_nodes = zeros(4,2);
for i = 1:size(NodeCoor,1)
    for j = 1:size(load_nodes,1)
        if load_nodes(j,:) == round(NodeCoor(i,:),3)
            lay_nodes(j,:) = i;
        end
    end
end

% Check if inserted rebar layer coordinates can be detected according to
% the inputted amount of elements
for i = 1:size(lay_nodes,1)
    if lay_nodes(i,:) == [0 0]
error('Selected amount of elements not appropriate for defined layer. Make sure the selected coordinates can be detected')
    end
end
% Check if nodes of layers are inserted counterclockwise
if lay_nodes(3) < lay_nodes(4)
    error('Nodes of layer areas must be inserted counterclockwise')
end
       
lay_elem_ind = zeros(1,4);
% Find four 'corner' elements of Rebar Layer area 1
for i = 1:size(ElemNode,1)   
    if ElemNode(i,1) == lay_nodes(1)      % Find Element having the specified
                                        % bottom left node
       lay_elem(i) = lay_id;    % Mark that element as being part of layer 1
       bl_elem_ind = i;                 % Store index of bottom left element
    end
    if ElemNode(i,2) == lay_nodes(2)      % Find Element having the specified
                                        % bottom right node
       lay_elem(i) = lay_id;    % Mark that element as being part of layer 1
       br_elem_ind = i;                 % Store index of bottom right element
    end
    if ElemNode(i,3) == lay_nodes(3)      % Find Element having the specified
                                        % top right node
       lay_elem(i) = lay_id;    % Mark that element as being part of layer 1
       tr_elem_ind = i;                 % Store index of top right element
    end
    if ElemNode(i,4) == lay_nodes(4)      % Find Element having the specified
                                        % top left node
       lay_elem(i) = lay_id;    % Mark that element as being part of layer 1       reb_lay_elem_ind(1) = i;     % Store index of element
       tl_elem_ind = i;                 % Store index of top left element
    end
end

%% Create a vector where indices of all elements in Layer 1 are stored %%
% 1) Case of distinct 'corner' elements
if bl_elem_ind ~= br_elem_ind && br_elem_ind ~= tr_elem_ind
lay_elem_ind = bl_elem_ind:br_elem_ind; % First row
for i = 1:nelem_y % All other rows
    lay_elem_ind(1+i,:) = [bl_elem_ind + i*nelem_x:...
                                         br_elem_ind + i*nelem_x];
    if any(any(lay_elem_ind == tr_elem_ind) == 1) && ...
            any(any(lay_elem_ind == tl_elem_ind) == 1) % until top left and right 
        break                                    % elements are found
    end
end
lay_elem_ind = sort(reshape(lay_elem_ind,1,[])); % Put the elements in the form of a vector


% 2) Case of one row of elements  
elseif bl_elem_ind ~= br_elem_ind && br_elem_ind == tr_elem_ind &&...
        bl_elem_ind == tl_elem_ind
          lay_elem_ind = bl_elem_ind:br_elem_ind; % First row
          
% 3) Case of one column of elements 
elseif bl_elem_ind == br_elem_ind && br_elem_ind ~= tr_elem_ind &&...
        tl_elem_ind == tr_elem_ind
          lay_elem_ind = bl_elem_ind:nelem_x:tl_elem_ind; % First column
          
% 4) Case of one single element 
elseif bl_elem_ind == br_elem_ind && br_elem_ind == tr_elem_ind &&...
        tl_elem_ind == tr_elem_ind
          lay_elem_ind = bl_elem_ind; % One element
end

% 5) Mark all elements in Layer
for i = 1:length(lay_elem_ind) 
    lay_elem(lay_elem_ind(i)) = lay_id;
end
 

%% Plot Layer Areas
% Plot Layer area
layerreb1 = fill([NodeCoor(ElemNode(bl_elem_ind,1),1),... % Layer 1
           NodeCoor(ElemNode(br_elem_ind,2),1),...
           NodeCoor(ElemNode(tr_elem_ind,3),1),...
           NodeCoor(ElemNode(tl_elem_ind,4),1)],...
           [NodeCoor(ElemNode(bl_elem_ind,1),2),...
           NodeCoor(ElemNode(br_elem_ind,2),2),...
           NodeCoor(ElemNode(tr_elem_ind,3),2),...
           NodeCoor(ElemNode(tl_elem_ind,4),2)],color,'facealpha',.2);
legend(layerreb1,'Rebar Layer')
        legend(gca,'off');
        legend('show')
hold on