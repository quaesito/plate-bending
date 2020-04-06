function [Neighbours,NeighbYieldIndex,NeighbYield] = NeighboursYielding...
    (NodeCoor,ElemNode,nelem_x,nelem_y,lambda,increment,YieldIndex)

%--------------------------------------------------------------------------
% Purpose:
%         Determine a few neighborhood relationships, and useful operative
%         matrices
% For all of them row represents the respective number of the element.
% Neighbours = [left neighb, right neighb, bottom neighb, top neighb]
% NeighbYieldIndex = [0 if element not neighbour of a yielding element
%                     1 if element is a neighbour of a yielding element]
% NeighbYield = [left neighb, right neighb, bottom neighb, top neighb]
%               0 if neighb element is not yielding
%               1 if neighb element is yielding
% Synopsis :
%           [Neighbours,NeighbYieldIndex,NeighbYield] = NeighboursYielding(NodeCoor,ElemNode,nelem_x,nelem_y,lambda,increment,YieldIndex)
% Variable Description:
%           NodeCoor - The nodal coordinates of the mesh
%           -----> NodeCoor = [node X Y] 
%           ElemNode - The nodal connectivity of the elements
%           -----> ElemNode = [node1 node2......]    
%           nelem_x - amount of element in x-direction
%           nelem_y - amount of element in y-direction
%           lambda - Load multiplication factor
%           increment - number of increment
%           YieldIndex - Vector stating whether element is yielding or not
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: Oct 2017
%--------------------------------------------------------------------------

nelem = length(ElemNode) ;                % number of elements
nnode = length(NodeCoor) ;                % total number of nodes in system
nnelem = size(ElemNode,2);                % number of nodes per element

% Initialization of the required matrices
X = zeros(nnelem,nelem) ;
Y = zeros(nnelem,nelem) ;

%% Derive coordinates X and Y of each node
for iel=1:nelem   
     for i=1:nnelem
     nd(i)=ElemNode(iel,i);         % extract connected node for (iel)-th element
     X(i,iel)=NodeCoor(nd(i),1);    % extract x value of the node
     Y(i,iel)=NodeCoor(nd(i),2);    % extract y value of the node
     end
end
    
%% Create Neighbours Matrix %%
Neighbours = zeros(nelem,4);
% Create a neighbours matrix where :
% Column 1 and 2, are, respectively, number of elements of neighbours to
% the left and to the right, and Column 3 and 4 the ones down and up.
for i = 1:nelem
    Neighbours(i,1) = i - 1;
    Neighbours(i,2) = i + 1;
    Neighbours(i,3) = i - nelem_x;
    Neighbours(i,4) = i + nelem_x;
    for j = 1:size(Neighbours,2)
       if Neighbours(i,j) <= 0
          Neighbours(i,j) = 1;
       end
    end
end

   % Discard unappropriate results
   % Set them = 1 for allowing the next group of loops to run
for i = 1:nelem_y
    Neighbours(i*nelem_x,2) = NaN;       % Falling out of mesh to the right
    Neighbours(i*nelem_x + 1,1) = NaN;   % Falling out of mesh to the left
end
    Neighbours(1,1) = NaN;               % Discard left neighbour for element 1
       Neighbours(end,:) = [];              % Discard last row generated by loop
        Neighbours(end,2) = NaN;          % Discard right neighbour for last element

for i = 1:nelem_x
    Neighbours(i,3) = NaN;                              % Falling out of mesh to the bottom
    Neighbours(length(Neighbours) - (i - 1),4) = NaN;   % Falling out of mesh to the top
end

%% Determine Coordinates of Yielded Elements, and of their Neighbours
% Create Operative Matrices
YieldX = zeros(nnelem,nelem);
YieldY = YieldX;
NeighbYieldX = YieldX;
NeighbYieldY = YieldX;
NeighbYieldIndex = zeros(length(YieldIndex),1);

     for i = 1:length(YieldIndex)
         if YieldIndex(i) == 1
            % Determine Coordinates of Yielded Elements
            YieldX(:,i) = X(:,i);
            YieldY(:,i) = Y(:,i);
                 % Determine Neighbours of Yielded Elements
                 for j = 1:size(X,1)
                     if isnan(Neighbours(i,j))
                     continue
                     else
                     NeighbYieldIndex(Neighbours(i,j)) = 1;
                     end 
                 end
         end
     end
     
     % Discard Yielded Indices from Neighbours of Yield Indices
     NeighbYieldIndex = NeighbYieldIndex - YieldIndex';

     for i = 1:length(YieldIndex)
         if NeighbYieldIndex(i) == 1
             % Determine Coordinates of Neighbours of Yielded Elements
             NeighbYieldX(:,i) = X(:,i);
             NeighbYieldY(:,i) = Y(:,i);
         end
     end
     
%% Detect if any Neighbour Elements is Over-Yielded
NeighbYield = zeros(nelem,size(Neighbours,2));
% NeighbYield matrix shows which neighbour elements have yielded
% Order: left, right, bottom, and top neighbours
% Yielded = 1, Not Yielded = 0, No Neighbour = NaN
for i = 1:nelem
    for j = 1:size(Neighbours,2)  
        if isnan(Neighbours(i,j))
           NeighbYield(i,j) = 0;
           
        else
        NeighbYield(i,j) = YieldIndex(Neighbours(i,j));
        end
    end
end

    % Discard unappropriate results
    % No Neighborhood = NaN
for i = 1:nelem_y
    %Neighbours(i*nelem_y,2) = NaN;       % Falling out of mesh to the right
    NeighbYield(i*nelem_y,2) = NaN;
    %Neighbours(i*nelem_y + 1,1) = NaN;   % Falling out of mesh to the left
    NeighbYield(i*nelem_y + 1,1) = NaN;
end
    %Neighbours(1,1) = NaN;               % Discard left neighbour for element 1
    NeighbYield(1,1) = NaN;
    %Neighbours(end,2) = NaN;               % Discard right neighbour for last
    NeighbYield(end,2) = NaN;
    %Neighbours(end,:) = [];              % Discard last row generated by loop
    NeighbYield(end,:) = [];
    
for i = 1:nelem_x
    %Neighbours(i,3) = NaN;                              % Falling out of mesh to the bottom
    NeighbYield(i,3) = NaN;
    %Neighbours(length(Neighbours) - (i - 1),4) = NaN;   % Falling out of mesh to the top
    NeighbYield(length(Neighbours) - (i - 1),4) = NaN;
end

% Discard non-zero entries in element at elastic status, if they are
% yielded or non-yielded neighbour of yielded elements
for i = 1:size(X,2)
    if any(YieldX(:,i) ~= 0) ||  any(NeighbYieldX(:,i) ~= 0)
        X(:,i) = zeros(size(X,1),1);
        Y(:,i) = X(:,i);
    end
end

%% Plot FEM mesh, Yielded Elements (in Blue), Non-Yielded Neighbours
% of yielded elements (in Red), Elements at elastic status (in Green)
     f1 = figure ;
     set(f1,'name',['Yielded Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}']','numbertitle','off','Color','w') ;
     h1 = fill(X,Y,'g');
     hold on
     h2 = fill(YieldX,YieldY,'b');
     hold on
     h3 = fill(NeighbYieldX,NeighbYieldY,'r');
     leg = legend([h1(1); h2(2); h3(3)],'Non-yielded elements','Yielded elements','Non-yielded neighbours of yielded elements');
     %title('Finite Element Mesh') ;
     daspect([1 1 1])
     set(gca,'fontsize', 20);
     axis on ;
     xlabel('x (m)'); ylabel('y (m)');
     %title(['Yielded Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
