function [PosSol, NegSol, NodeCoor] = ComputeMxyNielsenCriterion1(nmomel_x,Mppx,Mnpx,Mppy,Mnpy)
%--------------------------------------------------------------------------
%  Purpose:
%     Derive Linearized Nielsen's Yield Criterion 
%     for RC slabs
%
%  Synopsis:
%     ...
%
%  Variable Description:
%     ...
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: June 2017
%--------------------------------------------------------------------------

% Define Amount of elements along x- and y-axis (only even numbers are admissible)
nmomel_y = nmomel_x;                % Amount of elements on y-axis

if mod(nmomel_x,2) == 0
    error('Please choose even number of elements');
end

% Generate nodes in function of the desired amount of elements
nodes = createNodesPlasticMoments(Mppx, Mnpx, Mppy, Mnpy, nmomel_x, nmomel_y); % Node number + Node Coordinates
NodeCoor = nodes(:,2:3);           % Node Coordinates

% Create vectors with increasing bending moments in x- and y-directions
MxVec = NodeCoor(1:nmomel_x + 1,1);
MyVec = NodeCoor(1:nmomel_y + 1:end,2);

%% Linearized Yield Criterion

% Indices Vectors to be used in the loop
ix1 = 1:(nmomel_x + 1);             % Ascending order of indices
iy =  -1*sort(-(1:(nmomel_y + 1))); % Descending order of indices
% Create Operative Matrices
NegSol = zeros(length(ix1),1);
PosSol = NegSol;

% Loop deriving the roots of Mxy :
% Only the points along the bisectrice from (Mnpx,Mppy) to (Mppx,Mnpy) are
% needed
for k = 1:length(ix1)               % Ascending order of indices
syms MxyValues 
eqnLeftCone = (abs(Mnpx) + MxVec(ix1(k))).*(abs(Mnpy) + MyVec(iy(k))) == MxyValues.^2;
SolLeftCone = solve(eqnLeftCone, MxyValues);
NegSol(k) = double(SolLeftCone(1));
PosSol(k) = double(SolLeftCone(2));
end

end

