function [DistrMom] = DistributeMxy(RedAmMom,YieldIndex,Neighbours,DistrMom)
%----------------------------------------------------------
%  Purpose:
%     Distribute the amount of Mxy that was reduced when dragging back the
%     moments on the surface among all non-yielded nieghbours
%
%  Synopsis:
%     [DistrMom] = DistributeMxy(RedAmMom,YieldIndex,Neighbours,DistrMom)
%
%  Variable Description:
%     RedAmMom - [Mx; My; Mxy] Amount of Reduced Moment when dragging back 
%                Over-Yielded Moments of the surface
%     YieldIndex - 0 Not Yielding / 1 Yielding
%     Neighbours - [left neighb, right neighb, bottom neighb, top neighb]
%     DistrMom - [Mx; My; Mxy] of Distributed Moment (only distributed Mx
%                is included as input. Mx and My will be included as output
%                
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: Nov 2017
%--------------------------------------------------------------------------
% Find Indices (Elements number) of Elements whose Mxy Reduction was positive and negative
PosRedAmMxyIndex = find(RedAmMom(3,:) > 0); % - because if RedAmMom > 0, then the
NegRedAmMxyIndex = find(RedAmMom(3,:) < 0);

% Find Their Non-Yielded Neighbours
PosRedAmMxyNeighb = zeros(length(1));
PosRedAmMxy = zeros(length(PosRedAmMxyIndex),1);

for i = 1:length(PosRedAmMxyIndex)
    for j = 1:size(Neighbours,2)
        if isnan(Neighbours(PosRedAmMxyIndex(i),j))
        continue
        end
        % Check that the neighbour is not yielded
        if YieldIndex(Neighbours(PosRedAmMxyIndex(i),j)) == 0
            % Extract it in a vector
            PosRedAmMxyNeighb(end + 1) = Neighbours(PosRedAmMxyIndex(i),j);
            % Calculate Total Positive Reduced Amount of Mxy
            PosRedAmMxy(i) = RedAmMom(3,PosRedAmMxyIndex(i));
        end
    end
end
PosRedAmMxyNeighb(1) = [];

% Delete entries of Neighbours that are counted twice
PosRedAmMxyNeighb = sort(PosRedAmMxyNeighb);
for i = 1 : length(PosRedAmMxyNeighb) - 1
    if PosRedAmMxyNeighb(i) == PosRedAmMxyNeighb(i + 1);
        PosRedAmMxyNeighb(i) = 0;
    end
end
PosRedAmMxyNeighb = nonzeros(PosRedAmMxyNeighb);

NegRedAmMxyNeighb = zeros(length(1));
for i = 1:length(NegRedAmMxyIndex)
    for j = 1:size(Neighbours,2)
        if isnan(Neighbours(NegRedAmMxyIndex(i),j))
        continue
        end
        % Check that the neighbour is not yielded
        if YieldIndex(Neighbours(NegRedAmMxyIndex(i),j)) == 0
            % Extract it in a vector
            NegRedAmMxyNeighb(end + 1) = Neighbours(NegRedAmMxyIndex(i),j);
            % Calculate Total Negative Reduced Amount of Mxy
            NegRedAmMxy(i) = RedAmMom(3,NegRedAmMxyIndex(i));
        end
    end
end
NegRedAmMxyNeighb(1) = [];

% Delete entries of Neighbours that are counted twice
NegRedAmMxyNeighb = sort(NegRedAmMxyNeighb);
for i = 1 : length(NegRedAmMxyNeighb) - 1
    if NegRedAmMxyNeighb(i) == NegRedAmMxyNeighb(i + 1);
        NegRedAmMxyNeighb(i) = 0;
    end
end
NegRedAmMxyNeighb = nonzeros(NegRedAmMxyNeighb);

% Dispose the obtained value in the matrix with all Mom to be distributed

if exist('PosRedAmMxy','var')
DistrMom(3,PosRedAmMxyNeighb) = -sum(PosRedAmMxy)/(length(PosRedAmMxyNeighb));
end
if exist('NegRedAmMxy','var')
DistrMom(3,NegRedAmMxyNeighb) = -sum(NegRedAmMxy)/(length(NegRedAmMxyNeighb));
end

% Discard NaN Results
for i = 1:size(DistrMom,2)
    if isnan(DistrMom(3,i))
        DistrMom(3,i) = 0;
    end
end
