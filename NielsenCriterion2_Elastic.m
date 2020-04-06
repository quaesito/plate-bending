function [Mppx,Mnpx,Mppy,Mnpy] = NielsenCriterion2_Elastic(nmomel_x,Mppx,...
    Mnpx,Mppy,Mnpy,YieldMom,ElastMom,lambda,increment,MomP,Yield,reb_lay_elem,reb_lay) 

%--------------------------------------------------------------------------
%  Purpose:
%     Derive and Plot Quadratic Nielsen's Yield Criterion 
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


 %% Quadratic cone %%
% Create Operative Matrices
NegSol1 = zeros(length(MxVec),2); PosSol1 = NegSol1;
NegSol2 = NegSol1; PosSol2 = NegSol1;

% Loop deriving the roots of Mxy for every point
for ix1 = 1:(nmomel_x + 1)
    for iy = 1:(nmomel_y + 1)
syms MxyValues1 MxyValues2
eqnLeftCone = (abs(Mnpx) + MxVec(ix1)).*(abs(Mnpy) + MyVec(iy)) == MxyValues1.^2;
eqnRightCone = (abs(Mppx) - MxVec(ix1)).*(abs(Mppy) - MyVec(iy)) == MxyValues2.^2;
SolLeftCone = solve(eqnLeftCone, MxyValues1);
SolRightCone = solve(eqnRightCone, MxyValues2);
PosSol1(iy,ix1) = real(double(SolLeftCone(1)));
NegSol1(iy,ix1) = real(double(SolLeftCone(2)));
PosSol2(iy,ix1) = real(double(SolRightCone(1)));
NegSol2(iy,ix1) = real(double(SolRightCone(2)));
    end
end

% Derive Envelope of the solutions in each point
EnvelopeNegSol = min(NegSol1,NegSol2);
EnvelopePosSol = max(PosSol1,PosSol2);

%% Plot the obtained Yield Criterion
YieldCriterion = figure ;
set(YieldCriterion,'name',['Postprocessing - Iteration ' num2str(increment)...
    ', Rebar layer ' num2str(reb_lay)],'numbertitle','off') ;
surf(MxVec,MyVec,EnvelopePosSol);
hold on
surf(MxVec,MyVec,EnvelopeNegSol);
alpha 0
pbaspect([1 1 1])
xlabel('M_x [kNm/m]')
ylabel('M_y [kNm/m]')
zlabel('M_x_y [kNm/m]')
title(['Quadratic Yield Criterion - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
hold on

% Plot Over-Yielded and Elastic Moments

for iyield = 1:MomP
    if reb_lay_elem(iyield) == reb_lay
    hold on
    if Yield(iyield) == 1
YieldMomSc = scatter3(YieldMom(1,iyield),YieldMom(2,iyield),YieldMom(3,iyield),'MarkerEdgeColor','r','MarkerFaceColor',[1 0 0]);
    else % Yield(MomP) == 0
ElastMomSc = scatter3(ElastMom(1,iyield),ElastMom(2,iyield),ElastMom(3,iyield),'MarkerEdgeColor','g','MarkerFaceColor',[0 1 0]);
    end
end
end

% Legend
hold on
if exist('YieldMomSc','var') && exist('ElastMomSc','var')
leg = legend([YieldMomSc, ElastMomSc], {'Over-Yielded Moments','Non-yielded Moments'});
set(leg,'Location','NorthWest')
elseif exist('YieldMomSc','var') && exist('ElastMomSc','var') == 0
    leg = legend([YieldMomSc], {'Over-Yielded Moments'});
elseif exist('ElastMomSc','var') == 0
    leg = legend([ElastMomSc], {'Non-yielded Moments'});
end

% % Plot Yield Criterion with no envelope - Intersection of Planes
% YieldCriterion = figure ;
% set(YieldCriterion,'name','Postprocessing','numbertitle','off') ;
% surf(MxVec,MyVec,PosSol1);
% hold on 
% surfl(MxVec,MyVec,NegSol1);
% hold on
% surf(MxVec,MyVec,PosSol2);
% hold on
% surf(MxVec,MyVec,NegSol2);
% xlabel('M_x')
% ylabel('M_y')
% zlabel('M_x_y')
% title('Yield Criterion') ;

end

