function NielsenCriterion1_Elastic(nmomel_x,NegSol,PosSol,NodeCoor,...
    YieldMom,ElastMom,lambda,increment,MomP,Yield,reb_lay_elem,reb_lay) 
%--------------------------------------------------------------------------
%  Purpose:
%     Derive and Plot Linearized Nielsen's Yield Criterion 
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
nmomel_y = nmomel_x;               % Amount of elements on y-axis

%% Create Vertices and Faces to be plotted in the patch
% Indices Vectors to be used in the loop
i = 1:nmomel_x;
ix1 = nmomel_x*i + 1;
ix2 = nmomel_y*i + (nmomel_y + 1);

% Operative Matrices
Vertices1Neg = zeros(nmomel_x*3,3);
Vertices2Neg = Vertices1Neg;
Vertices1Pos = Vertices1Neg;
Vertices2Pos = Vertices1Pos;
Faces = zeros(nmomel_x,3);

%% Generate Vertices and Faces to be plotted in the patch
for v = 1:nmomel_x
    Faces(v,:) = (1:3)+3*(v-1);
    %% Vertices on left side of Bisectrice (1)
    % Negative Solutions
    Vertices1Neg(Faces(v,:),1) = [NodeCoor(1,1) NodeCoor(ix1(v),1) NodeCoor(ix2(v),1)];
    Vertices1Neg(Faces(v,:),2) = [NodeCoor(1,2) NodeCoor(ix1(v),2) NodeCoor(ix2(v),2)];
    Vertices1Neg(Faces(v,:),3) = [0 NegSol(v) NegSol(v+1)];
    % Positive Solutions
    Vertices1Pos(Faces(v,:),1) = [NodeCoor(1,1) NodeCoor(ix1(v),1) NodeCoor(ix2(v),1)];
    Vertices1Pos(Faces(v,:),2) = [NodeCoor(1,2) NodeCoor(ix1(v),2) NodeCoor(ix2(v),2)];
    Vertices1Pos(Faces(v,:),3) = [0 PosSol(v) PosSol(v+1)];
end

%% Vertices on right side of Bisectrice (2)
for v = 1:nmomel_x
    Faces(v,:) = (1:3)+3*(v-1);
    %% Vertices on left side of Bisectrice (1)
    % Negative Solutions
    Vertices2Neg(Faces(v,:),1) = [NodeCoor(end,1) NodeCoor(ix1(v),1) NodeCoor(ix2(v),1)];
    Vertices2Neg(Faces(v,:),2) = [NodeCoor(end,2) NodeCoor(ix1(v),2) NodeCoor(ix2(v),2)];
    Vertices2Neg(Faces(v,:),3) = [0 NegSol(v) NegSol(v+1)];
    % Positive Solutions
    Vertices2Pos(Faces(v,:),1) = [NodeCoor(end,1) NodeCoor(ix1(v),1) NodeCoor(ix2(v),1)];
    Vertices2Pos(Faces(v,:),2) = [NodeCoor(end,2) NodeCoor(ix1(v),2) NodeCoor(ix2(v),2)];
    Vertices2Pos(Faces(v,:),3) = [0 PosSol(v) PosSol(v+1)];
end

%% Plot Yield Criterion and State of Stresses
LinearizedYieldCriterion = figure;
set(LinearizedYieldCriterion,'name',['Postprocessing - Iteration '...
    num2str(increment) ', Rebar layer ' num2str(reb_lay)],'numbertitle','off') ;
p1neg = patch('Faces',Faces,'Vertices',Vertices1Neg,'FaceVertexCData',Vertices1Neg(:,3),'FaceColor','interp');
hold on
p1pos = patch('Faces',Faces,'Vertices',Vertices1Pos,'FaceVertexCData',Vertices1Pos(:,3),'FaceColor','interp');
hold on
p2neg = patch('Faces',Faces,'Vertices',Vertices2Neg,'FaceVertexCData',Vertices2Neg(:,3),'FaceColor','interp');
hold on
p2pos = patch('Faces',Faces,'Vertices',Vertices2Pos,'FaceVertexCData',Vertices2Pos(:,3),'FaceColor','interp');
alpha 0;
pbaspect([1 1 1])
grid on
xlabel('M_x [kNm/m]')
ylabel('M_y [kNm/m]')
zlabel('M_x_y [kNm/m]')
title(['Linearized Yield Criterion - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
set(gca,'FontSize', 19);

% Plot Over-Yielded and Elastic Moments
for iyield = 1:MomP
    if reb_lay_elem(iyield) == reb_lay
    hold on
    if Yield(iyield) == 1
YieldMomSc = scatter3(YieldMom(1,iyield),YieldMom(2,iyield),...
    YieldMom(3,iyield),'MarkerEdgeColor','r','MarkerFaceColor',[1 0 0]);
    else % Yield(MomP) == 0
ElastMomSc = scatter3(ElastMom(1,iyield),ElastMom(2,iyield),...
    ElastMom(3,iyield),'MarkerEdgeColor','g','MarkerFaceColor',[0 1 0]);
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
% set(leg,'Location','NorthWest')

% colormap('jet')
% colorbar

end