function PlotPrincipalMomentFieldonMesh(NodeCoor,ElemNode,PrincipalBendMoments,lambda,increment,inflim,suplim) 
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Principal Moment Fields on the mesh
% Synopsis :
%           PlotPrincipalMomentFieldonMesh(NodeCoor,ElemNode,Principal Moments) 
% Variable Description:
%           NodeCoor - The nodal coordinates of the mesh
%           -----> coordinates = [X Y] 
%           ElemNode - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]    
%           PrincipalBendMoments - Values of Bending moments (M1, M2, M12 in rows),
%                         per each element (ordered in columns).
%           lambda - Load factor
%           increment - iteration no.
%--------------------------------------------------------------------------

componentM1 = PrincipalBendMoments(1,:);
componentM2 = PrincipalBendMoments(2,:);
% componentM12 = PrincipalBendMoments(3,:);

nelem = length(ElemNode) ;               % number of elements
nnelem = size(ElemNode,2);               % number of nodes per element

% Initialization of the required matrices
X = zeros(nnelem,nelem) ;
Y = zeros(nnelem,nelem) ;
profileM1 = zeros(nnelem,length(componentM1));
profileM2 = profileM1; profileM12 = profileM1;
nd = zeros(nelem,nnelem);

for ielem = 1:nelem   
     for i = 1:nnelem
     nd(i) = ElemNode(ielem,i);         % extract connected node for (iel)-th element
     X(i,ielem) = NodeCoor(nd(i),1);    % extract x value of the node
     Y(i,ielem) = NodeCoor(nd(i),2);    % extract y value of the node
     end   
     profileM1(:,ielem) = componentM1(ielem) ;           % extract component value of the node 
     profileM2(:,ielem) = componentM2(ielem) ;           % extract component value of the node 
    % profileM12(:,ielem) = componentM12(ielem) ;         % extract component value of the node 
end
    
    % Plotting the FEM mesh and profile of the given components
     fM1 = figure ;
     set(fM1,'name',['M1 Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'],'numbertitle','off') ;
     fill(X,Y,profileM1)
     daspect([1 1 1])
     axis on ;
     xlabel('x (m)'); ylabel('y (m)');
     colorbar
     colormap(jet(256));
     caxis([suplim, inflim]);
     set(gca,'FontSize', 20);
     title(['M1 Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
     
     fM2 = figure ;
     set(fM2,'name',['M2 Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'],'numbertitle','off') ;
     fill(X,Y,profileM2)
     daspect([1 1 1])
     axis on ;
     xlabel('x (m)'); ylabel('y (m)');
     colorbar
     colormap(jet(256));
     caxis([suplim, inflim]);
     set(gca,'FontSize', 20);
     title(['M2 Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
     
%      fM12 = figure ;
%      set(fM12,'name','Postprocessing','numbertitle','off') ;
%      fill(X,Y,profileM12)
%      axis off ;
%      SetColorbar
%      title('M12 Field on Mesh')
end