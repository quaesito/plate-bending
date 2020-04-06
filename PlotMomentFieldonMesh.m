function PlotMomentFieldonMesh(NodeCoor,ElemNode,BendMoments,lambda,increment,inflim,suplim) 
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Moment Fields on the mesh
% Synopsis :
%           PlotMomentFieldonMesh(NodeCoor,ElemNode,BendMoments) 
% Variable Description:
%           NodeCoor - The nodal coordinates of the mesh
%           -----> coordinates = [X Y] 
%           ElemNode - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]    
%           BendMoments - Values of Bending moments (Mx, My, Mxy in rows),
%                         per each element (ordered in columns).
%           lambda - Load factor
%           increment - iteration no.
%--------------------------------------------------------------------------

componentMx = BendMoments(1,:);
componentMy = BendMoments(2,:);
componentMxy = BendMoments(3,:);

nelem = length(ElemNode) ;               % number of elements
nnelem = size(ElemNode,2);               % number of nodes per element

% Initialization of the required matrices
X = zeros(nnelem,nelem) ;
Y = zeros(nnelem,nelem) ;
profileMx = zeros(nnelem,length(componentMx));
profileMy = profileMx; profileMxy = profileMx;
nd = zeros(nelem,nnelem);

for ielem = 1:nelem   
     for i = 1:nnelem
     nd(i) = ElemNode(ielem,i);                   % extract connected node for (iel)-th element
     X(i,ielem) = NodeCoor(nd(i),1);              % extract x value of the node
     Y(i,ielem) = NodeCoor(nd(i),2);              % extract y value of the node
     end   
     profileMx(:,ielem) = componentMx(ielem) ;    % extract component value of the node 
     profileMy(:,ielem) = componentMy(ielem) ;    % extract component value of the node 
     profileMxy(:,ielem) = componentMxy(ielem) ;  % extract component value of the node 
end
    
    % Plotting the FEM mesh and profile of the given components
     fMx = figure ;
     set(fMx,'name',['Mx Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'],'numbertitle','off') ;
     fill(X,Y,profileMx)
     daspect([1 1 1])
     axis on ;
     xlabel('x (m)'); ylabel('y (m)');
     colorbar
     colormap(jet(256));
     title(colorbar,'[kNm/m]','FontSize',20);
     caxis([suplim, inflim]);
     title(['Mx Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
     set(gca,'FontSize', 20);
     
     fMy = figure ;
     set(fMy,'name',['My Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'],'numbertitle','off') ;
     fill(X,Y,profileMy)
     daspect([1 1 1])     
     axis on ;
     xlabel('x (m)'); ylabel('y (m)');
     colorbar
     colormap(jet(256));
     title(colorbar,'[kNm/m]','FontSize',20);
     caxis([suplim, inflim]);
     title(['My Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
     set(gca,'FontSize', 20);

     fMxy = figure ;
     set(fMxy,'name',['Mxy Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'],'numbertitle','off') ;
     fill(X,Y,profileMxy)
     daspect([1 1 1])
     axis on ;
     xlabel('x (m)'); ylabel('y (m)');
     colorbar
     colormap(jet(256));
     title(colorbar,'[kNm/m]','FontSize',20);
     caxis([suplim, inflim]);
     title(['Mxy Field on Mesh - Iteration ' num2str(increment) ' : \lambdaP = ' num2str(lambda(increment)) ' kN/m^{2}'])  
     set(gca,'FontSize', 20);

end