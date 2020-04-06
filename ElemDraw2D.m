function ElemDraw2D(ElemX,ElemY,plotpar,nspline,elnum,nodes,show)
%        ElemDraw2D(ElemX,ElemY,plotpar)
%		 ElemDraw2D(ElemX,ElemY,plotpar,nspline)
%		 ElemDraw2D(ElemX,ElemY,plotpar,nspline,elnum)
%-------------------------------------------------------------
% Draws the unformed element mesh for a number of 2D elements
% of the same type.
%
% INPUT  
%	ElemX		ElementCoordinates dim = nelem x nnode_elem
%	ElemY		Nelem indicates the number of items and nnode_elem
%               indicates number of nodes per element. The first
%               coordinates shall specify the nodes, then
%               follows coordinates of other nodes.
%	plotpar		Specifies parameters for the lines, for example
%               Plotpar = '-rx' If nodes have to be plotted.
%               The node statement MUST be in third position.
%               In the example it is 'x'.
%	nspline		Specifies how many points to draw along
%               each element side
%   elnum		Optional entry of item numbers
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: June 2017
%--------------------------------------------------------------------------

a = size(ElemX) ; b = size(ElemY) ;
 
if (a-b) == [0 0]
    nelem = a(1) ;
	nnode_elem = a(2) ;
else
   error('The size of the coordinates is incorrect');
end

NaNvec = NaN*ones(nelem,1) ; % Add NaN vector to delete a column when necessery

%% Triangular Elements 
if nnode_elem == 3 % CST-Element (Constant STrain)
    Xs = [ElemX(:,1) ElemX(:,2) ElemX(:,3) ElemX(:,1) NaNvec]' ;
    Ys = [ElemY(:,1) ElemY(:,2) ElemY(:,3) ElemY(:,1) NaNvec]' ;
    
	Xt = Xs ; Yt = Ys ;
    
elseif nnode_elem == 6 % LST-Element (Linear STrain)
    Xt = [ElemX(:,1) ElemX(:,4) ElemX(:,2) ElemX(:,5) ElemX(:,3) ElemX(:,6) ElemX(:,1) ] ;
    Yt = [ElemY(:,1) ElemY(:,4) ElemY(:,2) ElemY(:,5) ElemY(:,3) ElemY(:,6) ElemY(:,1) ] ;
	
	xiop = linspace(0,1,nspline) ;  % Isoparametric evaluation points starting from 0
	xined = linspace(1,0,nspline) ; % Isoparametric evaluation points starting from 1
	
	%-- Side 1 ----
	N1s1 = xined.*(2*xined-1) ; % Evaluation of N1 along side 1
	N2s1 = xiop .*(2*xiop -1) ; % Evaluation of N2 along side 1
	N4s1 = 4*xiop.*xined	  ; % Evaluation of N4 along side 1
	
	%-- Side 2 ----
	N2s2 = xined.*(2*xined-1) ; % Evaluation of N2 along side 2
	N3s2 = xiop .*(2*xiop -1) ; % Evaluation of N3 along side 2
	N5s2 = 4*xiop.*xined	  ; % Evaluation of N5 along side 2
	
	%-- Side 3 ----
	N3s3 = xined.*(2*xined-1) ; % Evaluation of N3 along side 3
	N1s3 = xiop .*(2*xiop -1) ; % Evaluation of N1 along side 3
	N6s3 = 4*xiop.*xined	  ; % Evaluation of N6 along side 3
	
	Xs = [Xt(:,1)*N1s1 + Xt(:,2)*N4s1 + Xt(:,3)*N2s1   NaNvec... side 1
		  Xt(:,3)*N2s2 + Xt(:,4)*N5s2 + Xt(:,5)*N3s2   NaNvec... side 2
		  Xt(:,5)*N3s3 + Xt(:,6)*N6s3 + Xt(:,7)*N1s3   NaNvec]' ; % side 3 ;

	Ys = [Yt(:,1)*N1s1 + Yt(:,2)*N4s1 + Yt(:,3)*N2s1   NaNvec... side 1
		  Yt(:,3)*N2s2 + Yt(:,4)*N5s2 + Yt(:,5)*N3s2   NaNvec... side 2
		  Yt(:,5)*N3s3 + Yt(:,6)*N6s3 + Yt(:,7)*N1s3   NaNvec]' ; % side 3 ;
	  
elseif nnode_elem == 10 % 10-nodes Elements QST (Quadratic STrain)
    Xt = [ElemX(:,1) ElemX(:,4) ElemX(:,5) ElemX(:,2) ElemX(:,6) ElemX(:,7) ElemX(:,3) ElemX(:,8) ElemX(:,9) ElemX(:,1) NaNvec] ;
    Yt = [ElemY(:,1) ElemY(:,4) ElemY(:,5) ElemY(:,2) ElemY(:,6) ElemY(:,7) ElemY(:,3) ElemY(:,8) ElemY(:,9) ElemY(:,1) NaNvec] ;
	
	xiop = linspace(0,1,nspline) ;  % Isoparametric evaluation points starting from 0
	xined = linspace(1,0,nspline) ; % Isoparametric evaluation points starting from 1
	
	disp('Note that ElemDraw2D has not yet been tested for QST elements') ;

	%-- Side 1 ----
	N1s1 = 0.5*xined.*(3*xined-1).*(3*xined-2) ; % Evaluation of N1 along side 1
	N2s1 = 0.5*xiop .*(3*xiop -1).*(3*xiop -2) ; % Evaluation of N2 along side 1
	N4s1 = 4.5*xiop .*xined.*(3*xined-1)	   ; % Evaluation of N4 along side 1
	N5s1 = 4.5*xined.*xiop .*(3*xiop -1)	   ; % Evaluation of N5 along side 1
	
	%-- Side 2 ----
	N2s2 = N1s1 ; N3s2 = N2s1 ; N6s2 = N4s1 ; N7s2 = N5s1 ;

	%-- Side 3 ----
	N3s3 = N2s2 ; N1s3 = N3s2 ; N8s3 = N6s2 ; N9s3 = N7s2 ;
	
	Xs = [Xt(:,1)*N1s1 + Xt(:,2)*N4s1 + Xt(:,3)*N5s1 + Xt(:,4)*N2s1   NaNvec... side 1
		  Xt(:,4)*N2s2 + Xt(:,5)*N6s2 + Xt(:,6)*N7s2 + Xt(:,7)*N3s2   NaNvec... side 2
		  Xt(:,7)*N3s3 + Xt(:,8)*N8s3 + Xt(:,9)*N9s3 + Xt(:,10)*N1s3  NaNvec]' ; % side 3 ;

	Ys = [Yt(:,1)*N1s1 + Yt(:,2)*N4s1 + Yt(:,3)*N5s1 + Yt(:,4)*N2s1   NaNvec... side 1
		  Yt(:,4)*N2s2 + Yt(:,5)*N6s2 + Yt(:,6)*N7s2 + Yt(:,7)*N3s2   NaNvec... side 2
		  Yt(:,7)*N3s3 + Yt(:,8)*N8s3 + Yt(:,9)*N9s3 + Yt(:,10)*N1s3  NaNvec]' ; % side 3 ;

%% Square Elements 
elseif nnode_elem == 4 % Four-nodes element
    Xs = [ElemX(:,1) ElemX(:,2) ElemX(:,3) ElemX(:,4) ElemX(:,1) NaNvec]' ;
    Ys = [ElemY(:,1) ElemY(:,2) ElemY(:,3) ElemY(:,4) ElemY(:,1) NaNvec]' ;
	
	Xt = Xs ; Yt = Ys ;
elseif nnode_elem == 8 % Eight-nodes element (serendipity element)

    Xt = [ElemX(:,1) ElemX(:,5) ElemX(:,2) ElemX(:,6) ElemX(:,3) ElemX(:,7) ElemX(:,4) ElemX(:,8) ElemX(:,1) NaNvec] ;
    Yt = [ElemY(:,1) ElemY(:,5) ElemY(:,2) ElemY(:,6) ElemY(:,3) ElemY(:,7) ElemY(:,4) ElemY(:,8) ElemY(:,1) NaNvec] ;
	
	xiop = linspace(-1,1,nspline) ;  % Isoparametric evaluation starting from -1
	xined = linspace(1,-1,nspline) ; % Isoparametric evaluation starting from +1
	
	%-- Side 1 ----
	N5s1 = 0.5*(1-xiop.^2)*2		  ; % Evaluation of N5 along side 1
	N1s1 = 0.25*(1-xiop)*2 - 0.5*N5s1 ; % Evaluation of N1 along side 1
	N2s1 = 0.25*(1+xiop)*2 - 0.5*N5s1 ; % Evaluation of N2 along side 1
	
	%-- Side 2 ----
	N6s2 = 0.5*2*(1-xiop.^2)		  ; % Evaluation of N6 along side 2
	N2s2 = 0.25*2*(1-xiop) - 0.5*N6s2 ; % Evaluation of N2 along side 2
	N3s2 = 0.25*2*(1+xiop) - 0.5*N6s2 ; % Evaluation of N3 along side 2
	
	%-- Side 3 ----
	N7s3 = 0.5*(1-xined.^2)*2		   ; % Evaluation of N7 along side 3
	N3s3 = 0.25*(1+xined)*2 - 0.5*N7s3 ; % Evaluation of N3 along side 3
	N4s3 = 0.25*(1-xined)*2 - 0.5*N7s3 ; % Evaluation of N4 along side 3

	%-- Side 4 ----
	N8s4 = 0.5*2*(1-xined.^2)		   ; % Evaluation of N8 along side 4
	N4s4 = 0.25*2*(1+xined) - 0.5*N8s4 ; % Evaluation of N4 along side 4
	N1s4 = 0.25*2*(1-xined) - 0.5*N8s4 ; % Evaluation of N1 along side 4
	
	Xs = [Xt(:,1)*N1s1 + Xt(:,2)*N5s1 + Xt(:,3)*N2s1   NaNvec... side 1
		  Xt(:,3)*N2s2 + Xt(:,4)*N6s2 + Xt(:,5)*N3s2   NaNvec... side 2
		  Xt(:,5)*N3s3 + Xt(:,6)*N7s3 + Xt(:,7)*N4s3   NaNvec... side 3
		  Xt(:,7)*N4s4 + Xt(:,8)*N8s4 + Xt(:,9)*N1s4   NaNvec]' ; % side 4 ;

	Ys = [Yt(:,1)*N1s1 + Yt(:,2)*N5s1 + Yt(:,3)*N2s1   NaNvec... side 1
		  Yt(:,3)*N2s2 + Yt(:,4)*N6s2 + Yt(:,5)*N3s2   NaNvec... side 2
		  Yt(:,5)*N3s3 + Yt(:,6)*N7s3 + Yt(:,7)*N4s3   NaNvec... side 3
		  Yt(:,7)*N4s4 + Yt(:,8)*N8s4 + Yt(:,9)*N1s4   NaNvec]' ; % side 4 ;
else
   error('This element type is not supported');
end

Xk = Xt' ; Yk = Yt' ;
Xnode = Xk(:) ;
Ynode = Yk(:) ;

%----- Plotting -------------------------------
f = figure ;
set(f,'name','Preprocessing - Mesh','numbertitle','off') ;axis('equal') ; hold on ;
plot(Xs(:),Ys(:),plotpar) ; % Elements plot
title('Finite Element Mesh','fontsize',14) ;

if size(plotpar,2) == 3 % Nodes plot
	plot(Xnode,Ynode,plotpar(3)) ;
end
	
if show == 1 % Element Number plots
	x0 = sum(ElemX')/nnode_elem ;
	y0 = sum(ElemY')/nnode_elem ; % Center coordinates of the elements
    for i = 1:nelem
        h = text(x0(i),y0(i),int2str(elnum(i)),'EdgeColor','red');
        set(h,'fontsize',14);
    end 

elseif show == 2
    % Element Number plots
	x0 = sum(ElemX')/nnode_elem ;
	y0 = sum(ElemY')/nnode_elem ; % Center coordinates of the elements
    for i = 1:nelem
        h = text(x0(i),y0(i),int2str(elnum(i)),'EdgeColor','red');
        set(h,'fontsize',14);
    end

    for inode = 1:nodes(end,1)
    % Node Number Plots
    h = text(nodes(inode,2),nodes(inode,3),int2str(nodes(inode,1)));
    set(h,'fontsize',14);
    end
end

xlabel('x (m)'); ylabel('y (m)');
axis off
 hold off 
end
%--------------------------end--------------------------------
