function bcdof = IntBoundaryCondition(typeBC,coordinates,dir,...
    coor_x1,coor_y1,coor_x2,coor_y2,nelem_x,color)

%--------------------------------------------------------------------------
%   Purpose:
%           To determine the boundary conditions degree of freedom
%   Synopsis:
%           bcdof = BoundaryCondition(typeBC,coordinates,loadstep)
%
%   Variable Description:
%           bcdof - boundary condition degree of freedom
%                   dofs are (w, thetax, thetay)
%           typeBC - string which gives the state of boundary condition
%           coordinates - geometric coordinates of nodes
%           
%--------------------------------------------------------------------------
if dir == 'x'
    L = find(coordinates(:,2) == coor_y1) ; % at y = coor (along X-axes)
elseif dir == 'y'
    L = find(coordinates(:,1) == coor_x1) ; % at x = coor (along Y-axes)
elseif dir == 'point'
    Ly = find(round(coordinates(:,2),3) == coor_y1) ; % at y = coor (along X-axes)
    Lx = find(round(coordinates(:,1),3) == coor_x1) ; % at y = coor (along X-axes)
    L = intersect(Lx,Ly);                    % find point support
elseif dir == 'patch'
    Ly1 = find(round(coordinates(:,2),3) == coor_y1) ; % at y = coor (along X-axes)
    Lx1 = find(round(coordinates(:,1),3) == coor_x1) ; % at x = coor (along Y-axes)
    Ly2 = find(round(coordinates(:,2),3) == coor_y2) ; % at y = coor (along X-axes)
    Lx2 = find(round(coordinates(:,1),3) == coor_x2) ; % at x = coor (along Y-axes)
    % Find 4 extreme nodes of the patch (COUNTERCLOCKWISE!!!)
    edgesL = intersect(Lx1,Ly1); 
    edgesL(end+1) = intersect(Lx2,Ly1); 
    edgesL(end+1) = intersect(Lx2,Ly2);  
    edgesL(end+1) = intersect(Lx1,Ly2); 
    % Find all nodes at the edges 
    edgesL = [edgesL,edgesL(1):edgesL(2),...
                     edgesL(2):nelem_x+1:edgesL(3),...
                     edgesL(4):edgesL(3),...
                     edgesL(1):nelem_x+1:edgesL(4)];
    % Find all internal nodes
    intL = [];
    for i = edgesL(1)+nelem_x+2:nelem_x+1:edgesL(4)-nelem_x % Iterate row by row and add all nodes in each row
        intL = [intL,i:i+(edgesL(2)-edgesL(1))-2];
    end
    edgesL = unique(edgesL(:));   % Delete Double Entries
    intL = unique(intL(:));   % Delete Double Entries
end

 
if strcmp(typeBC,'s-s-s-s')
    m = length(L) ;
         if dir == 'x'
            disp(['Simply supported along x-direction at ', num2str(coor_y1), ' m'])
         elseif dir == 'y'
            disp(['Simply supported along y-direction at ', num2str(coor_x1), ' m'])
         end
         
        dofL = zeros(1,2*m) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL(i2) = 3*L(i)-2 ;   % constraining w along x axis
        end   
        bcdof = dofL ;
        
        ss = plot(coordinates(L,1),coordinates(L,2),'--');
        legend(ss,'Support');
        legend(gca,'off');
        legend('show')
end
if  strcmp(typeBC,'c-c-c-c')
        m = length(L) ;
         if dir == 'x'
            disp(['Clamped along x-direction at ', num2str(coor_y1), ' m'])
         elseif dir == 'y'
            disp(['Clamped along y-direction at ', num2str(coor_x1), ' m'])
         end
         
        dofL = zeros(1,2*m) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL(i1) = 3*L(i)-1;  % constraining thetax 
            dofL(i1) = 3*L(i);    % constraining thetay 
            dofL(i2) = 3*L(i)-2 ; % constraining w 
        end   
        bcdof = dofL ;
        
        cc = plot(coordinates(L,1),coordinates(L,2),'x');
        legend(cc,'Clamped support');
        legend(gca,'off');
        legend('show')

end

if strcmp(typeBC,'points')
        m = length(L) ;
            disp(['Point Support at (',num2str(coor_x1),',',num2str(coor_y1),')'])
        dofL = zeros(1,2*m) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofL(i2) = 3*L(i)-2 ; % constraining w 
        end   

        bcdof = dofL ;
        
        ps = plot(coor_x1,coor_y1,'o');
        legend(ps,'Point support');
        legend(gca,'off');
        legend('show')
end

if strcmp(typeBC,'patchs, s-s-s-s')
            disp(['Patch Support spanning from (',...
                num2str(coor_x1),',',num2str(coor_y1),') to (',...
                num2str(coor_x2),',',num2str(coor_y2),')']);
        m = length(edgesL) ;
        % BCs at edges of patch
        dofedgesL = zeros(1,2*m) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofedgesL(i2) = 3*edgesL(i)-2 ; % constraining w along x axis
        end
        n = length(intL);
        % BCs within the patch
        dofintL = zeros(1,2*n) ;
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
%             dofintL(i1) = 3*intL(i)-1;  % constraining thetax 
%             dofintL(i1) = 3*intL(i);    % constraining thetay 
            dofintL(i2) = 3*intL(i)-2 ; % constraining w along 
        end   
        dofL = [dofedgesL,dofintL];
        dofL = unique(dofL);

        bcdof = dofL ;
        
        patch_supp = fill([coordinates(intersect(Lx1,Ly1),1),... % Layer 1
           coordinates(intersect(Lx2,Ly1),1),...
           coordinates(intersect(Lx2,Ly2),1),...
           coordinates(intersect(Lx1,Ly2),1)],...
           [coordinates(intersect(Lx1,Ly1),2),...
           coordinates(intersect(Lx2,Ly1),2),...
           coordinates(intersect(Lx2,Ly2),2),...
           coordinates(intersect(Lx1,Ly2),2)],color,'facealpha',.2);
        legend(patch_supp,'Patch support');
        legend(gca,'off');
        legend('show')
end

if strcmp(typeBC,'patchs, c-c-c-c')
            disp(['Patch Support spanning from (',...
                num2str(coor_x1),',',num2str(coor_y1),') to (',...
                num2str(coor_x2),',',num2str(coor_y2),')']);
        m = length(edgesL) ;
        % BCs at edges of patch
        dofedgesL = zeros(1,2*m) ;
        for i = 1:m
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofedgesL(i1) = 3*edgesL(i)-1;  % constraining thetax 
            dofedgesL(i1) = 3*edgesL(i);    % constraining thetay 
            dofedgesL(i2) = 3*edgesL(i)-2 ; % constraining w along
        end
        n = length(intL);
        % BCs within the patch
        dofintL = zeros(1,2*n) ;
        for i = 1:n
            i1 = 2*(i-1)+2 ;
            i2 = i1-1 ;
            dofintL(i1) = 3*intL(i)-1;  % constraining thetax 
            dofintL(i1) = 3*intL(i);    % constraining thetay 
            dofintL(i2) = 3*intL(i)-2 ; % constraining w along 
        end   
        dofL = [dofedgesL,dofintL];
        dofL = unique(dofL);
        
        bcdof = dofL ;
        
        %% Plot Layer Area          
        patch_clamp = fill([coordinates(intersect(Lx1,Ly1),1),... % Layer 1
           coordinates(intersect(Lx2,Ly1),1),...
           coordinates(intersect(Lx2,Ly2),1),...
           coordinates(intersect(Lx1,Ly2),1)],...
           [coordinates(intersect(Lx1,Ly1),2),...
           coordinates(intersect(Lx2,Ly1),2),...
           coordinates(intersect(Lx2,Ly2),2),...
           coordinates(intersect(Lx1,Ly2),2)],color,'facealpha',.2);
       legend(patch_clamp,'Patch clamp');
        legend(gca,'off');
        legend('show')
hold on
end
bcdof(any(bcdof == 0)) = []; % Remove any 0 entry
end