function [f] = Force(nnel,shape,P,ielem,load_elem,loadType)

%--------------------------------------------------------------------------
% Purpose :
%         Determines the force vector for the element 'iel'
% Synopsis :
%          [f] = Force(nnel,shape,P) 
% Variable Description:
%           f - element force vector
%           nnel - number of nodes per element
%           shape - Shape functions for the element
%           P - applied transverse pressure
%--------------------------------------------------------------------------

fef = shape*P ; % distributed pressure is converted into nodal forces, ...
...forces is applicable to transverse direction only i.e. related to w only
% Case 0) Uniform Transverse Load
if loadType == 0 
for i = 1:nnel
    i1=(i-1)*3+1;  
    i2=i1+1;
    i3=i2+1;
    f(i1,1) = fef(1) ;
    f(i2,1) = 0 ;
    f(i3,1) = 0 ;
end
% Case 1) One Patch load
elseif loadType == 1
    if load_elem(ielem) == 1
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    else
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = 0 ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    end
% Case 2) Two Patch loads
elseif loadType == 2
    if load_elem(ielem) == 1
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    elseif load_elem(ielem) == 2
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    else
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = 0 ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    end
% Case 3) Three Patch Loads
elseif loadType == 3
    if load_elem(ielem) == 1
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    elseif load_elem(ielem) == 2
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    elseif load_elem(ielem) == 3
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    else
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = 0 ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    end
% Case 4) Four Patch Loads
elseif loadType == 4
    if load_elem(ielem) == 1
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    elseif load_elem(ielem) == 2
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    elseif load_elem(ielem) == 3
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    elseif load_elem(ielem) == 4
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = fef(1) ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    else
        for i = 1:nnel
        i1=(i-1)*3+1;  
        i2=i1+1;
        i3=i2+1;
        f(i1,1) = 0 ;
        f(i2,1) = 0 ;
        f(i3,1) = 0 ;
        end
    end
    end
end
