function [kk,ff]=constraints(kk,ff,bcdof)

%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [K]{x}={F}
%
%  Synopsis:
%     [kk,ff]=constraints(kk,ff,bcdof)
%
%  Variable Description:
%     kk - system matrix before applying constraints 
%     ff - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%-----------------------------------------------------------
 
 n=length(bcdof);
 sdof=size(kk);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
       kk(j,c)=0;
    end

    kk(c,c)=1;
    ff(c)=0;
    % here we are making all elements zero
    % from those rows and columns corresponds to the boundary DOFs except kk(DOF,DOF) element.
 end

