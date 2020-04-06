function [Mpp,Mnp] = PlasticMoments(fcu,fy,B,H,C,phis,ns,phi,n)

%%%----------------------------------------------------------
%  Purpose:
%     Evaluate Maximum Sagging and Hogging Resisting Moment according to
%     code of practice HK Concrete 13
%
%  Synopsis:
%     [Mr] = ULSbending(fcu,B,H,C,phis,ns,phi,n)
%
%  Variable Description:
%     fcu       Concrete compressive strength [MPa]
%     B         Slab width [mm]
%     H         Slab height [mm]
%     C         Distance of rebars from bottom/top of the slab [mm]
%     phis      Rebars diameter on the lower edge [mm]
%     phi       Rebars diameter on the upper edge [mm]
%     ns        Amount of rebars on the lower edge [-]
%     n         Amount of rebars on the upper edge [-]
%--------------------------------------------------------------------------
% Michele De Filippo
% Department of Civil Engineering
% The Hong Kong University of Science and Technology
% Latest revision: June 2017
%--------------------------------------------------------------------------

%% Concrete characteristics %%
epscu = 3.5E-3;        % Concrete ultimate strain [-]
D = H - C;             % Distance of rebars from top of the slab [mm]
beta = 0.9;            % Stress-block height percentage on neutral axis [-]
k = beta/2;            % Compression force distance from top of the beam percentage on neutral axis [-]
            
%% Rebar characteristics %%
epss = 1.86E-3;           % Yielding Strain [-]
Es = fy/epss;             % Rebar Young's Modulus [MPa]
Ass = (pi*phis^2/4)*ns;   % Rebar area on the lower edge [mm^2]
As = (pi*phi^2/4)*n;      % Rebar area on the upper edge [mm^2]

B = B*1000;               % Convert B from m unit to mm unit

%% Create Iteration vector for xn
no = 1000;                     % Amount of iterations
xn_sag_it = zeros(no+1,1);     % Create starting vector
xn_hog_it = xn_sag_it;
eq_sag = xn_sag_it;
eq_hog = eq_sag;
Mr_sag = zeros(no,1);
Mr_hog = Mr_sag;
delta_xn = H/no;          % Define step xn
% Consider cross-section upside down for hogging moments
As_hog = Ass;
Ass_hog = As;

%% HK Concrete 13 Calculation %%

%% 1) DOUBLE REINFORCED CASE %%
if As ~= 0 && Ass ~= 0
    
% SAGGING RESISTING MOMENT %
for i = 2:no+1
    % xn from upper edge
    xn_sag_it(i,1) = xn_sag_it(i-1,1) + delta_xn;  % Next step of iteration
    if xn_sag_it(i) >= C      % xn exceeds c
        if epscu*(xn_sag_it(i,1) - C)/xn_sag_it(i,1) >= epss % upper rebars are yielding
            eq_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1) + As*0.87*fy -...
                Ass*0.87*fy)/1000;
            Mr_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1)*(D - k*xn_sag_it(i,1)) + ...
                As*0.87*fy*(D - C))*1E-6;
        else % upper rebars not yielding
            eq_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1) + ...
                As*(epscu*(xn_sag_it(i,1) - C)/xn_sag_it(i,1)*Es) - ...
                Ass*0.87*fy)/1000;
            Mr_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1)*(D - k*xn_sag_it(i,1)) + ...
                As*0.87*(epscu*(xn_sag_it(i,1) - C)/xn_sag_it(i,1)*Es)*(D - C))*1E-6;
         end
    else                 % xn within c
        if epscu*(C - xn_sag_it(i,1))/xn_sag_it(i,1) >= epss % upper rebars are yielding
            eq_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1) - As*0.87*fy -...
                Ass*0.87*fy)/1000 ;
            Mr_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1)*(D - k*xn_sag_it(i,1)) -...
                As*0.87*fy*(D - C))*1E-6;
        else % upper rebars not yielding
            eq_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1) - Ass*0.87*fy -...
                As*0.87*(epscu*(C - xn_sag_it(i,1))/xn_sag_it(i,1)*Es))/1000;
            Mr_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1)*(D - k*xn_sag_it(i,1)) -...
               As*0.87*(epscu*(C - xn_sag_it(i,1))/xn_sag_it(i,1)*Es)*(D - C))*1E-6;
        end
    end
end
eq_sag(1) = []; Mr_sag(1) = []; % Discard 0 values first entries

% HOGGING RESISTING MOMENT %
% Consider cross-section upside down
for j = 2:no+1
    % xn from upper edge
    xn_hog_it(j,1) = xn_hog_it(j-1,1) + delta_xn;  % Next step of iteration
    if xn_hog_it(j) >= C      % xn exceeds c
        if epscu*(xn_hog_it(j,1) - C)/xn_hog_it(j,1) >= epss % upper rebars are yielding
            eq_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1) + As_hog*0.87*fy -...
                Ass_hog*0.87*fy)/1000;
            Mr_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1)*(D - k*xn_hog_it(j,1)) + ...
                As_hog*0.87*fy*(D - C))*1E-6;
        else % upper rebars not yielding
            eq_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1) + ...
                As_hog*(epscu*(xn_hog_it(j,1) - C)/xn_hog_it(j,1)*Es) - ...
                Ass_hog*0.87*fy)/1000;
            Mr_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1)*(D - k*xn_hog_it(j,1)) + ...
                As_hog*0.87*(epscu*(xn_hog_it(j,1) - C)/xn_hog_it(j,1)*Es)*(D - C))*1E-6;
        end
    else                 % xn within c
        if epscu*(C - xn_hog_it(j,1))/xn_hog_it(j,1) >= epss % upper rebars are yielding
            eq_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1) - As_hog*0.87*fy -...
                Ass_hog*0.87*fy)/1000 ;
            Mr_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1)*(D - k*xn_hog_it(j,1)) -...
                As_hog*0.87*fy*(D - C))*1E-6;
        else % upper rebars not yielding
            eq_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1) - Ass_hog*0.87*fy -...
                As_hog*0.87*(epscu*(C - xn_hog_it(j,1))/xn_hog_it(j,1)*Es))/1000;
            Mr_hog(j) = (0.45*fcu*beta*B*xn_hog_it(j,1)*(D - k*xn_hog_it(j,1)) -...
               As_hog*0.87*(epscu*(C - xn_hog_it(j,1))/xn_hog_it(j,1)*Es)*(D - C))*1E-6;
        end
    end
end
eq_hog(1) = []; Mr_hog(1) = [];    % Discard 0 values first entries

[eq_sag,I_sag] = min(abs(eq_sag)); % Find index of absolute min value
eq_sag = [];                       % Discard its value
Mpp = Mr_sag(I_sag);               % Extract sagging resisting moment

[eq_hog,I_hog] = min(abs(eq_hog)); % Find index of absolute min value
eq_hog = [];                       % Discard its value
Mnp = -Mr_hog(I_hog);              % Extract sagging resisting moment

%% 2) SINGLE REINFORCED CASE - ONLY LOWER REBARS %%
elseif  As == 0 && Ass ~= 0
    % SAGGING RESISTING MOMENT % - Lower rebars supposed to ALWAYS yield
for i = 2:no+1
    % xn from upper edge
    xn_sag_it(i,1) = xn_sag_it(i-1,1) + delta_xn;  % Next step of iteration
    eq_sag(i) = (0.45*fcu*beta*B*xn_sag_it(i,1) - Ass*0.87*fy)/1000;
    Mr_sag(i) = 0.45*fcu*beta*B*xn_sag_it(i,1)*(D - k*xn_sag_it(i,1)/2)*1E-6;
end
    eq_sag(1) = []; Mr_sag(1) = [];    % Discard 0 values first entries   

    [eq_sag,I_sag] = min(abs(eq_sag)); % Find index of absolute min value
    eq_sag = [];                       % Discard its value
    Mpp = Mr_sag(I_sag);               % Extract sagging resisting moment
    
    % HOGGING RESISTING MOMENT % - assumed to be null
    Mnp = 0;
    
%% 3) SINGLE REINFORCED CASE - ONLY UPPER REBARS %%
elseif  As ~= 0 && Ass == 0
    % HOGGING RESISTING MOMENT % - Upper rebars supposed to ALWAYS yield
    for j = 2:no+1
    % xn from upper edge
    xn_hog_it(j,1) = xn_hog_it(j-1,1) + delta_xn;  % Next step of iteration
    eq_hog(j) = (0.45*fcu*beta*B*(H - xn_hog_it(j,1)) - As*0.87*fy)/1000;
    Mr_hog(j) = 0.45*fcu*beta*B*(H - xn_hog_it(j,1))*...
        (D - k*(H - xn_hog_it(j,1))/2)*1E-6;
    end
    
    eq_hog(1) = []; Mr_hog(1) = [];    % Discard 0 values first entries   

    [eq_hog,I_hog] = min(abs(eq_hog)); % Find index of absolute min value
    eq_hog = [];                       % Discard its value
    Mnp = -Mr_hog(I_hog);               % Extract sagging resisting moment
    
    % SAGGING RESISTING MOMENT % - assumed to be null
    Mpp = 0;
end
end
       
