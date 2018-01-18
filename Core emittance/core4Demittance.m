function [e4Dcore] = core4Demittance(xg_interp,xp_interp,yg_interp,yp_interp,sigmaxp_interp,sigmayp_interp,sigmaxpyp_interp,int_interp)
% This script reconstructs the projected 2D trace space x-xp and y-yp and
% calculates the core 4D emittance from charge distribution in the ellipse
% area. [Jorge GN]
% >> Input: use interpolated parameters from 'interpolation_PP.m' or 'interpolation_TEM.m'

fastest = true; % Choose fastest (discrete without density plot) or slowest (mapping with density plot) algorithm

if fastest
    [Int_4D,X_4D,Xp_4D,Y_4D,Yp_4D,r0_4D,S_4D] = ...
        Density_4D_discrete(xg_interp,xp_interp,yg_interp,yp_interp, ...
        sigmaxp_interp,sigmayp_interp,sigmaxpyp_interp,...
        int_interp,[100 100 100 100]);
else
    [Int_4D,X_4D,Xp_4D,Y_4D,Yp_4D,r0_4D,S_4D] = ...
        Density_4D(xg_interp,xp_interp,yg_interp,yp_interp, ...
        sigmaxp_interp,sigmayp_interp,sigmaxpyp_interp,...
        int_interp,[100 100 100 100]);
end

% corr = trapz(trapz(intx_interp.*xg_interp.*xp_interp))/trapz(trapz(intx_interp.*xg_interp.^2));  xp_uncorr= xp_interp-xg_interp.*corr;
% corr = trapz(trapz(inty_interp.*yg_interp.*yp_interp))/trapz(trapz(inty_interp.*yg_interp.^2));  yp_uncorr= yp_interp-yg_interp.*corr;
% [Int_4Du,X_4Du,Xp_4Du,Y_4Du,Yp_4Du,r0_4Du,S_4Du] = ...
%     Density_4D(xg_interp,xp_uncorr,yg_interp,yp_uncorr, ...
%     sigmaxp_interp,sigmayp_interp,sigmaxpyp_interp,...
%     int_interp,[100 100 100 100]);

%% Calculate 4D core
r_4D = [X_4D(:) Xp_4D(:) Y_4D(:) Yp_4D(:)]; r_4D = r_4D - repmat(r0_4D,[size(r_4D,1),1]);
figure(891); clf;
[e4Dcore,~,~] = core4D_distr(r_4D,S_4D,Int_4D); 

% r_4Du = [X_4Du(:) Xp_4Du(:) Y_4Du(:) Yp_4Du(:)]; r_4Du = r_4Du - repmat(r0_4Du,[size(r_4Du,1),1]);
% figure(892); clf;
% [e4Dcore,~,~] = core4D_distr(r_4Du,S_4Du,Int_4Du); 