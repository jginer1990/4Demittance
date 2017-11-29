%% 
% This script reconstructs the projected 2D trace space x-xp and y-yp and
% calculates the core 2D emittance from charge distribution in the ellipse
% area. [Jorge GN]
% >> Run this script after 'phasespace_shear_interpolation.m'

if strcmp(target,'TEM')
    int_interp = sqrt(intx_interp.*inty_interp);
end


[Int_4D,X_4D,Xp_4D,Y_4D,Yp_4D,r0_4D,S_4D] = ...
    Density_4D(xg_interp,xp_interp,yg_interp,yp_interp, ...
    sigmaxp_interp,sigmayp_interp,sigmaxpyp_interp,...
    int_interp,[100 100 100 100]);



%% Calculate 4D core
r_4D = [X_4D(:) Xp_4D(:) Y_4D(:) Yp_4D(:)]; r_4D = r_4D - repmat(r0_4D,[size(r_4D,1),1]);

figure(891); clf;
[e4Dcore,~,~] = core4D_distr(r_4D,S_4D,Int_4D); xlabel('\epsilon^{4D}_{n,subset} [m rad]');
