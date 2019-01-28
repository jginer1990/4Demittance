function [ex2Dcore,ey2Dcore] = core2Demittance(xg_interp,xp_interp,sigmaxp_interp,intx_interp,yg_interp,yp_interp,sigmayp_interp,inty_interp)
% This script reconstructs the projected 2D trace space x-xp and y-yp and
% calculates the core 2D emittance from charge distribution in the ellipse
% area. [Jorge GN]
% >> Input: use interpolated parameters from 'interpolation_PP.m' or 'interpolation_TEM.m'

% Variant 1: Using correlated x-xp, y-yp
figure(780); clf
subplot(1,2,1);
[Int_2Dxxp,X_2Dxxp,Xp_2Dxxp,r0_2Dxxp,S_2Dxxp] = Density_2D(xg_interp,xp_interp,sigmaxp_interp,intx_interp,[200 200]);
subplot(1,2,2);
[Int_2Dyyp,Y_2Dyyp,Yp_2Dyyp,r0_2Dyyp,S_2Dyyp] = Density_2D(yg_interp,yp_interp,sigmayp_interp,inty_interp,[200 200]); xlabel('y [mm]'); ylabel('yp [mrad]');


% Variant 2: Using uncorrelated xp and yp
% figure(780); clf
% subplot(121);
% corr = trapz(trapz(intx_interp.*xg_interp.*xp_interp))/trapz(trapz(intx_interp.*xg_interp.^2));  xp_uncorr= xp_interp-xg_interp.*corr;
% [Int_2Dxxp,X_2Dxxp,Xp_2Dxxp,r0_2Dxxp,S_2Dxxp] = Density_2D(xg_interp,xp_uncorr,sigmaxp_interp,intx_interp,[200 200]); ylabel('Uncorrelated xp [mrad]');
% subplot(122);
% corr = trapz(trapz(inty_interp.*yg_interp.*yp_interp))/trapz(trapz(inty_interp.*yg_interp.^2));  yp_uncorr= yp_interp-yg_interp.*corr;
% [Int_2Dyyp,Y_2Dyyp,Yp_2Dyyp,r0_2Dyyp,S_2Dyyp] = Density_2D(yg_interp,yp_uncorr,sigmayp_interp,inty_interp,[200 200]); xlabel('y [mm]'); ylabel('Uncorrelated yp [mrad]');


%% Calculate 2D core
rx_2D = [X_2Dxxp(:) Xp_2Dxxp(:)]; rx_2D = rx_2D - repmat(r0_2Dxxp,[size(rx_2D,1),1]);
ry_2D = [Y_2Dyyp(:) Yp_2Dyyp(:)]; ry_2D = ry_2D - repmat(r0_2Dyyp,[size(ry_2D,1),1]);

figure(890); clf;
subplot(2,1,1);
[ex2Dcore,~,~] = core2D_distr(rx_2D,S_2Dxxp,Int_2Dxxp); xlabel('\epsilon^{2D}_{nx,subset} [m rad]');
subplot(2,1,2);
[ey2Dcore,~,~] = core2D_distr(ry_2D,S_2Dyyp,Int_2Dyyp); xlabel('\epsilon^{2D}_{ny,subset} [m rad]');
