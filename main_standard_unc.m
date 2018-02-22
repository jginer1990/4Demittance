
%%This script selects particular files to analyse, loads an image, and sends it to be evaluated 
%%by the phase space routine. 

%% Initialization
clear
close all
%warning off;
% dbstop if all error

addpath('./Analysis_Tools/')
addpath('./Core emittance/')

%% Import parameters and files

Input_data_standard;

load(filename)

A = Z_image;

BG = mean(mean(A(1:30,1:30))); A=A-BG;

xvec = pxconv*(1:size(A,2));
yvec = pxconv*(1:size(A,1));

[X, Y] = meshgrid(xvec, yvec); %build screen coordinate matrices.
X0 = trapz(trapz(X.*A))/trapz(trapz(A)); %Mean value of x. A is 2D so need trapz(trapz()) for 2D integral. Assume symmetric distribution.
Y0 = trapz(trapz(Y.*A))/trapz(trapz(A)); %Mean value of y.
X = X-X0; %center coords.
Y = Y-Y0;
xvec = xvec-X0;
yvec = yvec-Y0;

[Xp,Yp]=meshgrid(1:size(A,2),1:size(A,1));
X0px=sum(A(:).*Xp(:))/sum(A(:));
Y0px=sum(A(:).*Yp(:))/sum(A(:));

try
    [S1,S2] = shear_fourier_v1(A);
catch
    [S1,S2] = shear_manual(A);
end

Asheared = shear_image(Xp,Yp,A,S1,S2,X0px,Y0px);

[locsx_sh,locsy_sh,minsx,minsy] = splitimage_TEM(Asheared,analysis);
[Locsx_sh,Locsy_sh] = meshgrid(locsx_sh,locsy_sh);
[Locsx,Locsy] = undo_shear(Locsx_sh,Locsy_sh,X0px,Y0px,S1,S2);
plot_screen_divided(A,Locsx,Locsy);
[xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,intx_unc,inty,inty_unc,sigmaxp,sigmayp,covx_scaled,covy_scaled] = phasespace_TEM_unc(A,X,Y,locsx_sh,locsy_sh,Locsx,Locsy,mask_prop,analysis);
[S,ex,ex_unc,ey,ey_unc,info] = beammatrix_TEM_unc(xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp,covx_scaled,covy_scaled,mask_prop);
[S_interp] = interpolation_TEM(Xp,Yp,info);

disp(length(locsx_sh))
disp(length(locsy_sh))

[ex,ey,e1,e2]=Emittance_2D_4D(S);
[ex_interp,ey_interp,e1_interp,e2_interp]=Emittance_2D_4D(S_interp);

    %% sort e1 & e2
sorted= sort([e1 e2]);
if ex<=ey
    e1=sorted(1); e2=sorted(2);
else
    e1=sorted(2); e2=sorted(1);
end
sorted= sort([e1_interp e2_interp]);
if ex_interp<=ey_interp
    e1_interp=sorted(1); e2_interp=sorted(2);
else
    e1_interp=sorted(2); e2_interp=sorted(1);
end   
