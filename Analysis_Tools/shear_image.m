function [Ashear] = shear_image(X,Y,A,Sx,Sy,X0,Y0)

XX= X;
YY= Y; % use same grid settings as X and Y (pixels)

Xshear = X + Sx*(Y-Y0);
Yshear = Y + Sy*(X-X0);

F=scatteredInterpolant(Xshear(:),Yshear(:),A(:)); 
Ashear = F(XX,YY);