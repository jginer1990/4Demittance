function [Ashear] = shear_image(X,Y,A,Sx,Sy,X0,Y0)
%SHEAR_IMAGE Shears image A with coordinates X,Y by shear factors Sx and Sy
%around X0 and Y0

Xshear = X + Sx*(Y-Y0);
Yshear = Y + Sy*(X-X0);

F=scatteredInterpolant(Xshear(:),Yshear(:),A(:)); 
Ashear = F(X,Y);

end