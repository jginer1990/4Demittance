function [ex, ey, e1, e2] = Emittance_2D_4D(Sig)

% Projected 2D emittance 
Sx = Sig(1:2,1:2);
ex = sqrt(det(Sx));

Sy = Sig(3:4,3:4);
ey = sqrt(det(Sy));

% Diagonalization 4D emittance

J = ...
    [0  1  0   0 ; ...
    -1   0   0   0 ; ...
    0   0   0   1 ; ...
    0   0   -1   0 ] ;
TJS = trace( (Sig*J)^2 );

e1 = 1/2 * sqrt( -TJS + sqrt(TJS^2-16*det(Sig)) )
e2 = 1/2 * sqrt( -TJS - sqrt(TJS^2-16*det(Sig)) )  

end