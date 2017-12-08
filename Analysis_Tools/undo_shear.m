function [Locsx,Locsy] = undo_shear(Locsx_sh,Locsy_sh,X0px,Y0px,S1,S2)
%SHEAR_IMAGE Inverse shear transformation

Locsx=X0px+1/(1-S1*S2)*((Locsx_sh-X0px) - S1*(Locsy_sh-Y0px));
Locsy=Y0px+1/(1-S1*S2)*(-S2*(Locsx_sh-X0px) +(Locsy_sh-Y0px));

end