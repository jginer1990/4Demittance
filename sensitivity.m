close all
clear

addpath './Analysis_Tools'
%load('origMatrix_64_24_cor_FBT1_PSF.mat')

origMatrix = [9.05e-8,2.21e-7,3.79e-10,4.40e-10;
   2.21e-7,5.49e-7,5.64e-10,6.17e-10;
    3.79e-10,5.64e-10,1.24e-8,1.59e-8;
    4.40e-10,6.17e-10,1.59e-8,2.04e-8];

xx = origMatrix(1,1);
xxp = origMatrix(1,2);
xpxp = origMatrix(2,2);
yy = origMatrix(3,3);
yyp = origMatrix(3,4);
ypyp = origMatrix(4,4);
xy = origMatrix(1,3);
xpy = origMatrix(2,3);
xyp = origMatrix(1,4);
xpyp = origMatrix(2,4);

varlist = [xx xxp xpxp yy yyp ypyp xy xpy xyp xpyp];

for ii = 1:length(varlist)
    var = varlist(ii);
    varvalues = linspace(var*0.9,var*1.1,200);
    
    for jj = 1:length(varvalues)
        v = varvalues(jj);
        
        newlist = varlist;
        newlist(ii) = v;
        
        M(1,1) = newlist(1);
        M(1,2) = newlist(2);
        M(2,2) = newlist(3);
        M(3,3) = newlist(4);
        M(3,4) = newlist(5);
        M(4,4) = newlist(6);
        M(1,3) = newlist(7);
        M(2,3) = newlist(8);
        M(1,4) = newlist(9);
        M(2,4) = newlist(10);

        M(2,1) = M(1,2);
        M(3,1) = M(1,3);
        M(4,1) = M(1,4);
        M(3,2) = M(2,3);
        M(4,2) = M(2,4);
        M(4,3) = M(3,4);
        
        [ex,ey,e1,e2] = Emittance_2D_4D(M);
        
        exlist(ii,jj) = ex;
        eylist(ii,jj) = ey;
        e1list(ii,jj) = e1;
        e2list(ii,jj) = e2;
        e4Dlist(ii,jj) = e1*e2;
        
        exlist(imag(exlist)~=0) = nan;
        eylist(imag(eylist)~=0) = nan;
        e1list(imag(e1list)~=0) = nan;
        e2list(imag(e2list)~=0) = nan;
        e4Dlist(imag(e4Dlist)~=0) = nan;
        
    end
end

%% Plot
xlist = linspace(-0.1,0.1,200);

figure
plot(xlist,exlist(1:3,:),'--')
legend('<x^2>','<xx''>','<x''x''>')
xlabel('Fractional error on matrix element')
ylabel('\epsilon_x [m rad]')

figure
plot(xlist,eylist(4:6,:),'--')
legend('<y^2>','<yy''>','<y''y''>')
xlabel('Fractional error on matrix element')
ylabel('\epsilon_y [m rad]')

figure
plot(xlist,e1list(1:10,:),'--')
legend('<x^2>','<xx''>','<x''x''>','<y^2>','<yy''>','<y''y''>','<xy>','<x''y>','<xy''>','<x''y''>')
xlabel('Fractional error on matrix element')
ylabel('\epsilon_1 [m rad]')

figure
plot(xlist,e2list(1:10,:),'--')
legend('<x^2>','<xx''>','<x''x''>','<y^2>','<yy''>','<y''y''>','<xy>','<x''y>','<xy''>','<x''y''>')
xlabel('Fractional error on matrix element')
ylabel('\epsilon_2 [m rad]')

figure
plot(xlist,e4Dlist(1:10,:),'--')
legend('<x^2>','<xx''>','<x''x''>','<y^2>','<yy''>','<y''y''>','<xy>','<x''y>','<xy''>','<x''y''>')
xlabel('Fractional error on matrix element')
ylabel('\epsilon_{4D} [(m rad)^2]')
