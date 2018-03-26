function [S,ex,ey,info] = beammatrix_TEM_unc(xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp,covx_scaled,covy_scaled,mask_prop)
%BEAMMATRIX_TEM_UNC computes moments and beam matrix

%%  Compute moments

x0 = sum(sum((xg) .* intx))/sum(sum(intx));
y0 = sum(sum((yg) .* inty))/sum(sum(inty));
xp0 = sum(sum((xp) .* intx))/sum(sum(intx));
yp0 = sum(sum((yp) .* inty))/sum(sum(inty));

S(1,1) = sum(sum((xg-x0).^2 .* intx))/sum(sum(intx));
S(3,3) = sum(sum((yg-y0).^2 .* inty))/sum(sum(inty));

S(1,2) = sum(sum((xg-x0) .* (xp-xp0) .* intx))/sum(sum(intx));
S(3,4) = sum(sum((yg-y0) .* (yp-yp0) .* inty))/sum(sum(inty));

S(2,2) = sum(sum(((xp-xp0).^2+sigmaxp.^2) .* intx))/sum(sum(intx));
S(4,4) = sum(sum(((yp-yp0).^2+sigmayp.^2) .* inty))/sum(sum(inty));

%% Calculate uncertainties
for ii=1:length(xb(:))
    Vx_ind_cell{ii} = covx_scaled{ii};
    % Remove correlation terms for test purposes only
%     V_ind_cell{ii}(1,2) = 0;
%     V_ind_cell{ii}(2,1) = 0;
%     V_ind_cell{ii}(1,3) = 0;
%     V_ind_cell{ii}(3,1) = 0;
%     V_ind_cell{ii}(2,3) = 0;
%     V_ind_cell{ii}(3,2) = 0;
end

for ii=1:length(yb(:))
    Vy_ind_cell{ii} = covy_scaled{ii};
end

Vx_unc = blkdiag(Vx_ind_cell{:});    
Vy_unc = blkdiag(Vy_ind_cell{:});

Jx = zeros(3,3*length(xb(:)));
sum_intx = sum(intx(:));
for ii = 1:length(xb(:))
    Jx(1,1+3*(ii-1)) = 0;
    Jx(1,2+3*(ii-1)) = 0;
    Jx(1,3+3*(ii-1)) = (xg(ii)^2-S(1,1))/sum_intx;
    Jx(2,1+3*(ii-1)) = 2*intx(ii)*(xs(ii)-xg(ii))/(mask_prop.driftLength^2*sum_intx);
    Jx(2,2+3*(ii-1)) = 2*intx(ii)*sigmaxp(ii)/sum_intx;
    Jx(2,3+3*(ii-1)) = ((xs(ii)-xg(ii))^2/mask_prop.driftLength^2 + sigmaxp(ii)^2 - S(2,2))/sum_intx;
    Jx(3,1+3*(ii-1)) = intx(ii)*xg(ii)/(mask_prop.driftLength * sum_intx);
    Jx(3,2+3*(ii-1)) = 0;
    Jx(3,3+3*(ii-1)) = (xg(ii)*(xs(ii)-xg(ii))/mask_prop.driftLength - S(1,2))/sum_intx;
end

Jy = zeros(3,3*length(yb(:)));
sum_inty = sum(inty(:));
for ii = 1:length(yb(:))
    Jy(1,1+3*(ii-1)) = 0;
    Jy(1,2+3*(ii-1)) = 0;
    Jy(1,3+3*(ii-1)) = (yg(ii)^2-S(3,3))/sum_inty;
    Jy(2,1+3*(ii-1)) = 2*inty(ii)*(ys(ii)-yg(ii))/(mask_prop.driftLength^2*sum_inty);
    Jy(2,2+3*(ii-1)) = 2*inty(ii)*sigmayp(ii)/sum_inty;
    Jy(2,3+3*(ii-1)) = ((ys(ii)-yg(ii))^2/mask_prop.driftLength^2 + sigmayp(ii)^2 - S(4,4))/sum_inty;
    Jy(3,1+3*(ii-1)) = inty(ii)*yg(ii)/(mask_prop.driftLength * sum_inty);
    Jy(3,2+3*(ii-1)) = 0;
    Jy(3,3+3*(ii-1)) = (yg(ii)*(ys(ii)-yg(ii))/mask_prop.driftLength - S(3,4))/sum_inty;
end


covMatrix_elements_x = Jx*Vx_unc*Jx';
% Remove correlation terms for test purposes only
%     covMatrix_elements(1,2) = 0;
%     covMatrix_elements(2,1) = 0;
%     covMatrix_elements(1,3) = 0;
%     covMatrix_elements(3,1) = 0;
%     covMatrix_elements(2,3) = 0;
%     covMatrix_elements(3,2) = 0;

covMatrix_elements_y = Jy*Vy_unc*Jy';

% var_x2 = 0;
% for ii = 1:length(xb(:))
%     contribx2 = ((xg(ii)^2 - S(1,1))/sum_int)^2 * V_ind_cell{ii}(3,3);
%     var_x2 = var_x2 + contribx2;
% end
%  

ex = sqrt(S(1,1)*S(2,2)-S(1,2)^2);
ey = sqrt(S(3,3)*S(4,4)-S(3,4)^2);

% Method 1
J_ex_sq = [S(2,2) S(1,1) -2*S(1,2)];
J_ey_sq = [S(4,4) S(3,3) -2*S(3,4)];

ex_sq_unc = sqrt(J_ex_sq*covMatrix_elements_x*J_ex_sq');
ey_sq_unc = sqrt(J_ey_sq*covMatrix_elements_y*J_ey_sq');

ex_unc_m1 = ex_sq_unc/(2*ex);
ey_unc = ey_sq_unc/(2*ey);

% Method 2
J_ex = [0.5*S(2,2)/ex 0.5*S(1,1)/ex -S(1,2)/ex];
ex_unc_m2 = sqrt(J_ex*covMatrix_elements_x*J_ex');


%% Interpolate vertex values from midpoint values

xcen = zeros(size(xg,1)+1,size(xg,2));
xcen(1,:) = xg(1,:);
xcen(end,:) = xg(end,:);
xcen(2:end-1,:) = (xg(1:end-1,:)+xg(2:end,:))/2;

xpcen = zeros(size(xp,1)+1,size(xp,2));
xpcen(1,:) = xp(1,:);
xpcen(end,:) = xp(end,:);
xpcen(2:end-1,:) = (xp(1:end-1,:)+xp(2:end,:))/2;

xscen = zeros(size(xs,1)+1,size(xs,2));
xscen(1,:) = xs(1,:);
xscen(end,:) = xs(end,:);
xscen(2:end-1,:) = (xs(1:end-1,:)+xs(2:end,:))/2;

intxcen = zeros(size(intx,1)+1,size(intx,2));
intxcen(1,:) = intx(1,:);
intxcen(end,:) = intx(end,:);
intxcen(2:end-1,:) = (intx(1:end-1,:)+intx(2:end,:))/2.*sign(intx(1:end-1,:)+intx(2:end,:)); % if any is zero, make intcen=0

ycen = zeros(size(yg,1),size(yg,2)+1);
ycen(:,1) = yg(:,1);
ycen(:,end) = yg(:,end);
ycen(:,2:end-1) = (yg(:,1:end-1)+yg(:,2:end))/2;

ypcen = zeros(size(yp,1),size(yp,2)+1);
ypcen(:,1) = yp(:,1);
ypcen(:,end) = yp(:,end);
ypcen(:,2:end-1) = (yp(:,1:end-1)+yp(:,2:end))/2;

yscen = zeros(size(ys,1),size(ys,2)+1);
yscen(:,1) = ys(:,1);
yscen(:,end) = ys(:,end);
yscen(:,2:end-1) = (ys(:,1:end-1)+ys(:,2:end))/2;

intycen = zeros(size(inty,1),size(inty,2)+1);
intycen(:,1) = inty(:,1);
intycen(:,end) = inty(:,end);
intycen(:,2:end-1) = (inty(:,1:end-1)+inty(:,2:end))/2.*sign(inty(:,1:end-1)+inty(:,2:end)); % if any is zero, make intcen=0

intcen = sqrt(intxcen.*intycen);
figure(99); hold on; plot(xbcen(intcen~=0),ybcen(intcen~=0),'go'); 
figure(199); hold on
for ii=1:numel(intcen)
    plot3(xbcen(ii)*[1 1],ybcen(ii)*[1 1],[0 intcen(ii)],'g-');
end

x0c = sum(sum(xcen.*intcen))./sum(sum(intcen));
y0c = sum(sum(ycen.*intcen))./sum(sum(intcen));
xp0c = sum(sum(xpcen.*intcen))./sum(sum(intcen));
yp0c = sum(sum(ypcen.*intcen))./sum(sum(intcen));

% Calculate correlations
S(1,3) = sum(sum((xcen-x0c) .* (ycen-y0c) .* intcen))        ./sum(sum(intcen));
S(1,4) = sum(sum((xcen-x0c) .* (ypcen-yp0c) .* intcen))    ./sum(sum(intcen));
S(2,3) = sum(sum((xpcen-xp0c) .* (ycen-y0c) .* intcen))    ./sum(sum(intcen));
S(2,4) = sum(sum((xpcen-xp0c) .* (ypcen-yp0c) .* intcen))./sum(sum(intcen));

% Symmetrize the beam matrix
S(2,1) = S(1,2);
S(3,1) = S(1,3);
S(4,1) = S(1,4);
S(3,2) = S(2,3);
S(4,2) = S(2,4);
S(4,3) = S(3,4);


%% Calculate uncertainties for coupling terms

NoEle = 3;
% Caclulate uncertainties on xcens
NoCols = size(xs,2);
Jxcencol = zeros(NoEle*size(xscen,1),NoEle*size(xs,1));

for ii = 1:NoEle
    Jxcencol(ii,ii) = 1;
    Jxcencol(size(Jxcencol,1)+1-ii,size(Jxcencol,2)+1-ii) = 1;
end

for ii = NoEle+1:size(Jxcencol,2)
    Jxcencol(ii,ii) = 0.5;
    Jxcencol(ii,ii-3) = 0.5;
end

Jxcencell = repmat({Jxcencol},1,NoCols);
Jxcen = blkdiag(Jxcencell{:}); % Make matrix for every column

Vx_unc_cen = Jxcen*Vx_unc*Jxcen';


% Caclulate uncertainties on ycens
NoRows = size(ys,1);
Jycenrow = zeros(NoEle*size(yscen,2),NoEle*size(ys,2));

for ii = 1:NoEle
    Jycenrow(ii,ii) = 1;
    Jycenrow(size(Jycenrow,1)+1-ii,size(Jycenrow,2)+1-ii) = 1;
end

for ii = NoEle+1:size(Jycenrow,2)
    Jycenrow(ii,ii) = 0.5;
    Jycenrow(ii,ii-3) = 0.5;
end

Jycencell = repmat({Jycenrow},1,NoRows);
Jycen = blkdiag(Jycencell{:}); % Make matrix for every row

Vy_unc_cen = Jycen*Vy_unc*Jycen';


% Calculate uncertainties on beam matrix elements

Vx_unc_cen_nosigma = Vx_unc_cen;
Vx_unc_cen_nosigma(2:3:end,:)=[]; % Remove sigma rows and columns from Vx_unc_cen
Vx_unc_cen_nosigma(:,2:3:end)=[]; % Remove sigma rows and columns from Vx_unc_cen
Vy_unc_cen_nosigma = Vy_unc_cen;
Vy_unc_cen_nosigma(2:3:end,:)=[]; % Remove sigma rows and columns from Vy_unc_cen
Vy_unc_cen_nosigma(:,2:3:end)=[]; % Remove sigma rows and columns from Vy_unc_cen

Vx_unc_cen_long = [Vx_unc_cen_nosigma zeros(size(Vx_unc_cen_nosigma))]; % Add zeros to end of row
Vy_unc_cen_long = [zeros(size(Vy_unc_cen_nosigma)) Vy_unc_cen_nosigma]; % Add zeros to end of row
V_unc_cen = vertcat(Vx_unc_cen_long,Vy_unc_cen_long);

% Reorder so we have x1,y1,Ix1,Iy1,x2,y2...
order = zeros(1,length(V_unc_cen));
order(1:2:end) = 1:length(V_unc_cen)/2;
order(2:2:end) = length(V_unc_cen)/2+1:length(V_unc_cen);
V_unc_cen_re = V_unc_cen(order,order);


Jxy = zeros(4,4*length(xcen(:))); % Derivative matrix for <xy>,<x'y'>,<xy'>,<x'y> with respect to xcen(i),Ix(i),ycen(i),Iy(i)

sum_intcen = sum(intcen(:));
dIdIx = 0.5*sqrt(intycen./intxcen);
dIdIy = 0.5*sqrt(intxcen./intycen);
sigmaxpyp = zeros(size(dIdIx));

for ii = 1:length(xcen(:))
    Jxy(1,1+4*(ii-1)) = 0;
    Jxy(1,2+4*(ii-1)) = 0;
    Jxy(1,3+4*(ii-1)) = (xcen(ii)*ycen(ii)-S(1,3))/sum_intcen * dIdIx(ii);
    Jxy(1,4+4*(ii-1)) = (xcen(ii)*ycen(ii)-S(1,3))/sum_intcen * dIdIy(ii);
    Jxy(2,1+4*(ii-1)) = intcen(ii)*(yscen(ii)-ycen(ii))/(mask_prop.driftLength^2*sum_intcen);
    Jxy(2,2+4*(ii-1)) = intcen(ii)*(xscen(ii)-xcen(ii))/(mask_prop.driftLength^2*sum_intcen);
    Jxy(2,3+4*(ii-1)) = ((xscen(ii)-xcen(ii))*(yscen(ii)-ycen(ii))/mask_prop.driftLength^2 + sigmaxpyp(ii)^2 - S(2,4))/sum_intcen * dIdIx(ii);
    Jxy(2,4+4*(ii-1)) = ((xscen(ii)-xcen(ii))*(yscen(ii)-ycen(ii))/mask_prop.driftLength^2 + sigmaxpyp(ii)^2 - S(2,4))/sum_intcen * dIdIy(ii);
    Jxy(3,1+4*(ii-1)) = 0;
    Jxy(3,2+4*(ii-1)) = intcen(ii)*xcen(ii)/(mask_prop.driftLength*sum_intcen);
    Jxy(3,3+4*(ii-1)) = (xcen(ii)*(yscen(ii)-ycen(ii))/mask_prop.driftLength - S(1,4))/sum_intcen * dIdIx(ii);
    Jxy(3,4+4*(ii-1)) = (xcen(ii)*(yscen(ii)-ycen(ii))/mask_prop.driftLength - S(1,4))/sum_intcen * dIdIy(ii);
    Jxy(4,1+4*(ii-1)) = intcen(ii)*ycen(ii)/(mask_prop.driftLength*sum_intcen);
    Jxy(4,2+4*(ii-1)) = 0;
    Jxy(4,3+4*(ii-1)) = (ycen(ii)*(xscen(ii)-xcen(ii))/mask_prop.driftLength - S(2,3))/sum_intcen * dIdIx(ii);
    Jxy(4,4+4*(ii-1)) = (ycen(ii)*(xscen(ii)-xcen(ii))/mask_prop.driftLength - S(2,3))/sum_intcen * dIdIy(ii);
end


covMatrix_elements_xy = Jxy*V_unc_cen_re*Jxy';


%% Uncertainties on intrinsic emittances

covMatrix_elements_combined = blkdiag(covMatrix_elements_x,covMatrix_elements_y,covMatrix_elements_xy); % order: <x^2>,<x'^2>,<xx'>,<y^2>,<y'^2>,<yy'>,<xy>,<x'y'>,<xy'>,<x'y> 

xx = S(1,1);
xpxp = S(2,2);
xxp = S(1,2);
yy = S(3,3);
ypyp = S(4,4);
yyp = S(3,4);
xy = S(1,3);
xpyp = S(2,4);
xyp = S(1,4);
xpy = S(2,3);

conveniencefun = sqrt((2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)^2 ...
-16*(xyp^2*xpy^2 - 2*xy*xyp*xpy*xpyp + xy^2*xpyp^2 - xyp^2*xpxp*yy + 2*xxp*xyp*xpyp*yy - xx*xpyp^2*yy ... 
+ 2*xy*xyp*xpxp*yyp - 2*xxp*xyp*xpy*yyp - 2*xxp*xy*xpyp*yyp + 2*xx*xpy*xpyp*yyp + xxp^2*yyp^2 - xx*xpxp*yyp^2 ...
- xy^2*xpxp*ypyp + 2*xxp*xy*xpy*ypyp - xx*xpy^2*ypyp - xxp^2*yy*ypyp + xx*xpxp*yy*ypyp));

denominator_de1 = 4*sqrt(-2*xxp^2 + 2*xx*xpxp - 4*xyp*xpy + 4*xy*xpyp - 2*yyp^2 + 2*yy*ypyp + conveniencefun);
denominator_de2 = 4*sqrt(-2*xxp^2 + 2*xx*xpxp - 4*xyp*xpy + 4*xy*xpyp - 2*yyp^2 + 2*yy*ypyp - conveniencefun);

de1_dxx = (2*xpxp + (-4*xpxp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp) ...
-16*(-xpyp^2*yy + 2*xpy*xpyp*yyp - xpxp*yyp^2 - xpy^2*ypyp + xpxp*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dxx = (2*xpxp - (-4*xpxp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp) ...
-16*(-xpyp^2*yy + 2*xpy*xpyp*yyp - xpxp*yyp^2 - xpy^2*ypyp + xpxp*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dxpxp = (2*xx + (-4*xx*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp) ...
-16*(-xyp^2*yy + 2*xy*xyp*yyp - xx*yyp^2 - xy^2*ypyp + xx*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dxpxp = (2*xx - (-4*xx*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp) ...
-16*(-xyp^2*yy + 2*xy*xyp*yyp - xx*yyp^2 - xy^2*ypyp + xx*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dxxp = (-4*xxp + (8*xxp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp) ...
- 16*(2*xyp*xpyp*yy - 2*xyp*xpy*yyp - 2*xy*xpyp*yyp + 2*xxp*yyp^2 + 2*xy*xpy*ypyp - 2*xxp*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dxxp = (-4*xxp - (8*xxp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp) ...
- 16*(2*xyp*xpyp*yy - 2*xyp*xpy*yyp - 2*xy*xpyp*yyp + 2*xxp*yyp^2 + 2*xy*xpy*ypyp - 2*xxp*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dyy = (2*ypyp + (-16*(-xyp^2*xpxp + 2*xxp*xyp*xpyp - xx*xpyp^2 - xxp^2*ypyp + xx*xpxp*ypyp) ...
- 4*ypyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dyy = (2*ypyp - (-16*(-xyp^2*xpxp + 2*xxp*xyp*xpyp - xx*xpyp^2 - xxp^2*ypyp + xx*xpxp*ypyp) ...
- 4*ypyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dypyp = (2*yy + (-16*(-xy^2*xpxp + 2*xxp*xy*xpy - xx*xpy^2 - xxp^2*yy + xx*xpxp*yy) ...
- 4*yy*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dypyp = (2*yy - (-16*(-xy^2*xpxp + 2*xxp*xy*xpy - xx*xpy^2 - xxp^2*yy + xx*xpxp*yy) ...
- 4*yy*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dyyp = (-4*yyp + (-16*(2*xy*xyp*xpxp - 2*xxp*xyp*xpy - 2*xxp*xy*xpyp + 2*xx*xpy*xpyp + 2*xxp^2*yyp - 2*xx*xpxp*yyp) ...
+ 8*yyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dyyp = (-4*yyp - (-16*(2*xy*xyp*xpxp - 2*xxp*xyp*xpy - 2*xxp*xy*xpyp + 2*xx*xpy*xpyp + 2*xxp^2*yyp - 2*xx*xpxp*yyp) ...
+ 8*yyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dxy = (4*xpyp + (-16*(-2*xyp*xpy*xpyp + 2*xy*xpyp^2 + 2*xyp*xpxp*yyp - 2*xxp*xpyp*yyp - 2*xy*xpxp*ypyp + 2*xxp*xpy*ypyp) ...
- 8*xpyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dxy = (4*xpyp - (-16*(-2*xyp*xpy*xpyp + 2*xy*xpyp^2 + 2*xyp*xpxp*yyp - 2*xxp*xpyp*yyp - 2*xy*xpxp*ypyp + 2*xxp*xpy*ypyp) ...
- 8*xpyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dxpyp = (4*xy + (-16*(-2*xy*xyp*xpy + 2*xy^2*xpyp + 2*xxp*xyp*yy - 2*xx*xpyp*yy - 2*xxp*xy*yyp + 2*xx*xpy*yyp) ...
- 8*xy*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dxpyp = (4*xy - (-16*(-2*xy*xyp*xpy + 2*xy^2*xpyp + 2*xxp*xyp*yy - 2*xx*xpyp*yy - 2*xxp*xy*yyp + 2*xx*xpy*yyp) ...
- 8*xy*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dxyp = (-4*xpy + (-16*(2*xyp*xpy^2 - 2*xy*xpy*xpyp - 2*xyp*xpxp*yy + 2*xxp*xpyp*yy + 2*xy*xpxp*yyp - 2*xxp*xpy*yyp) ...
+ 8*xpy*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de1;

de2_dxyp = (-4*xpy - (-16*(2*xyp*xpy^2 - 2*xy*xpy*xpyp - 2*xyp*xpxp*yy + 2*xxp*xpyp*yy + 2*xy*xpxp*yyp - 2*xxp*xpy*yyp) ...
+ 8*xpy*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp)) / (2*conveniencefun)) ...
/ denominator_de2;

de1_dxpy = (-4*xyp + (-16*(2*xyp^2*xpy - 2*xy*xyp*xpyp - 2*xxp*xyp*yyp + 2*xx*xpyp*yyp + 2*xxp*xy*ypyp - 2*xx*xpy*ypyp) ...
+ 8*xyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp))/ (2*conveniencefun)) ...
/ denominator_de1;

de2_dxpy = (-4*xyp - (-16*(2*xyp^2*xpy - 2*xy*xyp*xpyp - 2*xxp*xyp*yyp + 2*xx*xpyp*yyp + 2*xxp*xy*ypyp - 2*xx*xpy*ypyp) ...
+ 8*xyp*(2*xxp^2 - 2*xx*xpxp + 4*xyp*xpy - 4*xy*xpyp + 2*yyp^2 - 2*yy*ypyp))/ (2*conveniencefun)) ...
/ denominator_de2;

Jemittance = [de1_dxx de1_dxpxp de1_dxxp de1_dyy de1_dypyp de1_dyyp de1_dxy de1_dxpyp de1_dxyp de1_dxpy; ...
    de2_dxx de2_dxpxp de2_dxxp de2_dyy de2_dypyp de2_dyyp de2_dxy de2_dxpyp de2_dxyp de2_dxpy];

covMatrix_e1_e2 = Jemittance*covMatrix_elements_combined*Jemittance';

[ex,ey,e1,e2]=Emittance_2D_4D(S);

J_e4D = [e2 e1];

unc_e4D = sqrt(J_e4D*covMatrix_e1_e2*J_e4D');

%% Struct with all information of analysis
info.target='TEM';
info.xb=xb;
info.ybc=ybc;
info.yb=yb;
info.xbc=xbc;
info.xbcen=xbcen;
info.ybcen=ybcen;
info.sigmaxp=sigmaxp;
info.sigmayp=sigmayp;
info.sigmaxpyp=0;
info.xg=xg;
info.yg=yg;
info.xs=xs;
info.ys=ys;
info.xp=xp;
info.yp=yp;
info.x0=x0;
info.y0=y0;
info.x0c=x0c;
info.y0c=y0c;
info.xpcen=xpcen;
info.ypcen=ypcen;
info.xp0=xp0;
info.yp0=yp0;
info.xp0c=xp0c;
info.yp0c=yp0c;
info.intx=intx;
info.inty=inty;
info.xcen=xcen;
info.ycen=ycen;
info.intxcen=intxcen;
info.intycen=intycen;
info.intcen=intcen;
info.S=S;

display('*Success')

end