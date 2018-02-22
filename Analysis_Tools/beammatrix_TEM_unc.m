function [S,ex,ex_unc,ey,ey_unc,info] = beammatrix_TEM_unc(xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp,covx_scaled,covy_scaled,mask_prop)
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
J_ex_sq = [S(2,2) S(1,1) -2*S(1,2)];
J_ey_sq = [S(4,4) S(3,3) -2*S(3,4)];

ex_sq_var = J_ex_sq*covMatrix_elements_x*J_ex_sq';
ey_sq_var = J_ey_sq*covMatrix_elements_y*J_ey_sq';

ex = sqrt(S(1,1)*S(2,2)-S(1,2)^2);
ey = sqrt(S(3,3)*S(4,4)-S(3,4)^2);
ex_unc = sqrt(ex_sq_var)/(2*ex);
ey_unc = sqrt(ey_sq_var)/(2*ey);


%% Interpolate vertex values from midpoint values

clear ymid xmid intymid intxmid intxcen intycen ycen xcen
clear ypmid xpmid ypcen xpcen
ymid(:,2:size(yb,2)+1) = yg(:,:);
ymid(:,1) = ymid(:,2);
ymid(:,end+1) = ymid(:,end);
ycen = (ymid(:,1:end-1)+ymid(:,2:end)) /2;

ypmid(:,2:size(yb,2)+1) = yp(:,:);
ypmid(:,1) = ypmid(:,2);
ypmid(:,end+1) = ypmid(:,end);
ypcen = (ypmid(:,1:end-1)+ypmid(:,2:end)) /2;

intymid(:,2:size(yb,2)+1) = inty(:,:);
intymid(:,1) = intymid(:,2);
intymid(:,end+1) = intymid(:,end);
intycen = (intymid(:,1:end-1)+intymid(:,2:end)) /2 .*sign(intymid(:,1:end-1).*intymid(:,2:end)); % if any is zero, make intcen=0

xmid(:,2:size(xb,1)+1) = xg(:,:)';
xmid(:,1) = xmid(:,2);
xmid(:,end+1) = xmid(:,end);
xcen = ((xmid(:,1:end-1)+xmid(:,2:end))/2)';

xpmid(:,2:size(xb,1)+1) = xp(:,:)';
xpmid(:,1) = xpmid(:,2);
xpmid(:,end+1) = xpmid(:,end);
xpcen = ((xpmid(:,1:end-1)+xpmid(:,2:end))/2)';

intxmid(:,2:size(xb,1)+1) = intx(:,:)';
intxmid(:,1) = intxmid(:,2);
intxmid(:,end+1) = intxmid(:,end);
intxcen = ((intxmid(:,1:end-1)+intxmid(:,2:end))/2.*sign(intxmid(:,1:end-1).*intxmid(:,2:end)))'; % if any is zero, make intcen=0

intcen = sqrt(intycen.*intxcen);
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


% Struct with all information of analysis
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