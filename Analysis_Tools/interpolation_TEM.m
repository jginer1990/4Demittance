function [S_interp] = phasespace_shear_interpolation(X,Y,info)
%PHASESPACE_SHEAR_INTERPOLATION Find beam matrix by interpolation

xb = info.xb;
yb = info.yb;
ybc = info.ybc;
xg = info.xg;
yg = info.yg;
x0 = info.x0;
y0 = info.y0;
x0c = info.x0c;
y0c = info.y0c;
xp0 = info.xp0;
yp0 = info.yp0;
xp0c = info.xp0c;
yp0c  = info.yp0c;
xp = info.xp;
yp = info.yp;
xbc = info.xbc;
ybc = info.ybc;
xcen = info.xcen;
ycen = info.ycen;
xpcen = info.xpcen;
ypcen = info.ypcen;
sigmaxp = info.sigmaxp;
sigmayp = info.sigmayp;
intx = info.intx;
inty = info.inty;
intcen = info.intcen;

S_interp =zeros(4);

X_interp = min(X(1,:)):1:max(X(1,:));
Y_interp = min(Y(:,1)):1:max(Y(:,1));
[X_interp,Y_interp]=meshgrid(X_interp,Y_interp);

index = find(xb(:)~=0 & ybc(:)~=0 & intx(:)~=0);
F= scatteredInterpolant(xb(index),ybc(index),intx(index),'linear','none'); intx_interp = F(X_interp,Y_interp); intx_interp(isnan(intx_interp))=0;
F= scatteredInterpolant(xb(index),ybc(index),xg(index)-x0); xg_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xb(index),ybc(index),xp(index)-xp0,'linear','none'); xp_interp = F(X_interp,Y_interp); xp_interp(isnan(xp_interp))=0;
F= scatteredInterpolant(xb(index),ybc(index),sigmaxp(index),'linear','none'); sigmaxp_interp = F(X_interp,Y_interp); sigmaxp_interp(isnan(sigmaxp_interp))=1;
index = find(xbc(:)~=0 & yb(:)~=0 & inty(:)~=0);
F= scatteredInterpolant(xbc(index),yb(index),inty(index),'linear','none'); inty_interp = F(X_interp,Y_interp); inty_interp(isnan(inty_interp))=0;
F= scatteredInterpolant(xbc(index),yb(index),yg(index)-y0); yg_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xbc(index),yb(index),yp(index)-yp0,'linear','none'); yp_interp = F(X_interp,Y_interp); yp_interp(isnan(yp_interp))=0;
F= scatteredInterpolant(xbc(index),yb(index),sigmayp(index),'linear','none'); sigmayp_interp = F(X_interp,Y_interp); sigmayp_interp(isnan(sigmayp_interp))=1;
sigmaxpyp_interp=zeros(size(sigmaxp_interp));

x0_interp = sum(sum((xg_interp) .* intx_interp))/sum(sum(intx_interp));
y0_interp = sum(sum((yg_interp) .* inty_interp))/sum(sum(inty_interp));
xp0_interp = sum(sum((xp_interp) .* intx_interp))/sum(sum(intx_interp));
yp0_interp = sum(sum((yp_interp).* inty_interp))/sum(sum(inty_interp));

S_interp(1,1) = sum(sum((xg_interp-x0_interp).^2 .* intx_interp))/sum(sum(intx_interp));
S_interp(3,3) = sum(sum((yg_interp-y0_interp).^2 .* inty_interp))/sum(sum(inty_interp));
S_interp(1,2) = sum(sum((xg_interp-x0_interp) .* (xp_interp-xp0_interp) .* intx_interp))/sum(sum(intx_interp));
S_interp(3,4) = sum(sum((yg_interp-y0_interp) .* (yp_interp-yp0_interp) .* inty_interp))/sum(sum(inty_interp));
S_interp(2,2) = sum(sum(((xp_interp-xp0_interp).^2+sigmaxp_interp.^2) .* intx_interp))/sum(sum(intx_interp));
S_interp(4,4) = sum(sum(((yp_interp-yp0_interp).^2+sigmayp_interp.^2) .* inty_interp))/sum(sum(inty_interp));

% intxy_interp=sqrt(intx_interp.*inty_interp); [mx1,my1]= find(intxy_interp>0,1,'first'); [mx2,my2]= find(intxy_interp>0,1,'last');
% [Int_4D,X_4D,Xp_4D,Y_4D,Yp_4D,r0_4D,S_4D]=...
%     Density_4D_v3(xg_interp(mx1:mx2,my1:my2),xp_interp(mx1:mx2,my1:my2),yg_interp(mx1:mx2,my1:my2),yp_interp(mx1:mx2,my1:my2),...
%     sigmaxp_interp(mx1:mx2,my1:my2),sigmayp_interp(mx1:mx2,my1:my2),sigmaxpyp_interp(mx1:mx2,my1:my2),intxy_interp(mx1:mx2,my1:my2));




%%
clear ybmid ybcen xbmid xbcen
ybmid(:,2:size(yb,2)+1) = yb(:,:);
ybmid(:,1) = ybmid(:,2);
ybmid(:,end+1) = ybmid(:,end);
ybcen = sign(ybmid(:,1:end-1)).^2.*sign(ybmid(:,2:end)).^2.*(ybmid(:,1:end-1)+ybmid(:,2:end)) /2;

xbmid(:,2:size(xb,1)+1) = xb(:,:)';
xbmid(:,1) = xbmid(:,2);
xbmid(:,end+1) = xbmid(:,end);
xbcen = (sign(xbmid(:,1:end-1)).^2.*sign(xbmid(:,2:end)).^2.*(xbmid(:,1:end-1)+xbmid(:,2:end))/2)';

index = find(xbcen(:)~=0 & ybcen(:)~=0 & intcen(:)~=0);
F= scatteredInterpolant(xbcen(index),ybcen(index),intcen(index),'linear','none'); intcen_interp = F(X_interp,Y_interp); intcen_interp(isnan(intcen_interp))=0;
F= scatteredInterpolant(xbcen(index),ybcen(index),xcen(index)-x0c); xcen_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xbcen(index),ybcen(index),xpcen(index)-xp0c,'linear','none'); xpcen_interp = F(X_interp,Y_interp); xpcen_interp(isnan(xpcen_interp))=0;
F= scatteredInterpolant(xbcen(index),ybcen(index),ycen(index)-y0c); ycen_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xbcen(index),ybcen(index),ypcen(index)-yp0c,'linear','none'); ypcen_interp = F(X_interp,Y_interp); ypcen_interp(isnan(ypcen_interp))=0;


x0c_interp = sum(sum(xcen_interp.*intcen_interp))./sum(sum(intcen_interp));
y0c_interp = sum(sum(ycen_interp.*intcen_interp))./sum(sum(intcen_interp));
xp0c_interp = sum(sum(xpcen_interp.*intcen_interp))./sum(sum(intcen_interp));
yp0c_interp = sum(sum(ypcen_interp.*intcen_interp))./sum(sum(intcen_interp));

% Calculate correlations
S_interp(1,3) = sum(sum((xcen_interp-x0c_interp) .* (ycen_interp-y0c_interp) .* intcen_interp)) ./sum(sum(intcen_interp));
S_interp(1,4) = sum(sum((xcen_interp-x0c_interp) .* (ypcen_interp-yp0c_interp) .* intcen_interp)) ./sum(sum(intcen_interp));
S_interp(2,3) = sum(sum((xpcen_interp-xp0c_interp) .* (ycen_interp-y0c_interp) .* intcen_interp)) ./sum(sum(intcen_interp));
S_interp(2,4) = sum(sum((xpcen_interp-xp0c_interp) .* (ypcen_interp-yp0c_interp) .* intcen_interp)) ./sum(sum(intcen_interp));

% Symmetrize the beam matrix
S_interp(2,1) = S_interp(1,2);
S_interp(3,1) = S_interp(1,3);
S_interp(4,1) = S_interp(1,4);
S_interp(3,2) = S_interp(2,3);
S_interp(4,2) = S_interp(2,4);
S_interp(4,3) = S_interp(3,4);
