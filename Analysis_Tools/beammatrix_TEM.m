function [S,info] = beammatrix_TEM(xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp)
%BEAMMATRIX_TEM computes moments and beam matrix

%%  Compute moments

x0 = sum(sum((xg) .* intx))/sum(sum(intx));
y0 = sum(sum((yg) .* inty))/sum(sum(inty));
xp0 = sum(sum((xp) .* intx))/sum(sum(intx));
yp0 = sum(sum((yp).* inty))/sum(sum(inty));

S(1,1) = sum(sum((xg-x0).^2 .* intx))/sum(sum(intx));
S(3,3) = sum(sum((yg-y0).^2 .* inty))/sum(sum(inty));

S(1,2) = sum(sum((xg-x0) .* (xp-xp0) .* intx))/sum(sum(intx));
S(3,4) = sum(sum((yg-y0) .* (yp-yp0) .* inty))/sum(sum(inty));

S(2,2) = sum(sum(((xp-xp0).^2+sigmaxp.^2) .* intx))/sum(sum(intx));
S(4,4) = sum(sum(((yp-yp0).^2+sigmayp.^2) .* inty))/sum(sum(inty));

% Interpolate vertex values from midpoint values

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