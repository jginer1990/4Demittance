S_interp =zeros(4);

X_interp = min(X(1,:)):1:max(X(1,:));
Y_interp = min(Y(:,1)):1:max(Y(:,1));
[X_interp,Y_interp]=meshgrid(X_interp,Y_interp);

index = find(xb(:)~=0 & yb(:)~=0 & int(:)~=0);
F= scatteredInterpolant(xb(index),yb(index),int(index),'linear','none'); int_interp = F(X_interp,Y_interp); int_interp(isnan(int_interp))=0;
F= scatteredInterpolant(xb(index),yb(index),xg(index)-x0); xg_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xb(index),yb(index),xp(index)-xp0,'linear','none'); xp_interp = F(X_interp,Y_interp); xp_interp(isnan(xp_interp))=0;
F= scatteredInterpolant(xb(index),yb(index),sigmaxp(index),'linear','none'); sigmaxp_interp = F(X_interp,Y_interp); sigmaxp_interp(isnan(sigmaxp_interp))=1;
F= scatteredInterpolant(xb(index),yb(index),yg(index)-y0); yg_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xb(index),yb(index),yp(index)-yp0,'linear','none'); yp_interp = F(X_interp,Y_interp); yp_interp(isnan(yp_interp))=0;
F= scatteredInterpolant(xb(index),yb(index),sigmayp(index),'linear','none'); sigmayp_interp = F(X_interp,Y_interp); sigmayp_interp(isnan(sigmayp_interp))=1;
F= scatteredInterpolant(xb(index),yb(index),sigmaxpyp(index),'linear','none'); sigmaxpyp_interp = F(X_interp,Y_interp); sigmaxpyp_interp(isnan(sigmaxpyp_interp))=1;

x0_interp = sum(sum((xg_interp) .* int_interp))/sum(sum(int_interp));
y0_interp = sum(sum((yg_interp) .* int_interp))/sum(sum(int_interp));
xp0_interp = sum(sum((xp_interp) .* int_interp))/sum(sum(int_interp));
yp0_interp = sum(sum((yp_interp).* int_interp))/sum(sum(int_interp));

% Compute moments
S_interp(1,1) = sum(sum((xg_interp-x0_interp).^2 .* int_interp))/sum(sum(int_interp));
S_interp(3,3) = sum(sum((yg_interp-y0_interp).^2 .* int_interp))/sum(sum(int_interp));
S_interp(1,2) = sum(sum((xg_interp-x0_interp) .* (xp_interp-xp0_interp) .* int_interp))/sum(sum(int_interp));
S_interp(3,4) = sum(sum((yg_interp-y0_interp) .* (yp_interp-yp0_interp) .* int_interp))/sum(sum(int_interp));
S_interp(2,2) = sum(sum(((xp_interp-xp0_interp).^2+sigmaxp_interp.^2) .* int_interp))/sum(sum(int_interp));
S_interp(4,4) = sum(sum(((yp_interp-yp0_interp).^2+sigmayp_interp.^2) .* int_interp))/sum(sum(int_interp));

% Calculate correlations
S_interp(1,3) = sum(sum((xg_interp-x0_interp) .* (yg_interp-y0_interp) .* int_interp)) ./sum(sum(int_interp));
S_interp(1,4) = sum(sum((xg_interp-x0_interp) .* (yp_interp-yp0_interp) .* int_interp)) ./sum(sum(int_interp));
S_interp(2,3) = sum(sum((xp_interp-xp0_interp) .* (yg_interp-y0_interp) .* int_interp)) ./sum(sum(int_interp));
S_interp(2,4) = sum(sum((xp_interp-xp0_interp) .* (yp_interp-yp0_interp) .* int_interp)) ./sum(sum(int_interp));

% Symmetrize the beam matrix
S_interp(2,1) = S_interp(1,2);
S_interp(3,1) = S_interp(1,3);
S_interp(4,1) = S_interp(1,4);
S_interp(3,2) = S_interp(2,3);
S_interp(4,2) = S_interp(2,4);
S_interp(4,3) = S_interp(3,4);