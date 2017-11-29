X_interp = min(X(1,:)):4:max(X(1,:));
Y_interp = min(Y(:,1)):4:max(Y(:,1));
[X_interp,Y_interp]=meshgrid(X_interp,Y_interp);

index = find(xb(:)~=0 & ybc(:)~=0);
F= scatteredInterpolant(xb(index),ybc(index),intx(index),'linear','none'); intx_interp = F(X_interp,Y_interp); intx_interp(isnan(intx_interp))=0;
% F= scatteredInterpolant(xb(index),ybc(index),xg(index)-x0); xg_interp = F(X_interp,Y_interp);
% F= scatteredInterpolant(xb(index),ybc(index),xp(index)-xp0,'linear','none'); xp_interp = F(X_interp,Y_interp); xp_interp(isnan(xp_interp))=0;
% F= scatteredInterpolant(xb(index),ybc(index),sigmaxp(index),'linear','none'); sigmaxp_interp = F(X_interp,Y_interp); sigmaxp_interp(isnan(sigmaxp_interp))=1;
index = find(xbc(:)~=0 & yb(:)~=0);
F= scatteredInterpolant(xbc(index),yb(index),inty(index),'linear','none'); inty_interp = F(X_interp,Y_interp); inty_interp(isnan(inty_interp))=0;
% F= scatteredInterpolant(xbc(index),yb(index),yg(index)-y0); yg_interp = F(X_interp,Y_interp);
% F= scatteredInterpolant(xbc(index),yb(index),yp(index)-yp0,'linear','none'); yp_interp = F(X_interp,Y_interp); yp_interp(isnan(yp_interp))=0;
% F= scatteredInterpolant(xbc(index),yb(index),sigmayp(index),'linear','none'); sigmayp_interp = F(X_interp,Y_interp); sigmayp_interp(isnan(sigmayp_interp))=1;

X_interp = min(min(X_interp(:,find(sum(intx_interp,1)~=0)))):1:max(max(X_interp(:,find(sum(intx_interp,1)~=0))));
Y_interp = min(min(Y_interp(find(sum(inty_interp,2)~=0),:))):1:max(max(Y_interp(find(sum(inty_interp,2)~=0),:)));
[X_interp,Y_interp]=meshgrid(X_interp,Y_interp);
index = find(xb(:)~=0 & ybc(:)~=0);
F= scatteredInterpolant(xb(index),ybc(index),intx(index),'linear','none'); intx_interp = F(X_interp,Y_interp); intx_interp(isnan(intx_interp)| intx_interp<0)=0;
F= scatteredInterpolant(xb(index),ybc(index),xg(index)-x0); xg_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xb(index),ybc(index),xp(index)-xp0,'linear','none'); xp_interp = F(X_interp,Y_interp); xp_interp(isnan(xp_interp))=0;
F= scatteredInterpolant(xb(index),ybc(index),sigmaxp(index),'linear','none'); sigmaxp_interp = F(X_interp,Y_interp); sigmaxp_interp(isnan(sigmaxp_interp))=1;
index = find(xbc(:)~=0 & yb(:)~=0);
F= scatteredInterpolant(xbc(index),yb(index),inty(index),'linear','none'); inty_interp = F(X_interp,Y_interp); inty_interp(isnan(inty_interp) | inty_interp<0)=0;
F= scatteredInterpolant(xbc(index),yb(index),yg(index)-y0); yg_interp = F(X_interp,Y_interp);
F= scatteredInterpolant(xbc(index),yb(index),yp(index)-yp0,'linear','none'); yp_interp = F(X_interp,Y_interp); yp_interp(isnan(yp_interp))=0;
F= scatteredInterpolant(xbc(index),yb(index),sigmayp(index),'linear','none'); sigmayp_interp = F(X_interp,Y_interp); sigmayp_interp(isnan(sigmayp_interp))=1;

[Int_4D,X_4D,Xp_4D,Y_4D,Yp_4D,r0_4D,S_4D]=...
    Density_4D(xg_interp,xp_interp,yg_interp,yp_interp,...
    sigmaxp_interp,sigmayp_interp,sigmaxpyp_interp,sqrt(intx_interp.*inty_interp));

%%
figure(81); clf; display('* x - xp phase space')
xxpphasespace(xg_interp,xp_interp,sigmaxp_interp,intx_interp); 
cmap=colormap;
cmap(1,:)=[1 1 1];
% cmap(2,:)=[1 1 1];
colormap(cmap);
drawnow


figure(87); clf; display('* x - xp phase space uncorrelated with contours')
corr = trapz(trapz(intx_interp.*xg_interp.*xp_interp))/trapz(trapz(intx_interp.*xg_interp.^2));  xp_uncorr= xp_interp-xg_interp.*corr;
[x1,xp1,Ixxp]=xxpphasespace(xg_interp,xp_uncorr,sigmaxp_interp,intx_interp); 
cmap=colormap;
cmap(1,:)=[1 1 1];
% cmap(2,:)=[1 1 1];
colormap(cmap);
hold on
contour(x1/1e-3,xp1/1e-3,Ixxp,max(Ixxp(:))*[0.3 0.6],'m--','linewidth',2)
drawnow
%%
figure(82); clf; display('* y - yp phase space')
% corr = trapz(trapz(inty_interp.*yg_interp.*yp_interp))/trapz(trapz(inty_interp.*yg_interp.^2));  yp_interp= yp_interp-yg_interp.*corr;
yypphasespace(yg_interp,yp_interp,sigmayp_interp,inty_interp); drawnow
cmap=colormap;
cmap(1,:)=[1 1 1];
% cmap(2,:)=[1 1 1];
colormap(cmap);
drawnow

%%

figure(83); clf; display('* x - yp phase space')
xxpphasespace(xg_interp,yp_interp,sigmayp_interp,sqrt(intx_interp.*inty_interp));
xlabel('x [mm]'); ylabel('yp [mrad]'); drawnow

figure(84); clf; display('* y - xp phase space')
yypphasespace(yg_interp,xp_interp,sigmaxp_interp,sqrt(intx_interp.*inty_interp));
xlabel('y [mm]'); ylabel('xp [mrad]'); drawnow

figure(85); clf; display('* xp - yp phase space')
xpypphasespace(xp_interp,yp_interp,sigmaxp_interp,sigmayp_interp,sqrt(intx_interp.*inty_interp)); drawnow

figure(86); clf; display('* x - y phase space')
xyphasespace(xg_interp,yg_interp,sqrt(intx_interp.*inty_interp)); drawnow