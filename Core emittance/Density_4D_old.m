function [Int,X,Xp,Y,Yp,r0,Scov] = Density_4D(x,xp,y,yp,sxp,syp,sxpyp,int,lb,ub,Npoints)

if true
    sxp2=sxp; sxp2(sxp2==1)=0; syp2=syp; syp2(syp2==1)=0;
    figure(887); clf;
    subplot(241); imagesc(x); title('x')
    subplot(242); imagesc(xp); title('xp')
    subplot(243); imagesc(y); title('y')
    subplot(244); imagesc(yp); title('yp')
    subplot(245); imagesc(sxp2); title('sxp')
    subplot(246); imagesc(syp2); title('syp')
    subplot(247); imagesc(sxpyp); title('sxpyp')
    subplot(248); imagesc(int); title('int')
    drawnow
end


roi=find(int>0);
if nargin<=8 %default limits and number of points
    Npoints = 100*[1 1 1 1];
    lb= [min(x(roi)) min(xp(roi)) min(y(roi)) min(yp(roi))];
    ub= [max(x(roi)) max(xp(roi)) max(y(roi)) max(yp(roi))];
end

for i=1:length(Npoints) %limit number of points to 100
    if Npoints(i)>100
        Npoints(i)=100; display(['Number of points of variable ' num2str(i) 'is too large. Set to maximum 100.'])
    end
end

X_vec= linspace(lb(1),ub(1),Npoints(1)); dx= mean(mean(diff(x,1,2)));
Xp_vec= linspace(lb(2),ub(2),Npoints(2));
Y_vec= linspace(lb(3),ub(3),Npoints(3)); dy= mean(mean(diff(y,1,1)));
Yp_vec= linspace(lb(4),ub(4),Npoints(4));

[X,Xp,Y,Yp]=ndgrid(X_vec,Xp_vec,Y_vec,Yp_vec); 

x=x(roi);
xp=xp(roi);
y=y(roi);
yp=yp(roi);
sxp=sxp(roi);
syp=syp(roi);
sxpyp=sxpyp(roi);
int=int(roi);

Int=zeros(size(X)); 
fxxp=zeros([size(Int,1) size(Int,2) 1 1]);
fxy=zeros([size(Int,1) 1 size(Int,3) 1]);
fxyp=zeros([size(Int,1) 1 1 size(Int,4)]);
fxpy=zeros([1 size(Int,2) size(Int,3) 1]);
fxpyp=zeros([1 size(Int,2) 1 size(Int,4)]);
fyyp=zeros([1 1 size(Int,3) size(Int,4)]);

for i1=1:size(Int,1)
    for i2=1:size(Int,2);
        fxxp(i1,i2,1,1)=sum( ...
                    int(:).* ...
                    1/dx.*( heaviside(X(i1,i2,1,1)-(x(:)-dx/2))-heaviside(X(i1,i2,1,1)-(x(:)+dx/2)) ).* ...
                    1/sqrt(2*pi)./sxp(:).*exp(-(Xp(i1,i2,1,1)-xp(:)).^2./(2*sxp(:).^2)) ...
                    );
    end
end

for i1=1:size(Int,1)
    for i3=1:size(Int,3);
        fxy(i1,1,i3,1)=sum( ...
                    int(:).* ...
                    1/dx.*( heaviside(X(i1,1,i3,1)-(x(:)-dx/2))-heaviside(X(i1,1,i3,1)-(x(:)+dx/2)) ).* ...
                    1/dy.*( heaviside(Y(i1,1,i3,1)-(y(:)-dy/2))-heaviside(Y(i1,1,i3,1)-(y(:)+dy/2)) ) ...
                    );
    end
end

for i1=1:size(Int,1)
    for i4=1:size(Int,4);
        fxyp(i1,1,1,i4)=sum( ...
                    int(:).* ...
                    1/dx.*( heaviside(X(i1,1,1,i4)-(x(:)-dx/2))-heaviside(X(i1,1,1,i4)-(x(:)+dx/2)) ).* ...
                    1/sqrt(2*pi)./syp(:).*exp(-(Yp(i1,1,1,i4)-yp(:)).^2./(2*syp(:).^2)) ...
                    );
    end
end

for i2=1:size(Int,2)
    for i3=1:size(Int,3);
        fxpy(1,i2,i3,1)=sum( ...
                    int(:).* ...
                    1/dy.*( heaviside(Y(1,i2,i3,1)-(y(:)-dy/2))-heaviside(Y(1,i2,i3,1)-(y(:)+dy/2)) ).* ...
                    1/sqrt(2*pi)./sxp(:).*exp(-(Xp(1,i2,i3,1)-xp(:)).^2./(2*sxp(:).^2)) ...
                    );
    end
end

rho=sxpyp./sxp./syp;
for i2=1:size(Int,2)
    for i4=1:size(Int,4);
        fxpyp(1,i2,1,i4)=sum( ...
                    int(:).* ...
                    1./(2*pi.*sxp(:).*syp(:).*sqrt(1-rho(:).^2)).* ...
                    exp(-1./(2*(1-rho(:).^2)).*( ...
                        (Xp(1,i2,1,i4)-xp(:)).^2./(sxp(:).^2)+(Yp(1,i2,1,i4)-yp(:)).^2./(syp(:).^2) ...
                        -2*rho(:).*(Xp(1,i2,1,i4)-xp(:)).*(Yp(1,i2,1,i4)-yp(:))./(sxp(:).*syp(:)) ...
                        ) ) ...
                    );
    end
end

for i3=1:size(Int,3)
    for i4=1:size(Int,4);
        fyyp(1,1,i3,i4)=sum( ...
                    int(:).* ...
                    1/dy.*( heaviside(Y(1,1,i3,i4)-(y(:)-dy/2))-heaviside(Y(1,1,i3,i4)-(y(:)+dy/2)) ).* ...
                    1/sqrt(2*pi)./syp(:).*exp(-(Yp(1,1,i3,i4)-yp(:)).^2./(2*syp(:).^2)) ...
                    );
    end
end

Fxxp=squeeze(fxxp);
Fyyp=squeeze(fyyp);
Fxy=squeeze(fxy);
Fxyp=squeeze(fxyp);
Fxpy=squeeze(fxpy);
Fxpyp=squeeze(fxpyp);

fxxp= repmat(fxxp(:,:,1,1), [1 1 size(Int,3) size(Int,4)]);
fxy= repmat(fxy(:,1,:,1), [1 size(Int,2) 1 size(Int,4)]);
fxyp= repmat(fxyp(:,1,1,:), [1 size(Int,2) size(Int,3) 1]);
fxpy= repmat(fxpy(1,:,:,1), [size(Int,1) 1 1 size(Int,4)]);
fxpyp= repmat(fxpyp(1,:,1,:), [size(Int,1) 1 size(Int,3) 1]);
fyyp= repmat(fyyp(1,1,:,:), [size(Int,1) size(Int,2) 1 1]);

% Int = fxxp.*fxy.*fxyp.*fxpy.*fxpyp.*fyyp;
Int= fxxp.*fxy.*fxyp + fxxp.*fxpy.*fxpyp + fxy.*fxpy.*fyyp + fxyp.*fxpyp.*fyyp;

Int=Int/sum(Int(:));
r0=[sum(Int(:).*X(:)) , ...
    sum(Int(:).*Xp(:)) , ...
    sum(Int(:).*Y(:)) , ...
    sum(Int(:).*Yp(:)) , ...
    ];
Scov(1,1)=sum(Int(:).*(X(:)-r0(1)).^2) ;
Scov(2,2)=sum(Int(:).*(Xp(:)-r0(2)).^2) ;
Scov(3,3)=sum(Int(:).*(Y(:)-r0(3)).^2) ;
Scov(4,4)=sum(Int(:).*(Yp(:)-r0(4)).^2) ;
Scov(1,2)=sum(Int(:).*(X(:)-r0(1)).*(Xp(:)-r0(2))) ; Scov(2,1)=Scov(1,2);
Scov(1,3)=sum(Int(:).*(X(:)-r0(1)).*(Y(:)-r0(3))) ; Scov(3,1)=Scov(1,3);
Scov(1,4)=sum(Int(:).*(X(:)-r0(1)).*(Yp(:)-r0(4))) ; Scov(4,1)=Scov(1,4);
Scov(2,3)=sum(Int(:).*(Xp(:)-r0(2)).*(Y(:)-r0(3))) ; Scov(3,2)=Scov(2,3);
Scov(2,4)=sum(Int(:).*(Xp(:)-r0(2)).*(Yp(:)-r0(4))) ; Scov(4,2)=Scov(2,4);
Scov(3,4)=sum(Int(:).*(Y(:)-r0(3)).*(Yp(:)-r0(4))) ; Scov(4,3)=Scov(3,4);

%%
figure(888); clf
% Fxxp=sum(sum(Int,3),4);
% Fyyp=sum(sum(Int,1),2); Fyyp=squeeze(Fyyp);
% Fxy=sum(sum(Int,2),4); Fxy=squeeze(Fxy);
% Fxyp=sum(sum(Int,2),3); Fxyp=squeeze(Fxyp);
% Fxpy=sum(sum(Int,1),4); Fxpy=squeeze(Fxpy);
% Fxpyp=sum(sum(Int,1),3); Fxpyp=squeeze(Fxpyp);

subplot(321); imagesc(X_vec/1e-3,Xp_vec/1e-3,Fxxp); xlabel('x [mm]'); ylabel('xp [mrad]'); set(gca,'ydir','normal')
subplot(322); imagesc(Y_vec/1e-3,Yp_vec/1e-3,Fyyp); xlabel('y [mm]'); ylabel('yp [mrad]'); set(gca,'ydir','normal')
subplot(323); imagesc(X_vec/1e-3,Y_vec/1e-3,Fxy); xlabel('x [mm]'); ylabel('y [m]'); set(gca,'ydir','normal')
subplot(324); imagesc(X_vec/1e-3,Yp_vec/1e-3,Fxyp); xlabel('x [mm]'); ylabel('yp [mrad]'); set(gca,'ydir','normal')
subplot(325); imagesc(Xp_vec/1e-3,Y_vec/1e-3,Fxpy); xlabel('xp [mrad]'); ylabel('y [mm]'); set(gca,'ydir','normal')
subplot(326); imagesc(Xp_vec/1e-3,Yp_vec/1e-3,Fxpyp); xlabel('xp [mrad]'); ylabel('yp [mrad]'); set(gca,'ydir','normal')


dist={Fxxp,Fyyp,Fxy,Fxyp,Fxpy,Fxpyp};
arr={[X_vec;Xp_vec],[Y_vec;Yp_vec],[X_vec;Y_vec],[X_vec;Yp_vec],[Xp_vec;Y_vec],[Xp_vec;Yp_vec]};
for i=1:length(dist)
    [x1,x2]=meshgrid(arr{i}(1,:),arr{i}(2,:));
    mu(1)=sum(dist{i}(:).*x1(:))/sum(dist{i}(:));
    mu(2)=sum(dist{i}(:).*x2(:))/sum(dist{i}(:));
    
    co(1,1)=sum(dist{i}(:).*(x1(:)-mu(1)).^2)/sum(dist{i}(:));
    co(2,2)=sum(dist{i}(:).*(x2(:)-mu(2)).^2)/sum(dist{i}(:));
    co(1,2)=sum(dist{i}(:).*(x1(:)-mu(1)).*(x2(:)-mu(2)))/sum(dist{i}(:));
    co(2,1)=co(1,2);
    phi=1/2*atand(2*co(1,2)/(co(1,1)-co(2,2)));
    
    
    subplot(3,2,i); hold on
    t=0:20:360;
    x1sig=sqrt((co(1,1)*cosd(phi)^2-co(2,2)*sind(phi)^2)/(cosd(phi)^4-sind(phi)^4));
    x2sig=sqrt((co(1,1)*sind(phi)^2-co(2,2)*cosd(phi)^2)/(sind(phi)^4-cosd(phi)^4));
    a1=x1sig*cosd(t); a2=x2sig*sind(t);
    xellipse=mu(1)+a1*cosd(phi)-a2*sind(phi);
    yellipse=mu(2)+a1*sind(phi)+a2*cosd(phi);
    plot(xellipse/1e-3,yellipse/1e-3,'r-'); plot(mu(1)/1e-3,mu(2)/1e-3,'r+')
end
drawnow



end