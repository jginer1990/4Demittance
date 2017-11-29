function [Int,X,Xp,Y,Yp,r0,Scov] = Density_4D_interpolation(x,xp,y,yp,sxp,syp,sxpyp,int,lb,ub,Npoints)

rho=sxpyp./sxp./syp;

nb_sigmas = 3;
npoints_mesh=100;

Int_list=[];
R_list=[];
for i=1:numel(x)
    
    xp_mesh= linspace(xp(i)-nb_sigmas*sxp(i),xp(i)+nb_sigmas*sxp(i),npoints_mesh);
    yp_mesh= linspace(yp(i)-nb_sigmas*syp(i),yp(i)+nb_sigmas*syp(i),npoints_mesh);
    [xp_mesh,yp_mesh]=meshgrid(xp_mesh,yp_mesh);
    x_mesh=x(i)*ones(size(xp_mesh));
    y_mesh=y(i)*ones(size(xp_mesh));
    int_mesh= 1./(2*pi.*sxp(i).*syp(i).*sqrt(1-rho(i).^2)).* ...
                    exp(-1./(2*(1-rho(i).^2)).*( ...
                        (xp_mesh-xp(i)).^2./(sxp(i).^2)+(yp_mesh-yp(i)).^2./(syp(i).^2) ...
                        -2*rho(i).*(xp_mesh-xp(i)).*(yp_mesh-yp(i))./(sxp(i).*syp(i)) ...
                        ) );
    
    R_list=[R_list; [x_mesh(:) xp_mesh(:) y_mesh(:) yp_mesh(:)]];
    Int_list=[Int_list; int_mesh];
end

%%

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

X_vec= linspace(lb(1),ub(1),Npoints(1)); 
Xp_vec= linspace(lb(2),ub(2),Npoints(2));
Y_vec= linspace(lb(3),ub(3),Npoints(3)); 
Yp_vec= linspace(lb(4),ub(4),Npoints(4));
[X,Xp,Y,Yp]=ndgrid(X_vec,Xp_vec,Y_vec,Yp_vec); 


%%
figure(886); clf
Fxxp=sum(sum(Int,3),4);
Fyyp=sum(sum(Int,1),2); Fyyp=squeeze(Fyyp);
Fxy=sum(sum(Int,2),4); Fxy=squeeze(Fxy);
Fxyp=sum(sum(Int,2),3); Fxyp=squeeze(Fxyp);
Fxpy=sum(sum(Int,1),4); Fxpy=squeeze(Fxpy);
Fxpyp=sum(sum(Int,1),3); Fxpyp=squeeze(Fxpyp);

subplot(321); imagesc(X_vec/1e-3,Xp_vec/1e-3,Fxxp); xlabel('x [mm]'); ylabel('xp [mrad]'); set(gca,'ydir','normal')
subplot(322); imagesc(Y_vec/1e-3,Yp_vec/1e-3,Fyyp); xlabel('y [mm]'); ylabel('yp [mrad]'); set(gca,'ydir','normal')
subplot(323); imagesc(X_vec/1e-3,Y_vec/1e-3,Fxy); xlabel('x [mm]'); ylabel('y [m]'); set(gca,'ydir','normal')
subplot(324); imagesc(X_vec/1e-3,Yp_vec/1e-3,Fxyp); xlabel('x [mm]'); ylabel('yp [mrad]'); set(gca,'ydir','normal')
subplot(325); imagesc(Xp_vec/1e-3,Y_vec/1e-3,Fxpy); xlabel('xp [mrad]'); ylabel('y [mm]'); set(gca,'ydir','normal')
subplot(326); imagesc(Xp_vec/1e-3,Yp_vec/1e-3,Fxpyp); xlabel('xp [mrad]'); ylabel('yp [mrad]'); set(gca,'ydir','normal')


% dist={Fxxp,Fyyp,Fxy,Fxyp,Fxpy,Fxpyp};
% arr={[X_vec;Xp_vec],[Y_vec;Yp_vec],[X_vec;Y_vec],[X_vec;Yp_vec],[Xp_vec;Y_vec],[Xp_vec;Yp_vec]};
% for i=1:length(dist)
%     [x1,x2]=meshgrid(arr{i}(1,:),arr{i}(2,:));
%     mu(1)=sum(dist{i}(:).*x1(:))/sum(dist{i}(:));
%     mu(2)=sum(dist{i}(:).*x2(:))/sum(dist{i}(:));
%     
%     co(1,1)=sum(dist{i}(:).*(x1(:)-mu(1)).^2)/sum(dist{i}(:));
%     co(2,2)=sum(dist{i}(:).*(x2(:)-mu(2)).^2)/sum(dist{i}(:));
%     co(1,2)=sum(dist{i}(:).*(x1(:)-mu(1)).*(x2(:)-mu(2)))/sum(dist{i}(:));
%     co(2,1)=co(1,2);
%     phi=1/2*atand(2*co(1,2)/(co(1,1)-co(2,2)));
%     
%     
%     subplot(3,2,i); hold on
%     t=0:20:360;
%     x1sig=sqrt((co(1,1)*cosd(phi)^2-co(2,2)*sind(phi)^2)/(cosd(phi)^4-sind(phi)^4));
%     x2sig=sqrt((co(1,1)*sind(phi)^2-co(2,2)*cosd(phi)^2)/(sind(phi)^4-cosd(phi)^4));
%     a1=x1sig*cosd(t); a2=x2sig*sind(t);
%     xellipse=mu(1)+a1*cosd(phi)-a2*sind(phi);
%     yellipse=mu(2)+a1*sind(phi)+a2*cosd(phi);
%     plot(xellipse/1e-3,yellipse/1e-3,'r-'); plot(mu(1)/1e-3,mu(2)/1e-3,'r+')
% end
% drawnow

dist={Fxxp,Fyyp,Fxy,Fxyp,Fxpy,Fxpyp};
arr={r0([1 2]),r0([3 4]),r0([1 3]),r0([1 4]),r0([2 3]),r0([2 4])};
matr={Scov([1 2],[1 2]),Scov([3 4],[3 4]),Scov([1 3],[1 3]),Scov([1 4],[1 4]),Scov([2 3],[2 3]),Scov([2 4],[2 4])};
for i=1:length(dist)
    mu=arr{i};
    
    co=matr{i};
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