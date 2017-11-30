function [Int,X,Xp,r0,Scov] = Density_2D(x,xp,sxp,int,Npoints,lb,ub)

if false
    sxp2=sxp; sxp2(sxp2==1)=0;
    figure(887); clf;
    subplot(141); imagesc(x); title('x')
    subplot(142); imagesc(xp); title('xp')
    subplot(143); imagesc(sxp2); title('sx')
    subplot(144); imagesc(int); title('int')
    drawnow
end


roi=find(int>0);
if nargin<=4 %default limits and number of points
    Npoints = 250*[1 1];
    lb= [min(x(roi)) min(xp(roi))];
    ub= [max(x(roi)) max(xp(roi))];
elseif nargin<=5
    lb= [min(x(roi)) min(xp(roi))];
    ub= [max(x(roi)) max(xp(roi))];
end

for i=1:length(Npoints) %limit number of points to 100
    if Npoints(i)>250
        Npoints(i)=250; display(['Number of points of variable ' num2str(i) ' is too large. Set to maximum 250.'])
    end
end

X_vec= linspace(lb(1),ub(1),Npoints(1)); dX=X_vec(2)-X_vec(1);
Xp_vec= linspace(lb(2),ub(2),Npoints(2));

[X,Xp]=meshgrid(X_vec,Xp_vec); Int=zeros(size(X)); 

x=x(roi);
xp=xp(roi);
sxp=sxp(roi);
int=int(roi);

for i=1:numel(x)
    d=abs(X-x(i)); minval=min(d(:));
    aux1=zeros(size(d)); aux1(d==minval)=1/dX; %to do: check what happens when minval is too large
    aux2=1/(sqrt(2*pi)*sxp(i))* exp(-1/2*(Xp-xp(i)).^2/sxp(i)^2);
                    
    Int = Int + int(i).*aux1.*aux2;
%     display([num2str(i/numel(x)*100) '%'])
end

Int=Int/sum(Int(:));
r0=[sum(Int(:).*X(:)) , ...
    sum(Int(:).*Xp(:))];

Scov(1,1)=sum(Int(:).*(X(:)-r0(1)).^2) ;
Scov(2,2)=sum(Int(:).*(Xp(:)-r0(2)).^2) ;
Scov(1,2)=sum(Int(:).*(X(:)-r0(1)).*(Xp(:)-r0(2))) ; Scov(2,1)=Scov(1,2);

% figure(1086); clf
imagesc(X_vec/1e-3,Xp_vec/1e-3,Int'); xlabel('x [mm]'); ylabel('xp [mrad]'); set(gca,'ydir','normal')

phi=1/2*atand(2*Scov(1,2)/(Scov(1,1)-Scov(2,2)));

hold on
t=0:20:360;
x1sig=sqrt((Scov(1,1)*cosd(phi)^2-Scov(2,2)*sind(phi)^2)/(cosd(phi)^4-sind(phi)^4));
x2sig=sqrt((Scov(1,1)*sind(phi)^2-Scov(2,2)*cosd(phi)^2)/(sind(phi)^4-cosd(phi)^4));
a1=x1sig*cosd(t); a2=x2sig*sind(t);
xellipse=r0(1)+a1*cosd(phi)-a2*sind(phi);
yellipse=r0(2)+a1*sind(phi)+a2*cosd(phi);
plot(xellipse/1e-3,yellipse/1e-3,'r-'); plot(r0(1)/1e-3,r0(2)/1e-3,'r+')