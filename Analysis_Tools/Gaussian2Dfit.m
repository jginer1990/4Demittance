function [centroid,xsig,ysig,Sx,Sy,Sxy,phi_new,int] = Gaussian2Dfit(image, mask, ellipse)
%   @image the beam shot being analyzed
%   @mask region of intrest given as a 0 or 1 matrix (0 for false)
%   @ellipse=[xcenter,ycenter,xrms,yrms,sxy] an estimate rms ellipse, the rms dim should be close for the
%   fit to work
%   @ centroid: new centroid calculated from the gaussian fit of lineouts
%   @ xsig,ysig: fitted beam sizes at rotated angle without correlation 
%   @ Sx,Sy,Sxy: fitted beam sizes and correlation translated to the
%   original axes
%   @ phi_new: correlation angle of the 2Ddistribution corrected with fit

imgflag = false; % Plot fit

[row col]=find(mask);
rowmin=min(row);
rowmax=max(row);
colmin=min(col);
colmax=max(col);

%minimum number of points
if rowmin==0||rowmax==0||colmin==0||colmax==0
    disp('!Error fit 2D gaussian');
    quit;
end
if rowmax-rowmin <10
    rowmax=rowmax+15; rowmin=rowmin-15;
end
if colmax-colmin <10
    colmax=colmin+15; colmin=colmin-15;
end


image=image(rowmin:rowmax,colmin:colmax);
X=1:size(image,2);
Y=1:size(image,1);
CenterIm=[(colmin+colmax)/2,(rowmin+rowmax)/2];
CenterRot=[(min(X)+max(X))/2,(min(Y)+max(Y))/2];
c_x=ellipse(1)-CenterIm(1); c_y=ellipse(2)-CenterIm(2);
sx=ellipse(3); sy=ellipse(4); sxy=ellipse(5); phi=1/2*atand(2*sxy/(sx^2-sy^2));
center_x = CenterRot(1)+cosd(phi)*c_x-sind(phi)*c_y; 
center_y = CenterRot(2)+sind(phi)*c_x+cosd(phi)*c_y;
xrms=sqrt((sx^2*cosd(phi)^2-sy^2*sind(phi)^2)/(cosd(phi)^4-sind(phi)^4));
yrms=sqrt((sx^2*sind(phi)^2-sy^2*cosd(phi)^2)/(sind(phi)^4-cosd(phi)^4));

%Get lineouts of rotated image
rotimg=imrotate(image,-phi,'bilinear','crop');
lineoutx=sum(rotimg);
lineouty=sum(rotimg,2)';
heightx=max(lineoutx);
heighty=max(lineouty);

if imgflag
    figure(889); clf
    subplot(1,3,1); imagesc(X+colmin-1,Y+rowmin-1,image); title('ROI original'); set(gca,'ydir','normal')
    subplot(1,3,2); imagesc(X,Y,rotimg); title('ROI rotated'); set(gca,'ydir','normal')
    
    t=0:20:360;
    a1=xrms*cosd(t); a2=yrms*sind(t);
    xellipse=ellipse(1)+a1*cosd(phi)+a2*sind(phi);
    yellipse=ellipse(2)-a1*sind(phi)+a2*cosd(phi);
    subplot(1,3,1); hold on; plot(xellipse,yellipse,'g-'); plot(ellipse(1),ellipse(2),'g+')
    plot(CenterIm(1),CenterIm(2),'kx') %center image around which it rotates
    
    a1=center_x+xrms*cosd(t); a2=center_y+yrms*sind(t);
    subplot(1,3,2); hold on; plot(a1,a2,'g-'); plot(center_x,center_y,'g+')
    plot(CenterRot(1),CenterRot(2),'kx') %center image around which it rotates
    drawnow
end

% Fit model (gaussian)
F=@(x,xdata)(x(3)+x(4).*exp(-(1/2).*((((xdata-x(2))/x(1)).^2))));
x0=[xrms,center_x,0,heightx];
y0=[yrms,center_y,0,heighty];
lb=[0,0,-1,0];
ub=100*x0+[0,0,(heightx+heighty),0];
if isequal(lb,ub)
    ub=ub+[1,1,1,1];
end
options = optimset('Display','off');

ub=[length(X),length(X),heightx+1,2*heightx+1];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,X,lineoutx,lb,ub,options);
ub=[length(Y),length(Y),heighty+1,2*heighty+1];
[y,resnorm,~,exitflag,output] = lsqcurvefit(F,y0,Y,lineouty,lb,ub,options);

if imgflag
    figure(889); subplot(1,3,3)
    scatter(X,lineoutx,'MarkerEdgeColor','b');
    hold on;
    scatter(Y,lineouty,'MarkerEdgeColor','g');
    plot(X,F(x,X),'Color','b');
    plot(Y,F(y,Y),'Color','g');
    legend('xdata','ydata')
    drawnow
end


% centroid at original image (with xy correlation)
centroid=[(x(2)-CenterRot(1))*cosd(-phi) - (y(2)-CenterRot(2))*sind(-phi), ...
    (x(2)-CenterRot(1))*sind(-phi) + (y(2)-CenterRot(2))*cosd(-phi)] ...
    + CenterIm; % centroid [x,y] ([column,row])

% sigmas at rotated image (no correlation)
xsig=x(1);
ysig=y(1);

% sigmas at original image (with xy correlation)
Sx=sqrt(xsig^2*cosd(-phi)^2+ysig^2*sind(-phi)^2);
Sy=sqrt(ysig^2*cosd(-phi)^2+xsig^2*sind(-phi)^2);
Sxy=cosd(-phi)*sind(-phi)*(ysig^2-xsig^2);
phi_new=1/2*atand(2*Sxy/(Sx^2-Sy^2)); % Same as phi (cross-check)

if imgflag
    figure(889); 
  
    t=0:20:360;
    a1=xsig*cosd(t); a2=ysig*sind(t);
    xellipse=centroid(1)+a1*cosd(phi_new)+a2*sind(phi_new);
    yellipse=centroid(2)-a1*sind(phi_new)+a2*cosd(phi_new);
    subplot(131); hold on; plot(xellipse,yellipse,'r-'); plot(centroid(1),centroid(2),'r+')
    
    a1=x(2)+xsig*cosd(t); a2=y(2)+ysig*sind(t);
    subplot(132); hold on; plot(a1,a2,'r-'); plot(x(2),y(2),'r+')
    drawnow
end

int=mean([x(4)*sqrt(2*pi)*x(1)  y(4)*sqrt(2*pi)*y(1)]);

% Avoid large sigmas from bad fit.
if Sx>sx | Sy>sy
    centroid=ellipse(1:2);
    Sx= ellipse(3);
    Sy= ellipse(4);
    Sxy= ellipse(5);
    int= NaN;
end

end