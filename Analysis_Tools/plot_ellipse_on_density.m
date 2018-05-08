function plot_ellipse_on_density(xc,yc,sx,sy,sxy)

phi=1/2*atand(2*sxy/(sx^2-sy^2));
xrms=sqrt((sx^2*cosd(phi)^2-sy^2*sind(phi)^2)/(cosd(phi)^4-sind(phi)^4));
yrms=sqrt((sx^2*sind(phi)^2-sy^2*cosd(phi)^2)/(sind(phi)^4-cosd(phi)^4));

t=0:20:360;
a1=xrms*cosd(t); a2=yrms*sind(t);
xellipse=xc+a1*cosd(phi)+a2*sind(phi);
yellipse=yc-a1*sind(phi)+a2*cosd(phi);
plot(xellipse,yellipse,'r-');
end
    
