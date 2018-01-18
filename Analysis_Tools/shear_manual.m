function [S1,S2] = shear_manual(rhon)
%SHEAR_MANUAL Manually select points to find shear factors

figure(3)
imagesc(rhon)
set(gca,'YDir','normal')

disp('Trace VERTICAL bars')
[x1,y1]=ginput(1);
hold on;
plot(x1,y1,'rx')
[x2,y2]=ginput(1);
plot(x2,y2,'rx')
plot([x1 x2],[y1 y2],'r-')
disp('Trace HORIZONTAL bars')
[x3,y3]=ginput(1);
hold on;
plot(x3,y3,'mx')
[x4,y4]=ginput(1);
plot(x4,y4,'mx')
plot([x3 x4],[y3 y4],'m-')
drawnow

if abs(atand((y1-y2)/(x1-x2)))>abs(atand((y4-y3)/(x4-x3)))
    S1 = -(x1-x2)/(y1-y2);
    S2 = -(y4-y3)/(x4-x3);
else
    S2 = -(y1-y2)/(x1-x2);
    S1 = -(x4-x3)/(y4-y3);
end

end
