function [A_cropped,Acrop0] = crop(A,center,rx,ry)
%CROP crop image inside the elipse defined by center, rx and ry
 
    t = linspace(0,2*pi,100);
    xel=center(1)+rx*cos(t);
    yel=center(2)+ry*sin(t);
    for i= 1:length(xel)
        BGel(i)=A(round(yel(i)),round(xel(i))); % Background region just inside ellipse
    end

    hold on; plot(xel,yel,'y-')
    Xcrop = 1:size(A,2); Ycrop=1:size(A,1);
    [Xcrop,Ycrop]=meshgrid(Xcrop,Ycrop);
    Xcrop=Xcrop(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx);
    Ycrop=Ycrop(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx);
    Acrop= A(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx); Acrop2=Acrop;
    for i=1:size(Acrop,1)
        for j=1:size(Acrop,2)

            if (Xcrop(i,j)-center(1))^2/rx^2+(Ycrop(i,j)-center(2))^2/ry^2>1
                Acrop(i,j)=mean(BGel);
                Acrop0(i,j)=0;
            end
        end
    end

    A_cropped=Acrop-mean(BGel);
    A_cropped(A_cropped<0)=0; 
    figure(3); clf; imagesc(A_cropped); set(gca,'YDir','normal')

end

