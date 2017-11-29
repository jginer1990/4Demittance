function xyphasespace(xg,yg,intxy)


X= linspace(min(xg(:)),max(xg(:)),size(xg,2));
Y= linspace(min(yg(:)),max(yg(:)),size(yg,1));

[X,Y]=meshgrid(X,Y); Int=zeros(size(X));

gridspx= mean(mean(diff(xg,1,2)));
gridspy= mean(mean(diff(yg,1,1)));

count=0;
for ii=1:size(X,1)
    for jj=1:size(X,2)
        
        Int(ii,jj) = ...
            sum( ...
            intxy(:).* ...
            1/gridspx.*( heaviside(X(ii,jj)-(xg(:)-gridspx/2))-heaviside(X(ii,jj)-(xg(:)+gridspx/2)) ).* ...
            1/gridspy.*( heaviside(Y(ii,jj)-(yg(:)-gridspy/2))-heaviside(Y(ii,jj)-(yg(:)+gridspy/2)) ) ...
            );
        count=count+1;
    end
%     display([num2str(count/numel(X)*100) '%'])
%     imagesc(X(1,:)/1e-3,Y(:,1)/1e-3,Int); axis('xy'); xlabel('x [mm]'); ylabel('xp [mrad]');  drawnow

end

imagesc(X(1,:)/1e-3,Y(:,1)/1e-3,Int); axis('xy'); xlabel('x [mm]'); ylabel('y [mm]'); 