function [Xv,Xpv,Int] = xxpphasespace(xg,xp,sigmaxp,intx)


X= linspace(min(xg(:)),max(xg(:)),size(xg,2));
Xp= linspace(min(xp(:)),max(xp(:)),200);

[X,Xp]=meshgrid(X,Xp); Int=zeros(size(X));

gridsp= mean(mean(diff(xg,1,2)));

count=0;
for ii=1:size(X,1)
    for jj=1:size(X,2)
        
        Int(ii,jj) = ...
            sum( ...
            intx(:).* ...
            1/gridsp.*( heaviside(X(ii,jj)-(xg(:)-gridsp/2))-heaviside(X(ii,jj)-(xg(:)+gridsp/2)) ).* ...
            1/sqrt(2*pi)./sigmaxp(:).*exp(-(Xp(ii,jj)-xp(:)).^2./(2*sigmaxp(:).^2)) ...
            );
        count=count+1;
    end
%     display([num2str(count/numel(X)*100) '%'])
%     imagesc(X(1,:)/1e-3,Xp(:,1)/1e-3,Int); axis('xy'); xlabel('x [mm]'); ylabel('xp [mrad]');  drawnow

end
Xv=X(1,:); Xpv=Xp(:,1)';
imagesc(X(1,:)/1e-3,Xp(:,1)/1e-3,Int); axis('xy'); xlabel('x [mm]'); ylabel('xp [mrad]'); 