function yypphasespace(yg,yp,sigmayp,inty)


Y= linspace(min(yg(:)),max(yg(:)),size(yg,1));
Yp= linspace(min(yp(:)),max(yp(:)),200);

[Y,Yp]=meshgrid(Y,Yp); Int=zeros(size(Y));

gridsp= mean(mean(diff(yg,1,1)));

count=0;
for ii=1:size(Y,1)
    for jj=1:size(Y,2)
        
        Int(ii,jj) = ...
            sum( ...
            inty(:).* ...
            1/gridsp.*( heaviside(Y(ii,jj)-(yg(:)-gridsp/2))-heaviside(Y(ii,jj)-(yg(:)+gridsp/2)) ).* ...
            1/sqrt(2*pi)./sigmayp(:).*exp(-(Yp(ii,jj)-yp(:)).^2./(2*sigmayp(:).^2)) ...
            );
        count=count+1;
    end
%     display([num2str(count/numel(Y)*100) '%'])
%     imagesc(Y(1,:)/1e-3,Yp(:,1)/1e-3,Int); axis('xy'); xlabel('y [mm]'); ylabel('yp [mrad]'); drawnow
    
end

imagesc(Y(1,:)/1e-3,Yp(:,1)/1e-3,Int); axis('xy'); xlabel('y [mm]'); ylabel('yp [mrad]'); 