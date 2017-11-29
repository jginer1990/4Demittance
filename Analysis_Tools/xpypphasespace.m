function xpypphasespace(xp,yp,sigmaxp,sigmayp,int)


Xp= linspace(min(xp(:)),max(xp(:)),200);
Yp= linspace(min(yp(:)),max(yp(:)),200);

[Xp,Yp]=meshgrid(Xp,Yp); Int=zeros(size(Xp));

count=0;
for ii=1:size(Xp,1)
    for jj=1:size(Xp,2)
        
        Int(ii,jj) = ...
            sum( ...
            int(:).* ...
            1/sqrt(2*pi)./sigmaxp(:).*exp(-(Xp(ii,jj)-xp(:)).^2./(2*sigmaxp(:).^2)).* ...
            1/sqrt(2*pi)./sigmayp(:).*exp(-(Yp(ii,jj)-yp(:)).^2./(2*sigmayp(:).^2)) ...
            ); %To do: Needs to take into account sigmaxpyp!!
        count=count+1;
    end
%     display([num2str(count/numel(Xp)*100) '%'])
end

imagesc(Xp(1,:)/1e-3,Yp(:,1)/1e-3,Int); axis('xy'); xlabel('xp [mrad]'); ylabel('yp [mrad]'); 