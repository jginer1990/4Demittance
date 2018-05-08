function [S,info] = phasespace_PP(rhon,X,Y,Xn,Yn,X0px,Y0px,S1,S2,mask_prop,locsx,locsy)
%PHASESPACE_PP analyse image and compute parameters


%% Determine beamlets centroids and rms ellipse parameters


xb = zeros(length(locsy)-1,length(locsx)-1);
yb=xb; sx=xb; sy=xb; sxy=xb; int=xb;
for i=1:length(locsy)-1
    
    for j=1:length(locsx)-1
        
        Xroi_sh = locsx(j):locsx(j+1);
        Yroi_sh = locsy(i):locsy(i+1);
        [Xroi_sh,Yroi_sh]=meshgrid(Xroi_sh,Yroi_sh);
        Xroi = round( X0px+ 1/(1-S1*S2)*( (Xroi_sh-X0px)-S1*(Yroi_sh-Y0px) ) );
        Yroi = round( Y0px+ 1/(1-S1*S2)*( -S2*(Xroi_sh-X0px)+(Yroi_sh-Y0px) ) );
        Iroi =zeros(size(Xroi));
        mask=zeros(size(rhon));
        for ii=1:size(Xroi,1)
            for jj=1:size(Xroi,2)
                if Xroi(ii,jj)>0 & Xroi(ii,jj)<=size(rhon,2) & Yroi(ii,jj)>0 & Yroi(ii,jj)<=size(rhon,1)
                    Iroi(ii,jj) = rhon(Yroi(ii,jj),Xroi(ii,jj));
                    mask(Yroi(ii,jj),Xroi(ii,jj)) = 1;
                end
            end
        end
        
        int(i,j) = sum(Iroi(:));
        if int(i,j)>0
            xb(i,j) = sum(Iroi(:).*Xroi(:))/int(i,j);
            sx(i,j) = sqrt( sum(Iroi(:).*(Xroi(:)-xb(i,j)).^2)/int(i,j) );
            yb(i,j) = sum(Iroi(:).*Yroi(:))/int(i,j);
            sy(i,j) = sqrt( sum(Iroi(:).*(Yroi(:)-yb(i,j)).^2)/int(i,j) );
            sxy(i,j) = sum(Iroi(:).*(Xroi(:)-xb(i,j)).*(Yroi(:)-yb(i,j)))/int(i,j);
        end
%         figure(888); clf; h1= surf(Xroi,Yroi,Iroi); set(h1, 'edgecolor','none'); hold on;
%         h2=surf(Xroi,Yroi,int(i,j)*funexp2D(Xroi,Yroi,xb(i,j),yb(i,j),sx(i,j)^2,sy(i,j)^2,sxy(i,j))); set(h2, 'edgecolor','r', 'facecolor','none')
%         pause
        
        % Correction of the centroids and sigmas by fitting to a
        % 2D-Gaussian distribution:
        ellipse=[xb(i,j),yb(i,j),sx(i,j),sy(i,j),sxy(i,j)];
        rhon1=rhon; rhon1(mask==0)=0; %filter only beamlet
        [centroid,~,~,sx(i,j),sy(i,j),sxy(i,j),~,int(i,j)] = Gaussian2Dfit(rhon1, mask, ellipse); 
        if isnan(int(i,j)); int(i,j)=sum(Iroi(:)); end
        xb(i,j)=centroid(1); yb(i,j)=centroid(2);
%         figure(888); clf; h1= surf(Xroi,Yroi,Iroi); set(h1, 'edgecolor','none'); hold on;
%         h2=surf(Xroi,Yroi,int(i,j)*funexp2D(Xroi,Yroi,centroid(1),centroid(2),sx(i,j)^2,sy(i,j)^2,sxy(i,j))); set(h2, 'edgecolor','r', 'facecolor','none')
%         pause
    end
end

% figure(99); hold on; plot(xb,yb,'rx');
figure(99); hold on;
for i=1:numel(xb)
    if ~isnan(int(i)) & int(i)>0
        plot_ellipse_on_density(xb(i),yb(i),sx(i),sy(i),sxy(i));
    end
end




%% Plot intensities

figure(199); clf
if true
    imagesc(rhon); set(gca,'Xdir','Normal','Ydir','Normal'); hold on;
    for ii=1:numel(int)
        plot3(xb(ii)*[1 1],yb(ii)*[1 1],[0 int(ii)],'r-');
    end
    zlim([0 max(int(:))]);
else
    h1=surf(X,Y,rhon); set(h1, 'edgecolor','none'); hold on
    intfit=zeros(size(rhon));
    for i=1:size(xb,1)
        for j=1:size(xb,2)
            intfit=intfit+int(i,j)*funexp2D(X,Y,xb(i,j),yb(i,j),sx(i,j)^2,sy(i,j)^2,sxy(i,j));
        end
    end
    h2= surf(X,Y,intfit); set(h2, 'edgecolor','r', 'facecolor','none')
    zlim([0 max(intfit(:))]);
end
view(-30,45); xlim([min(locsx(:)) max(locsx(:))]); ylim([min(locsy(:)) max(locsy(:))]); 
xlabel('x [px]'); ylabel('y [px]'); zlabel('intensity')


%% Compute values

xvec = Xn(1, :);
yvec = Yn(:, 1);
mmppx =  mean(diff(xvec)); % meter per px. (Xn in [m])
mmppy =  mean(diff(yvec));

xg = repmat((1:size(xb,2))*mask_prop.gridSpacing,size(xb,1),1); % this gives the positions of the bars at the grid intersections [m].
yg = repmat((1:size(yb,1))*mask_prop.gridSpacing,size(yb,2),1)';

xs = interp1(1:length(xvec), xvec, xb,'linear','extrap'); % "xs, ys" will be used to compute x' and y' (xp and yp) [m]
ys = interp1(1:length(yvec), yvec, yb,'linear','extrap'); % find the positions [m] of the bar centers at the screen.

xp = (xs-xg)/mask_prop.driftLength; %convert to angles.
yp = (ys-yg)/mask_prop.driftLength;

sigmaxp = sx*mmppx/mask_prop.driftLength; % convert to [m]
sigmayp = sy*mmppy/mask_prop.driftLength;
sigmaxpyp = sxy*mmppx*mmppy/mask_prop.driftLength^2;


%%  Compute moments

x0 = sum(sum((xg) .* int))/sum(sum(int));
y0 = sum(sum((yg) .* int))/sum(sum(int));
xp0 = sum(sum((xp) .* int))/sum(sum(int));
yp0 = sum(sum((yp).* int))/sum(sum(int));

S(1,1) = sum(sum((xg-x0).^2 .* int))/sum(sum(int));
S(3,3) = sum(sum((yg-y0).^2 .* int))/sum(sum(int));

S(1,2) = sum(sum((xg-x0) .* (xp-xp0) .* int))/sum(sum(int));
S(3,4) = sum(sum((yg-y0) .* (yp-yp0) .* int))/sum(sum(int));

S(2,2) = sum(sum(((xp-xp0).^2+sigmaxp.^2) .* int))/sum(sum(int));
    S22_centroids= sum(sum(((xp-xp0).^2) .* int))/sum(sum(int));
    S22_sigmas= sum(sum((sigmaxp.^2) .* int))/sum(sum(int));
S(4,4) = sum(sum(((yp-yp0).^2+sigmayp.^2) .* int))/sum(sum(int));
    S44_centroids= sum(sum(((yp-yp0).^2) .* int))/sum(sum(int));
    S44_sigmas= sum(sum((sigmayp.^2) .* int))/sum(sum(int));

% Calculate correlations
S(1,3) = sum(sum((xg-x0) .* (yg-y0) .* int)) ./sum(sum(int));
S(1,4) = sum(sum((xg-x0) .* (yp-yp0) .* int))./sum(sum(int));
S(2,3) = sum(sum((xp-xp0) .* (yg-y0) .* int))./sum(sum(int));
S(2,4) = sum(sum(((xp-xp0).*(yp-yp0)+sigmaxpyp) .* int))./sum(sum(int));
    S24_centroids= sum(sum(((xp-xp0).*(yp-yp0)) .* int))./sum(sum(int));
    S24_sigmas= sum(sum((sigmaxpyp) .* int))./sum(sum(int));

% Symmetrize the beam matrix
S(2,1) = S(1,2);
S(3,1) = S(1,3);
S(4,1) = S(1,4);
S(3,2) = S(2,3);
S(4,2) = S(2,4);
S(4,3) = S(3,4);

% Struct with all information of analysis
info.target='PP';
info.xb=xb;
info.yb=yb;
info.sigmaxp=sigmaxp;
info.sigmayp=sigmayp;
info.sigmaxpyp=sigmaxpyp;
info.xg=xg;
info.yg=yg;
info.xs=xs;
info.ys=ys;
info.xp=xp;
info.yp=yp;
info.x0=x0;
info.y0=y0;
info.xp0=xp0;
info.yp0=yp0;
info.int=int;
info.S=S;
info.S22_centroids=S22_centroids;
info.S22_sigmas=S22_sigmas;
info.S44_centroids=S44_centroids;
info.S44_sigmas=S44_sigmas;
info.S24_centroids=S24_centroids;
info.S24_sigmas=S24_sigmas;

display('*Success')

end
