
%% Determines 4D phase space
rhon = A;
Xn = X;
Yn = Y;

[X,Y]=meshgrid(1:size(rhon,2),1:size(rhon,1));
X0px=sum(rhon(:).*X(:))/sum(rhon(:));
Y0px=sum(rhon(:).*Y(:))/sum(rhon(:));

%% Find angles using Radon transform
% try
%     theta = 1:220;
%     [Rad,xp] = radon(rhon,theta);
%     xpFilt = xp(xp>-100&xp<100);
%     RadFilt = Rad(xp>-100&xp<100,:); %RadFilt(RadFilt<0)=0;
%     % figure(91); imagesc(theta,xpFilt,RadFilt); set(gca,'YDir','normal')
%     fcon=max(RadFilt)-min(RadFilt);
%     % figure(92); plot(theta,fcon)
%     [radpeaks, radlocstheta]= findpeaks(fcon,'MinPeakDistance',30);
%     [sortedValues, sortedIndices] = sort(radpeaks, 'descend');
%     radlocsthetaSorted = radlocstheta(sortedIndices);
%     if 180-abs(radlocsthetaSorted(1)-radlocsthetaSorted(2))<10
%         distance(1) = min(radlocsthetaSorted(1)-theta(1),theta(end)-radlocsthetaSorted(1));
%         distance(2) = min(radlocsthetaSorted(2)-theta(1),theta(end)-radlocsthetaSorted(2));
%         [~,i]= max(distance); twoangles(1)=radlocsthetaSorted(i);
%         twoangles(2)=radlocsthetaSorted(3);
%     else
%         twoangles=radlocsthetaSorted(1:2);
%     end
%     
%     for i=1:2
%         range= twoangles(i)-5:twoangles(i)+5;
%         par= polyfit(range,fcon(range),2); twoangles_fit(i)=-par(2)/2/par(1);
%     end
%     
%     angle2 = twoangles_fit(twoangles_fit>=45&twoangles_fit<135);
%     angle1 = twoangles_fit(twoangles_fit~=angle2);
%     S1 = -tand(angle1);
%     S2 = tand(angle2-90);
%     
%     Asheared = shear_image(X,Y,rhon,S1,S2,X0px,Y0px);
%     figure(2)
%     imagesc(Asheared)
%     set(gca,'YDir','normal')
%     
%     errorshear = false;
% catch
%     errorshear = true;
%     display('!Radon transformation failed.')
% end

%% Find angles with 2D-Fourier transform

shear_fourier_v2; 


%% Shear image manually
if errorshear
    manual = 1;
else
    manual = 0;%input('*Find angles manually? (0 or 1): ');
end 
if manual
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
    
    Asheared = shear_image(X,Y,rhon,S1,S2,X0px,Y0px);
    figure(2)
    imagesc(Asheared)
    set(gca,'YDir','normal')
elseif errorshear
    Asheared = rhon;
    S1 = 0; S2 = 0;
end

%% Project onto x and y axes

intenx = sum(Asheared);
intenx = intenx/max(intenx);
inteny = sum(Asheared,2);
inteny = inteny/max(inteny);

%% Determine trough locations

intenxsmooth = sgolayfilt(intenx,7,31);
intenysmooth = sgolayfilt(inteny,7,31);
[peaksx, locspeaksx]= findpeaks(intenxsmooth, 'minpeakdistance', minpeakdistance,'minpeakheight', minpeakheight_x);
[peaksy, locspeaksy]= findpeaks(intenysmooth', 'minpeakdistance', minpeakdistance,'minpeakheight', minpeakheight_y);

locsx=locspeaksx(1:(end-1))/2+locspeaksx(2:end)/2; 
locsx=round([locspeaksx(1)-(locsx(1)-locspeaksx(1)) locsx locspeaksx(end)+(locspeaksx(end)-locsx(end))]);
if locsx(1)<1; locsx(1)=[]; end;
if locsx(end)>length(intenx); locsx(end)=[]; end;
locsy=locspeaksy(1:(end-1))/2+locspeaksy(2:end)/2; 
locsy=round([locspeaksy(1)-(locsy(1)-locspeaksy(1)) locsy locspeaksy(end)+(locspeaksy(end)-locsy(end))]);
if locsy(1)<1; locsy(1)=[]; end;
if locsy(end)>length(inteny); locsy(end)=[]; end;

% intenxsmooth = sgolayfilt(intenx,7,61);
% intenysmooth = sgolayfilt(inteny,7,61);
% [minsx, locsx]= findpeaks(-intenxsmooth, 'minpeakdistance', minpeakdistance,'minpeakheight', maxtroughheight_x);
% [minsy, locsy]= findpeaks(-intenysmooth', 'minpeakdistance', minpeakdistance,'minpeakheight', maxtroughheight_y);

% idx = zeros(size(locsx,2),1);
% idy = zeros(size(locsy,2),1);
% idx(1) =0;
% idx(end) = 0;
% idy(1) = 0;
% idy(end) = 0;
% totcounts = sum(intenx);
% 
% for j = 2:size(locsx,2)-1
%     if (sum(intenx(locsx(j):locsx(j+1))) < threshint*totcounts) | abs(locsx(j+1)-locsx(j))>maxpeakdistance
%         idx(j) = 0;
%     else
%         idx(j) = 1;
%         jlast=j+1;
%     end
% end
% idx(jlast) = 1;
%
% totcounts = sum(inteny);
% for j = 2:size(locsy,2)-1
%     if (sum(inteny(locsy(j):locsy(j+1))) < threshint*totcounts) | abs(locsy(j+1)-locsy(j))>maxpeakdistance
%         idy(j) = 0;
%     else
%         idy(j) = 1;
%         jlast=j+1;
%     end
% end
% idy(jlast) =1;
% 
% peaksx(~idx) = [];
% locsx(~idx) = [];
% peaksy(~idy) = [];
% locsy(~idy) = [];

figure(9); clf
subplot(211); plot(1:length(intenx),intenx,'g-',1:length(intenxsmooth),intenxsmooth,'b-'); xlabel('x[px]'); ylabel('Projected intensity');
hold on; plot(locspeaksx,peaksx,'r^',locsx,zeros(size(locsx)),'rv'); hold on; plot([1 length(intenxsmooth)], minpeakheight_x*[1 1],'r-'); hold off
subplot(212); plot(1:length(inteny),inteny,'g-',1:length(intenysmooth),intenysmooth,'b-'); xlabel('y[px]'); ylabel('Projected intensity');
hold on; plot(locspeaksy,peaksy,'r^',locsy,zeros(size(locsy)),'rv'); hold on; plot([1 length(intenysmooth)], minpeakheight_y*[1 1],'r-'); hold off
    
% Plot divided screen image
figure(2)
set(gcf,'Name','Screen sheared image')
imagesc(Asheared); axis('xy')
siz = size(Asheared);
hold on
for ii=1:length(locsx)
    plot(locsx(ii)*[1 1], [1 siz(1)], 'y'); 
end

for ii=1:length(locsy)
    plot([1 siz(2)],locsy(ii)*[1 1], 'y');
end
hold off

if length(locsx)<3 || length(locsy)<3
    error('Not enough troughs in data.')
end

%% Undo shear in Locsx,Locsy

[Locsx_sh,Locsy_sh] = meshgrid(locsx,locsy);

Locsx=X0px+1/(1-S1*S2)*((Locsx_sh-X0px) -S1*(Locsy_sh-Y0px));
Locsy=Y0px+1/(1-S1*S2)*(-S2*(Locsx_sh-X0px) +(Locsy_sh-Y0px));

figure(99); clf
imagesc(rhon); 
set(gcf,'Name','Screen image')
set(gca,'YDir','normal')
hold on
for i=1:size(Locsx,1)
    plot([Locsx(i,1) Locsx(i,end)],[Locsy(i,1) Locsy(i,end)],'y-')
end
for i=1:size(Locsx,2)
    plot([Locsx(1,i) Locsx(end,i)],[Locsy(1,i) Locsy(end,i)],'y-')
end


%% Determine beamlets centroids and rms ellipse parameters

funexp2D=@(x,y,mux,muy,sxx,syy,sxy)...
    1/(2*pi*sqrt(sxx*syy-sxy^2))*exp(-1/(2*(1-sxy^2/sxx/syy))*((x-mux).^2/sxx+(y-muy).^2/syy-2*sxy/sxx/syy*(x-mux).*(y-muy)));

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
                Iroi(ii,jj) = rhon(Yroi(ii,jj),Xroi(ii,jj));
                mask(Yroi(ii,jj),Xroi(ii,jj)) = 1;
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

figure(99); hold on; plot(xb,yb,'rx');




%% Plot intensities

figure(199); clf
if false
    imagesc(rhon); set(gca,'Xdir','Normal','Ydir','Normal'); hold on;
    for ii=1:numel(int)
        plot3(xb(ii)*[1 1],yb(ii)*[1 1],[0 int(ii)],'r-');
    end
else
    h1=surf(X,Y,rhon); set(h1, 'edgecolor','none'); hold on
    intfit=zeros(size(rhon));
    for i=1:size(xb,1)
        for j=1:size(xb,2)
            intfit=intfit+int(i,j)*funexp2D(X,Y,xb(i,j),yb(i,j),sx(i,j)^2,sy(i,j)^2,sxy(i,j));
        end
    end
    h2= surf(X,Y,intfit); set(h2, 'edgecolor','r', 'facecolor','none')
end
view(-30,45); zlim([0 max(intfit(:))]); xlim([min(Locsx(:)) max(Locsx(:))]); ylim([min(Locsy(:)) max(Locsy(:))]); 
xlabel('x [px]'); ylabel('y [px]'); zlabel('intensity')


%% Compute values

xvec = Xn(1, :);
yvec = Yn(:, 1);
mmppx =  mean(diff(xvec)); % meter per px. (Xn in [m])
mmppy =  mean(diff(yvec));

xg = repmat((1:size(xb,2))*gridSpacing,size(xb,1),1); % this gives the positions of the bars at the grid intersections [m].
yg = repmat((1:size(yb,1))*gridSpacing,size(yb,2),1)';

xs = interp1(1:length(xvec), xvec, xb,'linear','extrap'); % "xs, ys" will be used to compute x' and y' (xp and yp) [m]
ys = interp1(1:length(yvec), yvec, yb,'linear','extrap'); % find the positions [m] of the bar centers at the screen.

xp = (xs-xg)/driftLength; %convert to angles.
yp = (ys-yg)/driftLength;

sigmaxp = sx*mmppx/driftLength; % convert to [m]
sigmayp = sy*mmppy/driftLength;
sigmaxpyp = sxy*mmppx*mmppy/driftLength^2;


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

% Interpolate vertex values from midpoint values

% Symmetrize the beam matrix
S(2,1) = S(1,2);
S(3,1) = S(1,3);
S(4,1) = S(1,4);
S(3,2) = S(2,3);
S(4,2) = S(2,4);
S(4,3) = S(3,4);

display('*Success')
% S

% Struct with all information of analysis
info.target=target;
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
info.int=int;
info.S=S;
info.S22_centroids=S22_centroids;
info.S22_sigmas=S22_sigmas;
info.S44_centroids=S44_centroids;
info.S44_sigmas=S44_sigmas;
info.S24_centroids=S24_centroids;
info.S24_sigmas=S24_sigmas;
