function [S,info] = phasespace_shear_PP(rhon,X,Y,Xn,Yn,X0px,Y0px,S1,S2,mask_prop,locsx,locsy)
%PHASESPACE_SHEAR analyse image and compute parameters


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

shear_fourier;

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

intenxsmooth = sgolayfilt(intenx,7,61);
intenysmooth = sgolayfilt(inteny,7,61);
[minsx, locsx]= findpeaks(-intenxsmooth, 'minpeakdistance', minpeakdistance,'minpeakheight', maxtroughheight_x);
[minsy, locsy]= findpeaks(-intenysmooth', 'minpeakdistance', minpeakdistance,'minpeakheight', maxtroughheight_y);

idx = zeros(size(minsx,2),1);
idy = zeros(size(minsy,2),1);
idx(1) =0;
idx(end) = 0;
idy(1) = 0;
idy(end) = 0;
totcounts = sum(intenx);

for j = 2:size(minsx,2)-1
    if (sum(intenx(locsx(j):locsx(j+1))) < threshint*totcounts) | abs(locsx(j+1)-locsx(j))>maxpeakdistance
        idx(j) = 0;
    else
        idx(j) = 1;
        jlast=j+1;
    end
end
idx(jlast) = 1;

totcounts = sum(inteny);
for j = 2:size(minsy,2)-1
    if (sum(inteny(locsy(j):locsy(j+1))) < threshint*totcounts) | abs(locsy(j+1)-locsy(j))>maxpeakdistance
        idy(j) = 0;
    else
        idy(j) = 1;
        jlast=j+1;
    end
end
idy(jlast) =1;

minsx(~idx) = [];
locsx(~idx) = [];
minsy(~idy) = [];
locsy(~idy) = [];

figure(9); clf
subplot(211); plot(1:length(intenx),-intenx,'g-',1:length(intenxsmooth),-intenxsmooth,'b-'); xlabel('x[px]'); ylabel('Projected intensity');
hold on; plot(locsx,minsx,'r^'); hold on; plot([1 length(intenxsmooth)], maxtroughheight_x*[1 1],'r-'); hold off
subplot(212); plot(1:length(inteny),-inteny,'g-',1:length(intenysmooth),-intenysmooth,'b-'); xlabel('y[px]'); ylabel('Projected intensity');
hold on; plot(locsy,minsy,'r^'); hold on; plot([1 length(intenysmooth)], maxtroughheight_y*[1 1],'r-'); hold off
    
% Plot divided screen image
figure(2)
set(gcf,'Name','Screen sheared image')
imagesc(Asheared); axis('xy')
siz = size(Asheared);
hold on
for ii=1:length(minsx)
    plot(locsx(ii)*[1 1], [1 siz(1)], 'y'); 
end

for ii=1:length(minsy)
    plot([1 siz(2)],locsy(ii)*[1 1], 'y');
end
hold off

if length(minsx)<3 || length(minsy)<3
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


%% Determine pixel width of bars on screen

avpeaksepx = mean(mean(diff(Locsx,1,2)));
avpeaksepy = mean(mean(diff(Locsy,1,1)));

barwidth_x_pix = round(avpeaksepx/pitch_to_bar_width_ratio); % barwidth at screen in pixels
barwidth_y_pix = round(avpeaksepy/pitch_to_bar_width_ratio);

ypos = locsy;
xpos = locsx;

%% Loop over y widths to determine parameters in x
xb=zeros(length(ypos)-1,length(xpos)); sx=xb; intx=xb; ybc=xb;
for j=1:length(ypos)-1; % Loop over y
    intenx = sum(Asheared(ypos(j):ypos(j+1),:)); % Select intensity values between two horizontal bars
%    intenx = intenx/max(intenx); % Normalize intensity

    %% Plot intensity
    figure(100); clf;
    subplot(length(ypos),1,j)
    plot(1:length(intenx), intenx, xpos, -minsx, 'o'); % Plot projections and troughs

    %% Fit erf in x
    for i=1:length(xpos); % Loop over x
        flag = 0;
        xin=round((Locsx(j,i)+Locsx(j+1,i))/2);
        yin=round((Locsy(j,i)+Locsy(j+1,i))/2); ybc(j,i)=yin;
        regionx = xin-round(0.5*avpeaksepx):xin+round(0.5*avpeaksepx); % Define region from midpoint to midpoint in x
        regiony = yin-round(0.3*avpeaksepy):yin+round(interval_pc/2/100*avpeaksepy); % Define region along y with 60% of the points between bars
        regionx(regionx<1 | regionx>size(rhon,2))=[]; regiony(regiony<1 | regiony>size(rhon,1))=[];
        if any(regionx<0) || any(regionx>size(intenx,2))
            flag = 1;
        end
        if flag == 0
            Ixroi = sum(rhon(regiony,regionx),1);
            try
                [sx(j,i), xb(j,i), intx(j,i)] = fittingtest(regionx, Ixroi, barwidth_x_pix, driftLength); % Fit erf of region between two peaks
            catch
                Ixroi = sgolayfilt(Ixroi,7,21);
                [sx(j,i), xb(j,i), intx(j,i)] = fittingtest(regionx, Ixroi, barwidth_x_pix, driftLength); % Fit erf of region between two peaks
            end
        end
    end
end        

%% Loop over x widths to determine parameters in y
yb=zeros(length(ypos),length(xpos)-1); sy=yb; inty=yb; xbc=yb;
for j=1:length(xpos)-1; % Loop over x
    temp = Asheared(:,xpos(j):xpos(j+1));
    inteny = sum(temp,2);
    %% Plot intensity
    figure(200)
    subplot(length(xpos),1,j)
    plot(1:length(inteny), inteny, locsy, -minsy, 'o'); % Plot projections and troughs
    %% Fit erf in y
    for i=1:length(ypos); % Loop over y
        flag = 0;
        xin=round((Locsx(i,j)+Locsx(i,j+1))/2); xbc(i,j)=xin;
        yin=round((Locsy(i,j)+Locsy(i,j+1))/2);
        regiony = yin-round(0.5*avpeaksepy):yin+round(0.5*avpeaksepy);
        regionx = xin-round(0.3*avpeaksepx):xin+round(interval_pc/2/100*avpeaksepx);
        regionx(regionx<1 | regionx>size(rhon,2))=[]; regiony(regiony<1 | regiony>size(rhon,1))=[];
        if any(regiony<0) || any(regiony>size(inteny,1))
            flag = 1;
        end
        if flag == 0
            Iyroi = sum(rhon(regiony,regionx),2);
            try
                [sy(i,j), yb(i,j), inty(i,j)] = fittingtest(regiony, Iyroi', barwidth_y_pix, driftLength); % Fit erf of region between two peaks.
            catch
                Iyroi = sgolayfilt(Iyroi,7,21);
                [sy(i,j), yb(i,j), inty(i,j)] = fittingtest(regiony, Iyroi', barwidth_y_pix, driftLength); % Fit erf of region between two peaks.
            end
        end
    end 
end


%% Plot locations

figure(99); hold on;  plot(xb(intx~=0),ybc(intx~=0),'r>'); plot(xbc(inty~=0),yb(inty~=0),'r^')

clear ybmid ybcen xbmid xbcen
ybmid(:,2:size(yb,2)+1) = yb(:,:);
ybmid(:,1) = ybmid(:,2);
ybmid(:,end+1) = ybmid(:,end);
ybcen = (ybmid(:,1:end-1)+ybmid(:,2:end)) /2;

xbmid(:,2:size(xb,1)+1) = xb(:,:)';
xbmid(:,1) = xbmid(:,2);
xbmid(:,end+1) = xbmid(:,end);
xbcen = ((xbmid(:,1:end-1)+xbmid(:,2:end))/2)';

%% Plot intensities

figure(199); clf
imagesc(rhon); set(gca,'Xdir','Normal','Ydir','Normal'); hold on;
for ii=1:numel(intx)
    plot3(xb(ii)*[1 1],ybc(ii)*[1 1],[0 intx(ii)],'y-');
end
for ii=1:numel(inty)
    plot3(xbc(ii)*[1 1],yb(ii)*[1 1],[0 inty(ii)],'r-');
end
view(-30,45)
xlabel('x [px]'); ylabel('y [px]'); zlabel('intensity')

%% Find xb positions based on minima

% xb_mins = zeros(length(ypos)-1,length(xposcent)-1);
% for j=1:length(ypos)-1; % Loop over y
%     intenx = sum(rhon(ypos(j):ypos(j+1),:)); % Select intensity values between two horizontal bars
%     intenx = intenx/max(intenx); % Normalize intensity
%     for i=1:length(xposcent)-1; % Loop over x
%         regionx = xposcent(i):xposcent(i+1);
%         Ixroi = intenx(regionx);
%         Ixroi = sgolayfilt(Ixroi,7,21);
%         [~,minindex] = min(Ixroi);
%         xb_mins(j,i) = regionx(minindex);
%     end    
% end    
% 
% %% Find yb positions based on minima
% yb_mins = zeros(length(yposcent)-1,length(xpos)-1);
% for j=1:length(xpos)-1; 
%     temp = rhon(:,xpos(j):xpos(j+1));
%     inteny = sum(temp,2);
%     inteny = inteny/max(inteny);
%     for i=1:length(yposcent)-1; 
%         regiony = yposcent(i):yposcent(i+1);
%         Iyroi = inteny(regiony);
%         Iyroi = sgolayfilt(Iyroi,7,21);
%         [~,minindex] = min(Iyroi);
%         yb_mins(i,j) = regiony(minindex);
%     end    
% end    

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

sigmaxp = sx*mmppx; % convert to [m] (we take into account the shear transformation)
sigmayp = sy*mmppy;


%%  Compute moments

x0 = sum(sum((xg) .* intx))/sum(sum(intx));
y0 = sum(sum((yg) .* inty))/sum(sum(inty));
xp0 = sum(sum((xp) .* intx))/sum(sum(intx));
yp0 = sum(sum((yp).* inty))/sum(sum(inty));

S(1,1) = sum(sum((xg-x0).^2 .* intx))/sum(sum(intx));
S(3,3) = sum(sum((yg-y0).^2 .* inty))/sum(sum(inty));

S(1,2) = sum(sum((xg-x0) .* (xp-xp0) .* intx))/sum(sum(intx));
S(3,4) = sum(sum((yg-y0) .* (yp-yp0) .* inty))/sum(sum(inty));

S(2,2) = sum(sum(((xp-xp0).^2+sigmaxp.^2) .* intx))/sum(sum(intx));
S(4,4) = sum(sum(((yp-yp0).^2+sigmayp.^2) .* inty))/sum(sum(inty));

% Interpolate vertex values from midpoint values

clear ymid xmid intymid intxmid intxcen intycen ycen xcen
clear ypmid xpmid ypcen xpcen
ymid(:,2:size(yb,2)+1) = yg(:,:);
ymid(:,1) = ymid(:,2);
ymid(:,end+1) = ymid(:,end);
ycen = (ymid(:,1:end-1)+ymid(:,2:end)) /2;

ypmid(:,2:size(yb,2)+1) = yp(:,:);
ypmid(:,1) = ypmid(:,2);
ypmid(:,end+1) = ypmid(:,end);
ypcen = (ypmid(:,1:end-1)+ypmid(:,2:end)) /2;

intymid(:,2:size(yb,2)+1) = inty(:,:);
intymid(:,1) = intymid(:,2);
intymid(:,end+1) = intymid(:,end);
intycen = (intymid(:,1:end-1)+intymid(:,2:end)) /2 .*sign(intymid(:,1:end-1).*intymid(:,2:end)); % if any is zero, make intcen=0

xmid(:,2:size(xb,1)+1) = xg(:,:)';
xmid(:,1) = xmid(:,2);
xmid(:,end+1) = xmid(:,end);
xcen = ((xmid(:,1:end-1)+xmid(:,2:end))/2)';

xpmid(:,2:size(xb,1)+1) = xp(:,:)';
xpmid(:,1) = xpmid(:,2);
xpmid(:,end+1) = xpmid(:,end);
xpcen = ((xpmid(:,1:end-1)+xpmid(:,2:end))/2)';

intxmid(:,2:size(xb,1)+1) = intx(:,:)';
intxmid(:,1) = intxmid(:,2);
intxmid(:,end+1) = intxmid(:,end);
intxcen = ((intxmid(:,1:end-1)+intxmid(:,2:end))/2.*sign(intxmid(:,1:end-1).*intxmid(:,2:end)))'; % if any is zero, make intcen=0

intcen = sqrt(intycen.*intxcen);
figure(99); hold on; plot(xbcen(intcen~=0),ybcen(intcen~=0),'go'); 
figure(199); hold on
for ii=1:numel(intcen)
    plot3(xbcen(ii)*[1 1],ybcen(ii)*[1 1],[0 intcen(ii)],'g-');
end

x0c = sum(sum(xcen.*intcen))./sum(sum(intcen));
y0c = sum(sum(ycen.*intcen))./sum(sum(intcen));
xp0c = sum(sum(xpcen.*intcen))./sum(sum(intcen));
yp0c = sum(sum(ypcen.*intcen))./sum(sum(intcen));

% Calculate correlations
S(1,3) = sum(sum((xcen-x0c) .* (ycen-y0c) .* intcen))        ./sum(sum(intcen));
S(1,4) = sum(sum((xcen-x0c) .* (ypcen-yp0c) .* intcen))    ./sum(sum(intcen));
S(2,3) = sum(sum((xpcen-xp0c) .* (ycen-y0c) .* intcen))    ./sum(sum(intcen));
S(2,4) = sum(sum((xpcen-xp0c) .* (ypcen-yp0c) .* intcen))./sum(sum(intcen));

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
info.ybc=ybc;
info.yb=yb;
info.xbc=xbc;
info.xbcen=xbcen;
info.ybcen=ybcen;
info.sigmaxp=sigmaxp;
info.sigmayp=sigmayp;
info.sigmaxpyp=0;
info.xg=xg;
info.yg=yg;
info.xs=xs;
info.ys=ys;
info.xp=xp;
info.yp=yp;
info.intx=intx;
info.inty=inty;
info.intxcen=intxcen;
info.intycen=intycen;
info.intcen=intcen;
info.S=S;
% info.S22_centroids=S22_centroids;
% info.S22_sigmas=S22_sigmas;
% info.S44_centroids=S44_centroids;
% info.S44_sigmas=S44_sigmas;
% info.S24_centroids=S24_centroids;
% info.S24_sigmas=S24_sigmas;