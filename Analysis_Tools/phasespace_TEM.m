function [S,info] = phasespace_TEM(rhon,Xn,Yn,locsx,locsy,Locsx,Locsy,minsx,minsy,mask_prop,analysis)
%PHASESPACE_TEM analyse image and compute parameters

avpeaksepx = mean(mean(diff(Locsx,1,2)));
avpeaksepy = mean(mean(diff(Locsy,1,1)));

barwidth_x_pix = round(avpeaksepx/mask_prop.pitch_to_bar_width_ratio); % barwidth at screen in pixels
barwidth_y_pix = round(avpeaksepy/mask_prop.pitch_to_bar_width_ratio);

ypos = locsy;
xpos = locsx;

%% Loop over y widths to determine parameters in x
xb=zeros(length(ypos)-1,length(xpos)); sx=xb; intx=xb; ybc=xb;
for j=1:length(ypos)-1; % Loop over y
    intenx = sum(rhon(ypos(j):ypos(j+1),:)); % Select intensity values between two horizontal bars
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
        regiony = yin-round(0.3*avpeaksepy):yin+round(analysis.interval_pc/2/100*avpeaksepy); % Define region along y with 60% of the points between bars
        regionx(regionx<1 | regionx>size(rhon,2))=[]; regiony(regiony<1 | regiony>size(rhon,1))=[];
        if any(regionx<0) || any(regionx>size(intenx,2))
            flag = 1;
        end
        if flag == 0
            Ixroi = sum(rhon(regiony,regionx),1);
            try
                [sx(j,i), xb(j,i), intx(j,i)] = fittingtest(regionx, Ixroi, barwidth_x_pix, mask_prop.driftLength); % Fit erf of region between two peaks
            catch
                Ixroi = sgolayfilt(Ixroi,7,21);
                [sx(j,i), xb(j,i), intx(j,i)] = fittingtest(regionx, Ixroi, barwidth_x_pix, mask_prop.driftLength); % Fit erf of region between two peaks
            end
        end
    end
end        

%% Loop over x widths to determine parameters in y
yb=zeros(length(ypos),length(xpos)-1); sy=yb; inty=yb; xbc=yb;
for j=1:length(xpos)-1; % Loop over x
    temp = rhon(:,xpos(j):xpos(j+1));
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
        regionx = xin-round(0.3*avpeaksepx):xin+round(analysis.interval_pc/2/100*avpeaksepx);
        regionx(regionx<1 | regionx>size(rhon,2))=[]; regiony(regiony<1 | regiony>size(rhon,1))=[];
        if any(regiony<0) || any(regiony>size(inteny,1))
            flag = 1;
        end
        if flag == 0
            Iyroi = sum(rhon(regiony,regionx),2);
            try
                [sy(i,j), yb(i,j), inty(i,j)] = fittingtest(regiony, Iyroi', barwidth_y_pix, mask_prop.driftLength); % Fit erf of region between two peaks.
            catch
                Iyroi = sgolayfilt(Iyroi,7,21);
                [sy(i,j), yb(i,j), inty(i,j)] = fittingtest(regiony, Iyroi', barwidth_y_pix, mask_prop.driftLength); % Fit erf of region between two peaks.
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

xg = repmat((1:size(xb,2))*mask_prop.gridSpacing,size(xb,1),1); % this gives the positions of the bars at the grid intersections [m].
yg = repmat((1:size(yb,1))*mask_prop.gridSpacing,size(yb,2),1)';

xs = interp1(1:length(xvec), xvec, xb,'linear','extrap'); % "xs, ys" will be used to compute x' and y' (xp and yp) [m]
ys = interp1(1:length(yvec), yvec, yb,'linear','extrap'); % find the positions [m] of the bar centers at the screen.

xp = (xs-xg)/mask_prop.driftLength; %convert to angles.
yp = (ys-yg)/mask_prop.driftLength;

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


% Struct with all information of analysis
info.target='TEM';
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
info.x0=x0;
info.y0=y0;
info.x0c=x0c;
info.y0c=y0c;
info.xpcen=xpcen;
info.ypcen=ypcen;
info.xp0=xp0;
info.yp0=yp0;
info.xp0c=xp0c;
info.yp0c=yp0c;
info.intx=intx;
info.inty=inty;
info.xcen=xcen;
info.ycen=ycen;
info.intxcen=intxcen;
info.intycen=intycen;
info.intcen=intcen;
info.S=S;

display('*Success')

end
