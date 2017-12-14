function [xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp] = phasespace_TEM(rhon,Xn,Yn,locsx,locsy,Locsx,Locsy,mask_prop,analysis)
%PHASESPACE_TEM analyse image and compute parameters

avpeaksepx = mean(diff(locsx));
avpeaksepy = mean(diff(locsy));

barwidth_x_pix = round(avpeaksepx/mask_prop.pitch_to_bar_width_ratio); % barwidth at screen in pixels
barwidth_y_pix = round(avpeaksepy/mask_prop.pitch_to_bar_width_ratio);

%% Loop over y widths to determine parameters in x
xb=zeros(length(locsy)-1,length(locsx)); sx=xb; intx=xb; ybc=xb;
for j=1:length(locsy)-1; % Loop over y
    intenx = sum(rhon(locsy(j):locsy(j+1),:)); % Select intensity values between two horizontal bars
%    intenx = intenx/max(intenx); % Normalize intensity

    %% Plot intensity
    figure(100); clf;
    subplot(length(locsy),1,j)
    plot(1:length(intenx), intenx); % Plot projections and troughs

    %% Fit erf in x
    for i=1:length(locsx); % Loop over x
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
yb=zeros(length(locsy),length(locsx)-1); sy=yb; inty=yb; xbc=yb;
for j=1:length(locsx)-1; % Loop over x
    temp = rhon(:,locsx(j):locsx(j+1));
    inteny = sum(temp,2);
    %% Plot intensity
    figure(200)
    subplot(length(locsx),1,j)
    plot(1:length(inteny), inteny); % Plot projections and troughs
    %% Fit erf in y
    for i=1:length(locsy); % Loop over y
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

end
