function [locsx,locsy] = splitimage_PP(Asheared,analysis)
%SPLITIMAGE_PP Finds peaks and then splits image into sections around beamlets for pepperpot analysis

%% Project onto x and y axes

intenx = sum(Asheared);
intenx = intenx/max(intenx);
inteny = sum(Asheared,2);
inteny = inteny/max(inteny);

%% Find peaks

intenxsmooth = sgolayfilt(intenx,7,31);
intenysmooth = sgolayfilt(inteny,7,31);
[peaksx, locspeaksx]= findpeaks(intenxsmooth, 'minpeakdistance', analysis.minpeakdistance,'minpeakheight', analysis.minpeakheight_x);
[peaksy, locspeaksy]= findpeaks(intenysmooth', 'minpeakdistance', analysis.minpeakdistance,'minpeakheight', analysis.minpeakheight_y);

locsx=locspeaksx(1:(end-1))/2+locspeaksx(2:end)/2; 
locsx=round([locspeaksx(1)-(locsx(1)-locspeaksx(1)) locsx locspeaksx(end)+(locspeaksx(end)-locsx(end))]);
if locsx(1)<1; locsx(1)=[]; end;
if locsx(end)>length(intenx); locsx(end)=[]; end;
locsy=locspeaksy(1:(end-1))/2+locspeaksy(2:end)/2; 
locsy=round([locspeaksy(1)-(locsy(1)-locspeaksy(1)) locsy locspeaksy(end)+(locspeaksy(end)-locsy(end))]);
if locsy(1)<1; locsy(1)=[]; end;
if locsy(end)>length(inteny); locsy(end)=[]; end;


figure(9); clf
subplot(211); plot(1:length(intenx),intenx,'g-',1:length(intenxsmooth),intenxsmooth,'b-'); xlabel('x[px]'); ylabel('Projected intensity');
hold on; plot(locspeaksx,peaksx,'r^',locsx,zeros(size(locsx)),'rv'); hold on; plot([1 length(intenxsmooth)], analysis.minpeakheight_x*[1 1],'r-'); hold off
subplot(212); plot(1:length(inteny),inteny,'g-',1:length(intenysmooth),intenysmooth,'b-'); xlabel('y[px]'); ylabel('Projected intensity');
hold on; plot(locspeaksy,peaksy,'r^',locsy,zeros(size(locsy)),'rv'); hold on; plot([1 length(intenysmooth)], analysis.minpeakheight_y*[1 1],'r-'); hold off
    
%% Plot divided screen image
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
end

