function [locsx,locsy,minsx,minsy] = splitimage_TEM(Asheared,analysis)
%SPLITIMAGE_TEM Finds peaks and then splits image into sections around beamlets for TEM grid analysis

if(isfield(analysis,'splitsmoothparam'))
    sp = analysis.splitsmoothparam;
else
    sp = 61;
end


%% Project onto x and y axes

intenx = sum(Asheared);
intenx = intenx/max(intenx);
inteny = sum(Asheared,2);
inteny = inteny/max(inteny);

%% Determine trough locations


intenxsmooth = sgolayfilt(intenx,7,sp);
intenysmooth = sgolayfilt(inteny,7,sp);
[minsx, locsx]= findpeaks(-intenxsmooth, 'minpeakdistance', analysis.minpeakdistance,'minpeakheight', analysis.maxtroughheight_x);
[minsy, locsy]= findpeaks(-intenysmooth', 'minpeakdistance', analysis.minpeakdistance,'minpeakheight', analysis.maxtroughheight_y);

idx = zeros(size(minsx,2),1);
idy = zeros(size(minsy,2),1);
idx(1) =0;
idx(end) = 0;
idy(1) = 0;
idy(end) = 0;
totcounts = sum(intenx);

for j = 2:size(minsx,2)-1
    if (sum(intenx(locsx(j):locsx(j+1))) < analysis.threshint*totcounts) | abs(locsx(j+1)-locsx(j))>analysis.maxpeakdistance
        idx(j) = 0;
    else
        idx(j) = 1;
        jlast=j+1;
    end
end
idx(jlast) = 1;

totcounts = sum(inteny);
for j = 2:size(minsy,2)-1
    if (sum(inteny(locsy(j):locsy(j+1))) < analysis.threshint*totcounts) | abs(locsy(j+1)-locsy(j))>analysis.maxpeakdistance
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

figure(92); clf
subplot(211); plot(1:length(intenx),-intenx,'g-',1:length(intenxsmooth),-intenxsmooth,'b-'); xlabel('x[px]'); ylabel('Projected intensity');
hold on; plot(locsx,minsx,'r^'); hold on; plot([1 length(intenxsmooth)], analysis.maxtroughheight_x*[1 1],'r-'); hold off
subplot(212); plot(1:length(inteny),-inteny,'g-',1:length(intenysmooth),-intenysmooth,'b-'); xlabel('y[px]'); ylabel('Projected intensity');
hold on; plot(locsy,minsy,'r^'); hold on; plot([1 length(intenysmooth)], analysis.maxtroughheight_y*[1 1],'r-'); hold off


% Plot divided screen image
figure(2)
set(gcf,'Name','Screen sheared image')
imagesc(Asheared); axis('xy')
siz = size(Asheared);
hold on
for ii=1:length(minsx)
    plot(locsx(ii)*[1 1], [1 siz(1)], 'y:'); 
end

for ii=1:length(minsy)
    plot([1 siz(2)],locsy(ii)*[1 1], 'y:');
end
plot(1:length(intenxsmooth),intenxsmooth*100,'r-')
plot(intenysmooth*100,1:length(intenysmooth),'r-')
hold off

if length(minsx)<3 || length(minsy)<3
    error('Not enough troughs in data.')
end
end

