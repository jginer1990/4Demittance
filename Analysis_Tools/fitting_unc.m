function [sigma,x0,intx,covx] = fitting_unc(x,y,barwidth,driftLength,analysis,divline)

if(isfield(analysis,'sigma_initguess'))
    sigma_initguess = analysis.sigma_initguess;
else
    sigma_initguess = 5;
end

%% Search for minimum
[minval, ~] = min(y);
xminval = x(y==minval);
        
%% Check for multiple minima
if length(xminval)~=1 && mean(diff(xminval))<10
    xminval = mean(xminval); %If there are multiple minima, and they are close, just average them.
elseif length(xminval) ~=1 && min(minval)~=0
    error('Multiple minima, some nonzero.'); % if not, something's wrong.
end
subplot(1,2,1)
hold on
plot(x,y)
plot(xminval,minval,'*')

%% Fit a line to determine overall slope of beam shape
xstart = mean(x(1:4));
xend = mean(x(end-3:end));
ystart = mean(y(1:4));
yend = mean(y(end-3:end));

gradinitguess = (yend-ystart)/(xend-xstart);
baseinitguess = ystart-gradinitguess*xstart;

%% Fit a double erf function

warning('') % Clear last warning message

modelFunerf = @(p,x) p(5)+p(4)*x+ p(3).*(-erf((x-p(1)+ barwidth/2)/(sqrt(2)*driftLength*p(2)))+erf((x-p(1)-barwidth/2)/(sqrt(2)*driftLength*p(2))));
startingVals = [xminval,sqrt(sigma_initguess),max(y),gradinitguess,baseinitguess]; % initial guess gfit, change if necessary.

[erfcoefs,R,J,covx,MSE,ErrorModelInfo] = nlinfit(x, y, modelFunerf, startingVals); % fit
% covx = covx(1:3,1:3);
subplot(1,2,2);
plot(x,y);
hold on
plot(x,modelFunerf(erfcoefs,x),'c--');

sigma = erfcoefs(2);
x0 = erfcoefs(1);
intx = erfcoefs(5)+erfcoefs(4)*x0;

J = [...
    1 0 0 0 0; ...
    0 1 0 0 0; ...
    erfcoefs(4) 0 0 erfcoefs(1) 1] ;
covx = J*covx*J';


%% Sanity cuts

[warnMsg, warnId] = lastwarn;

if (x0<min(x) | x0>max(x) | driftLength*sigma> max(x)-min(x)) | intx<0 | ~isempty(warnMsg);
    x0=divline;
    intx=gradinitguess*x0+baseinitguess;
    sigma=0;
    covx = zeros(size(covx));
    
    if (x0<min(x) | x0>max(x)) | intx<0
        x0=divline;
        intx=0;
    end
end