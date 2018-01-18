function [sigma,x0,intx] = fittingtest(x,y,barwidth,driftLength)

%% Search for minimum
[minval, ~] = min(y);
xminval = x(y==minval);
        
%% Check for multiple minima
if length(xminval)~=1 && mean(diff(xminval))<10
    xminval = mean(xminval); %If there are multiple minima, and they are close, just average them.
elseif length(xminval) ~=1 && min(minval)~=0
    error('Multiple minima, some nonzero.'); % if not, something's wrong.
end
subplot(1, 2,1)
hold on
plot(x,y)
plot(xminval,minval,'*')

%% Fit a line to determine overall slope of beam shape
xstart = mean(x(1:4));
xend = mean(x(end-3:end));
ystart = mean(y(1:4));
yend = mean(y(end-3:end));

gradinitguess = (yend-ystart)/(xend-xstart);

%% Fit a double erf function
%modelFunerf = @(p,x) p(4)+ p(3).*(-erf((x-p(1)+ p(5))/(sqrt(2)*driftLength*p(2)))+erf((x-p(1)-p(5))/(sqrt(2)*driftLength*p(2)))); % double erf model function.
modelFunerf = @(p,x) p(5)+p(4)*x+ p(3).*(-erf((x-p(1)+ barwidth/2)/(sqrt(2)*driftLength*p(2)))+erf((x-p(1)-barwidth/2)/(sqrt(2)*driftLength*p(2))));
startingVals = [xminval,5,max(y),gradinitguess,0]; % initial guess gfit, change if necessary.

erfcoefs = nlinfit(x, y, modelFunerf, startingVals); % fit
subplot(1,2,2);
plot(x,y);
hold on
plot(x,modelFunerf(erfcoefs,x),'c--');

sigma = erfcoefs(2);
x0 = erfcoefs(1);
intx = erfcoefs(5)+erfcoefs(4)*xminval;

if (sigma < 0 | intx<0)
    sigma = 5;
    x0 = xminval;
    intx = 0;
end


%% Sanity cuts

if (x0<min(x) | x0>max(x) | sigma> max(x)-min(x))
    x0=0;
    intx=0;
end