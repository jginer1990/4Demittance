function [sigma,x0,intx] = fitting_lowemit(x,y,barwidth,driftLength,analysis)

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

%% Fit a double erf function
sigma = 0;
x0 = xminval;
intx = mean([ystart yend]);


%% Sanity cuts

if (x0<min(x) | x0>max(x)| intx<0)
    x0=0;
    intx=0;
end