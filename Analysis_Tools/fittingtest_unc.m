function [sigma,x0,intx,covx] = fittingtest_unc(x,y,barwidth,driftLength,analysis)

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

%% Discard (if applicable) points at both ends of the data to avoid overlap effects
if(isfield(analysis,'fitroi_pc'))
    fitroi_pc = analysis.fitroi_pc/100;
    roisize=length(x);
    x(round((1-(1-fitroi_pc)/2)*roisize):end)=[]; x(1:round((1-fitroi_pc)/2*roisize))=[];
    y(round((1-(1-fitroi_pc)/2)*roisize):end)=[]; y(1:round((1-fitroi_pc)/2*roisize))=[];
end


%% Fit a double erf function
% modelFunerf = @(p,x) (p(4)+p(3)*(x-p(1))).*1/2.*(erfc((x-p(1)+ barwidth/2)/(sqrt(2)*p(2)))+erfc(-(x-p(1)-barwidth/2)/(sqrt(2)*p(2)))) + ...
%     1/2*sqrt(2/pi)*p(3).*p(2).*(exp(-(x-p(1)-barwidth/2).^2./(sqrt(2)*p(2)).^2) - exp(-(x-p(1)+barwidth/2).^2./(sqrt(2)*p(2)).^2)) ;
modelFunerf = @(p,x) (p(4)+p(3)*(x-p(1))).*1/2.*(erfc((x-p(1)+ barwidth/2)/(sqrt(2)*p(2)))+erfc(-(x-p(1)-barwidth/2)/(sqrt(2)*p(2))));
startingVals = [xminval,sigma_initguess,gradinitguess,baseinitguess]; % initial guess gfit, change if necessary.

% [erfcoefs,R,J,covx,MSE,ErrorModelInfo] = nlinfit(x, y, modelFunerf, startingVals); % fit
lowb=[min(x),1,-Inf,0];
upb=[max(x),(max(x)-min(x))*0.75,Inf,2*max(y)];
[erfcoefs,~,residual,~,~,~,J] = lsqcurvefit(modelFunerf, startingVals, x, y, lowb,upb); % fit
covx=inv(J.'*J)*var(residual);

subplot(1,2,2);
plot(x,y);
hold on
plot(x,modelFunerf(erfcoefs,x),'c--',x,erfcoefs(4)+erfcoefs(3)*(x-erfcoefs(1)),'c:',erfcoefs(1),erfcoefs(4),'co');

sigma = erfcoefs(2)/driftLength;
x0 = erfcoefs(1);
intx = erfcoefs(4);

J = [...
    1 0 0 0 ; ...
    0 driftLength 0 0 ; ...
    0 0 0 1] ;
covx = J*covx*J';


%% Sanity cuts

if (x0<min(x) | x0>max(x) | driftLength*sigma> max(x)-min(x)) | intx<0
    x0=xminval;
    intx=gradinitguess*xminval+baseinitguess;
    sigma=0;
    
    if (x0<min(x) | x0>max(x)) | intx<0
        x0=0;
        intx=0;
    end
end