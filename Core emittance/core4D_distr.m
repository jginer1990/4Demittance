function [e4Dcore, em, fraction] = core4D_distr(r_4D,S_4D,Int)
num=numel(Int);


plot_gaussian=true;
maxfrac = .8; % maximum fraction of charge to make distrubution profile

full_em=sqrt(det(S_4D));
invS_4D=inv(S_4D);

% Normalized distance with respect to center
R_4D=zeros(num,1);
progress=5;
for i=1:num
    R_4D(i)=sqrt( r_4D(i,:)*invS_4D*r_4D(i,:)' );
    if i/num*100>=progress
        display(['4D core emittance (1/2):  ' num2str(progress) '%'])
        progress=progress+5;
    end
end

% Sort coordinates by normalized distance
[~,isort]=sort(R_4D); 
r_4D_sorted = r_4D(isort,:);
Int0=Int(isort); Int0 = Int0/sum(Int0);
cumInt0=cumsum(Int0); %cumInt0=cumInt0/cumInt0(end);


N0=100; npoints=150; % start and resolution of the profile
% step=round(find(cumInt0<maxfrac,1,'last')/npoints)*ones(npoints,1);  % use linear spacing
step=(0:npoints).^2; step=step/npoints^2*(find(cumInt0<maxfrac,1,'last')-N0)+N0; step=round(diff(step)); % use quadratic spacing


% Build volume-charge profile arrays
i=1; k=N0;
fraction=zeros(npoints,1); em=fraction; fraction_gaussian=fraction; em_gaussian=fraction;
progress=1;
while cumInt0(k)<maxfrac
    em(i) = emittance_fraction(r_4D_sorted(1:k,:),Int0(1:k));
    fraction(i) = cumInt0(k);

    if plot_gaussian
        fraction_gaussian(i) = 1-(1+R_4D(isort(k))^2/2)*exp(-R_4D( isort(k) )^2/2);
        em_gaussian(i)=full_em*(1- R_4D(isort(k))^4/8/(exp(R_4D(isort(k))^2/2)-(1+R_4D(isort(k))^2/2)))^2;
    end
    k=k+step(min(i,npoints)); i=i+1; 
    if cumInt0(k)/maxfrac*100>=progress
        display(['4D core emittance (2/2):  ' num2str(progress) '%'])
        progress=progress+1;
    end
end
fraction(fraction==0)=[]; em(fraction==0)=[];

% Plot in figure
% figure(890); clf;
np=fraction;
Vn=em;
if plot_gaussian; plot(em_gaussian,fraction_gaussian,'m--','linewidth',2); hold on; plot([0 full_em*2/9],[0 1],'m:','linewidth',1);end
plot(Vn,np,'b.-','linewidth',2); ylim([0 1]); xlim([0 Vn(end)])
xlabel('\epsilon^{4D}_{subset} [m^2 rad^2]');
ylabel('Fraction charge')
set(gca,'fontsize',12)

% Fit profile and extract slope (core brightness) and core emittance (volume at fraction=1)
roi=find(np<0.1); % Fit only first points
if length(roi)<20; roi=1:20; end;
fun = @(m,x) m(1)*x+m(2)*x.^2; % fit function
funlin= @(m,x) m(1)*x; 
par0=[np(roi(end)) 0]; Lpar = [0 -Inf]; Upar=[Inf 0];
opts = optimset('Display','off');
% [par,resid,J,~,~,~]= nlinfit(Vn(roi)/Vn(roi(end)),np(roi),fun,par0,'weights',1./np(roi).^2);
[par,~,resid,~,~,~,J]= lsqcurvefit(fun,par0,Vn(roi)/Vn(roi(end)),np(roi),Lpar,Upar,opts);
ci = nlparci(par,resid,'jacobian',J);
        
hold on; 
%plot(Vn,fun(par,Vn/Vn(end)),'r-');
plot(Vn,funlin(par(1),Vn/Vn(roi(end))),'g-');
e4Dcore=Vn(roi(end))/par(1);
plot([e4Dcore e4Dcore],[0 1],'k:');
text(e4Dcore,0,[' (\epsilon^{4D}_{core})^{1/2}=' num2str(sqrt(e4Dcore)/1e-9,'%.1f') ' [nm rad]'],...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
