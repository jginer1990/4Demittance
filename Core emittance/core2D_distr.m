function [e2Dcore, em, fraction] = core2D_distr(r_2D,S_2D,Int)

plot_gaussian=true;

maxfrac = .99; % maximum fraction of charge to make distrubution profile
N0=50; step=1; % start and resolution of the profile

full_em = sqrt(det(S_2D));
invS_2D=inv(S_2D);

% Normalized distance with respect to center
R_2D=zeros(numel(Int),1);
for i=1:numel(Int)
    R_2D(i)=sqrt( r_2D(i,:)*invS_2D*r_2D(i,:)' );
end

% Sort coordinates by normalized distance
[~,isort]=sort(R_2D); 
r_2D_sorted = r_2D(isort,:);
Int0=Int(isort); Int0 = Int0/sum(Int0);
cumInt0=cumsum(Int0); %cumInt0=cumInt0/cumInt0(end);


% Build volume-charge profile arrays
i=1; k=N0;
fraction=[]; em=[]; fraction_gaussian=[];
while cumInt0(k)<maxfrac
    em(i) = emittance_fraction(r_2D_sorted(1:k,:),Int0(1:k));
    fraction(i) = cumInt0(k);
    if plot_gaussian
        fraction_gaussian(i) = 1-exp(-R_2D( isort(k) )^2/2);
        em_gaussian(i)=full_em*(1- R_2D(isort(k))^2/2/(exp(R_2D(isort(k))^2/2)-1));
    end
    k=k+step; i=i+1;
end



% Plot in figure
% figure(890); clf;
np=fraction;
Vn=em;
if plot_gaussian; plot(em_gaussian,fraction_gaussian,'m--','linewidth',2); hold on; plot([0 full_em/2],[0 1],'m:','linewidth',1);end
plot(Vn,np,'b.'); ylim([0 1]); xlim([0 Vn(end)]); 
xlabel('\epsilon^{2D}_{subset} [m rad]');
ylabel('Fraction charge')
set(gca,'fontsize',12)

hold on; 
% Fit profile and extract slope (core brightness) and core emittance (volume at fraction=1)
roi=find(np<0.1); % Fit only first points
if length(roi)<10; roi=1:10; end;
fun = @(m,x) m(1)*x+m(2)*x.^2; % fit function
funlin= @(m,x) m(1)*x; 
par0=[np(roi(end)) 0]; Lpar = [0 -Inf]; Upar=[Inf 0];
opts = optimset('Display','off');
% [par,resid,J,~,~,~]= nlinfit(Vn(roi)/Vn(roi(end)),np(roi),fun,par0,'weights',1./np(roi).^2);
[par,~,resid,~,~,~,J]= lsqcurvefit(fun,par0,Vn(roi)/Vn(roi(end)),np(roi),Lpar,Upar,opts);
ci = nlparci(par,resid,'jacobian',J);
        
%plot(Vn,fun(par,Vn/Vn(end)),'r-');
plot(Vn,funlin(par(1),Vn/Vn(roi(end))),'g-');
e2Dcore=Vn(roi(end))/par(1);
plot([e2Dcore e2Dcore],[0 1],'k:');
text(e2Dcore,0,[' \epsilon^{2D}_{core}=' num2str(e2Dcore/1e-9,'%.1f') ' [nm rad]'],...
    'VerticalAlignment','bottom','HorizontalAlignment','left')