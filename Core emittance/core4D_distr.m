function [e4Dcore, em, fraction] = core4D_distr(r_4D,S_4D,Int)

global gamma beta

maxfrac = .9; % maximum fraction of charge to make distrubution profile
N0=100; step=100; % start and resolution of the profile

invS_4D=inv(S_4D);

% Normalized distance with respect to center
R_4D=[];
for i=1:numel(Int)
    R_4D(i)=sqrt( r_4D(i,:)*invS_4D*r_4D(i,:)' );
end

% Sort coordinates by normalized distance
[~,isort]=sort(R_4D); 
r_4D_sorted = r_4D(isort,:);
Int0=Int(isort); 
cumInt0=cumsum(Int0); cumInt0=cumInt0/cumInt0(end);


% Build volume-charge profile arrays
i=1; k=N0;
fraction=[]; em=[];
while cumInt0(k)<maxfrac
    em(i) = gamma*beta*sqrt( abs( det( cov(r_4D_sorted(1:k,:)) ) )  );
    fraction(i) = cumInt0(k);
    k=k+step; i=i+1; 
    display([num2str(cumInt0(k)/cumInt0(end)/maxfrac*100,'%.2f') '%'])
end

% Plot in figure
% figure(890); clf;
np=fraction;
Vn=em;
plot(Vn,np,'b.'); ylim([0 1]); xlim([0 Vn(end)])
xlabel('\epsilon^{4D}_{subset} [m^2 rad^2]');
ylabel('Fraction charge')
set(gca,'fontsize',12)

% Fit profile and extract slope (core brightness) and core emittance (volume at fraction=1)
roi=find(np<0.1); % Fit only first points
fun = @(m,x) m(1)*x+m(2)*x.^2; % fit function
funlin= @(m,x) m(1)*x; 
par0=[np(end) 0]; Lpar = [0 -Inf]; Upar=[Inf 0];
opts = optimset('Display','off');
[par,resid,J,~,~,~]= nlinfit(Vn(roi)/Vn(end),np(roi),fun,par0,'weights',1./np(roi).^2);
ci = nlparci(par,resid,'jacobian',J);
        
hold on; 
%plot(Vn,fun(par,Vn/Vn(end)),'r-');
plot(Vn,funlin(par(1),Vn/Vn(end)),'g-');
e4Dcore=Vn(end)/par(1);
plot([e4Dcore e4Dcore],[0 1],'k:');
text(e4Dcore,0,[' \epsilon^{4D}_{n,core}=' num2str(sqrt(e4Dcore)/1e-9,'%.1f') ' [nm rad]'],...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
