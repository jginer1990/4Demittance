function [e2Dcore, em, fraction] = core2D_distr(r_2D,S_2D,Int)

global gamma beta

maxfrac = 1.0; % maximum fraction of charge to make distrubution profile
N0=80; step=1; % start and resolution of the profile

invS_2D=inv(S_2D);

% Normalized distance with respect to center
R_2D=[];
for i=1:numel(Int)
    R_2D(i)=sqrt( r_2D(i,:)*invS_2D*r_2D(i,:)' );
end

% Sort coordinates by normalized distance
[~,isort]=sort(R_2D); 
r_2D_sorted = r_2D(isort,:);
Int0=Int(isort); 
cumInt0=cumsum(Int0); cumInt0=cumInt0/cumInt0(end);


% Build volume-charge profile arrays
i=1; k=N0;
fraction=[]; em=[];
while cumInt0(k)<maxfrac
    em(i) = gamma*beta*sqrt( abs( det( cov(r_2D_sorted(1:k,:)) ) )  );
    fraction(i) = cumInt0(k);
    k=k+step; i=i+1;
end

% Plot in figure
% figure(890); clf;
np=fraction;
Vn=em;
plot(Vn,np,'b.'); ylim([0 1]); xlim([0 Vn(end)])
xlabel('\epsilon^{2D}_{subset} [m rad]');
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
e2Dcore=Vn(end)/par(1);
plot([e2Dcore e2Dcore],[0 1],'k:');
text(e2Dcore,0,[' \epsilon^{2D}_{n,core}=' num2str(e2Dcore/1e-9,'%.1f') ' [nm rad]'],...
    'VerticalAlignment','bottom','HorizontalAlignment','left')
