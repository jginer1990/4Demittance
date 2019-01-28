
searchmatfile = ['Data_' keywords{1} '*.mat']; 

filesdata = dir([folder searchmatfile]);

clear ea eb exn eyn e1n e2n en4D exn_interp eyn_interp e1n_interp e2n_interp en4D_interp exn2Dcore eyn2Dcore en4Dcore 
for i=1:length(filesdata)
    mydata=load([folder filesdata(i).name]); 
    
    exn(i)=mydata.res.exn;
    eyn(i)=mydata.res.eyn;
    e1n(i)=mydata.res.e1n;
    e2n(i)=mydata.res.e2n;
    en4D(i)=sqrt(abs(det(mydata.res.S)))*(gamma*beta)^2;
    exn_interp(i)=mydata.res.exn_interp;
    eyn_interp(i)=mydata.res.eyn_interp;
    e1n_interp(i)=mydata.res.e1n_interp;
    e2n_interp(i)=mydata.res.e2n_interp;
    en4D_interp(i)=sqrt(abs(det(mydata.res.S_interp)))*(gamma*beta)^2;
    exn2Dcore(i)=mydata.res.exn2Dcore;
    eyn2Dcore(i)=mydata.res.eyn2Dcore;
    en4Dcore(i)=mydata.res.en4Dcore;
%     counts(i)=mydata.counts;
end
n_set = length(exn);

%%
figure(152); clf

select = [1 2 7 5 6];
% select = [1 2 7];

for i =1:length(select)

plotemittance = select(i);
switch plotemittance
    case 1
        ea = exn;
        eb = eyn;
        aleg='\epsilon_{n,x}'; bleg='\epsilon_{n,y}';
        ti='Projected emittance (no interpolation)'; 
    case 2
%         ea = e1n;
%         eb = e2n;      
        ea = min(e1n,e2n);
        eb = max(e1n,e2n);  
        aleg='\epsilon_{n,1}'; bleg='\epsilon_{n,2}';
        ti='Intrinsic emittance (no interpolation)';
    case 3
        ea = exn_interp;
        eb = eyn_interp;        
        aleg='\epsilon_{n,x}'; bleg='\epsilon_{n,y}';
        ti='Projected emittance (interpolation)';
    case 4
        ea = e1n_interp;
        eb = e2n_interp;
        aleg='\epsilon_{n,1}'; bleg='\epsilon_{n,2}';
        ti='Intrinsic emittance (interpolation)';
    case 5
        ea = exn2Dcore;
        eb = eyn2Dcore;
        aleg='\epsilon^{2D}_{n,x}'; bleg='\epsilon^{2D}_{n,y}';
        ti='2D Core Emittance';
    case 6
        ea = sqrt(exn2Dcore.*eyn2Dcore);
        eb = sqrt(en4Dcore);
        aleg='(\epsilon^{2D}_{n,x}\cdot\epsilon^{2D}_{n,y})^{1/2}'; bleg='(\epsilon^{4D}_{n})^{1/2}';
        ti='2D vs 4D Core Emittance';
    case 7
        ea = sqrt(exn.*eyn);
        eb = sqrt(en4D);
        aleg='(\epsilon^{2D}_{n,x}\cdot\epsilon^{2D}_{n,y})^{1/2}'; bleg='(\epsilon^{4D}_{n})^{1/2}';
        ti='2D vs 4D Emittance (no interpolation)';
    case 8
        ea = sqrt(exn_interp.*eyn_interp);
        eb = sqrt(en4D_interp);
        aleg='(\epsilon^{2D}_{n,x}\cdot\epsilon^{2D}_{n,y})^{1/2}'; bleg='(\epsilon^{4D}_{n})^{1/2}';
        ti='2D vs 4D Core Emittance (interpolation)';
    case 9
        ea = counts/1e15;
        eb = zeros(size(counts));
        aleg='charge counts'; bleg='';
        ti='charge';
    otherwise
        break;
end

subplot(length(select),1,i)
plot(1:n_set,ea/1e-9,'b>'); hold on;
plot(1:n_set,eb/1e-9,'r^');
limy=get(gca,'ylim'); limy(1)=0; ylim(limy);
xlim([0 n_set+1])
onlyreal= find(real(ea)~=0 & ~isnan(ea));
plot([0 n_set+1],mean(ea(onlyreal))/1e-9*[1 1],'b--')
text(0,mean(ea(onlyreal))/1e-9,...
    [num2str(mean(ea(onlyreal))/1e-9,'%.1f') '\pm' num2str(std(ea(onlyreal))/1e-9,'%.1f')]...
    ,'VerticalAlignment','top')
onlyreal= find(real(eb)~=0 & ~isnan(eb));
plot([0 n_set+1],mean(eb(onlyreal))/1e-9*[1 1],'r--')
text(0,mean(eb(onlyreal))/1e-9,...
    [num2str(mean(eb(onlyreal))/1e-9,'%.1f') '\pm' num2str(std(eb(onlyreal))/1e-9,'%.1f')]...
    ,'VerticalAlignment','top')

xlabel('Image #')
ylabel('Normalized emittance [nm rad]')
legend(aleg,bleg,'location','southeast')
title(ti)
end