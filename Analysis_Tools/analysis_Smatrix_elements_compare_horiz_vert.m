clear keyw
keyw =  ...
    'TEM300_S1_1p2_S2_0p37_' ;
%     'PP300_S1_1p2_S2_0p37_' ;
%     'TEM500_S1_1p2_S2_0p37_' ;


clear S S_interp
element = [1 2 3 4 6 7 8 11 12 16]';
figure_number = [1 2 4 5 3 5 6 1 2 3];
element_name = {'\Sigma_{xx,yy}','\Sigma_{xx'',yy''}','\Sigma_{xy}','\Sigma_{xy'',x''y}','\Sigma_{x''x'',y''y''}',...
    '\Sigma_{xy'',x''y}','\Sigma_{x''y''}','\Sigma_{xx,yy}','\Sigma_{xx'',yy''}','\Sigma_{x''x'',y''y''}'};
figure(153); clf

searchmatfile = ['Data_' keyw '*.mat'];

filesdata = dir([folder searchmatfile]);


for i=1:length(filesdata)
    mydata=load([folder filesdata(i).name]);
    
    S{m}(:,i) = mydata.res.S(element);
    S_interp{m}(:,i) = mydata.res.S_interp(element);
    
    
end
n_set = size(S{m},2);

for j=1:10
    subplot(3,2,figure_number(j));
    plot(1:n_set,S{m}(j,:)/1e-9,'o-'); hold on;
    avg=mean(S{m}(j,:));
    err=std(S{m}(j,:));
    text(n_set,S{m}(j,end)/1e-9,[num2str(avg/1e-9,'%.1f') '\pm' num2str(err/1e-9,'%.1f')])
end


for j=1:10
    subplot(3,2,figure_number(j));
    ylabel([element_name{j} '[\times10^{-9}]']);
    %     y_lim= get(gca,'ylim');   [~,i0]=min(abs(y_lim)); y_lim(i0)=0; ylim(y_lim);
    x_lim=get(gca,'xlim'); plot(x_lim,[0 0],'k-');
end
