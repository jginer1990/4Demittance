clear keyw
keyw = keywords{1};

clear S S_interp
element = [1 2 3 4 6 7 8 11 12 16]';
element_name = {'\langle x^2\rangle [m^2]','\langle xx''\rangle [m rad]',...
    '\langle xy\rangle [m^2]','\langle xy''\rangle [m rad]','\langle x''^2\rangle [rad^2]',...
    '\langle x''y\rangle [m rad]','\langle x''y''\rangle [rad^2]','\langle y^2\rangle [m^2]',...
    '\langle yy''\rangle [m rad]','\langle y''^2\rangle [rad^2]'};
figure(153); clf
    
    searchmatfile = ['Data_' keyw '*.mat'];
    
    filesdata = dir([folder searchmatfile]);
    
    
    for i=1:length(filesdata)
        mydata=load([folder filesdata(i).name]);
        
        S(:,i) = mydata.res.S(element);
        S_interp(:,i) = mydata.res.S_interp(element);
        
        
    end
    n_set = size(S,2);
    
    for j=1:10
        subplot(5,2,j); 
        plot(1:n_set,S(j,:),'o-'); hold on;
        avg=mean(S(j,:));
        err=std(S(j,:));
        text(n_set,S(j,end),[num2str(avg/1e-9,'%.1f') '\pm' num2str(err/1e-9,'%.1f')])
    end
    


for j=1:10
    subplot(5,2,j);
    ylabel(element_name{j});
%     y_lim= get(gca,'ylim');   [~,i0]=min(abs(y_lim)); y_lim(i0)=0; ylim(y_lim);
    x_lim=get(gca,'xlim'); plot(x_lim,[0 0],'k-');
end
% legend(keyw)