%% Compute average sigmaxp and sigmayp, and effect on emittance

searchmatfile = ['Data_' keywords{1} '*.mat'];

filesdata = dir([folder searchmatfile]);

clear ex ey sxp_mean syp_mean ex_cov ey_cov ex_sigma ey_sigma Cxx CxX CXX Mx Cyy CyY CYY My
    
for i=1:length(filesdata)
    mydata=load([folder filesdata(i).name]);
%     display([num2str(i) '/' num2str(length(filesdata))  ': '  mydata.fname '.tif']);
    
    gamma = mydata.saveparameters.gamma;
    beta = mydata.saveparameters.beta;
    driftLength = mydata.saveparameters.driftLength;
    pxconv=mydata.saveparameters.pxconv;
    
    xb= mydata.analysis.xb;
    yb= mydata.analysis.yb;
    xg = mydata.analysis.xg;
    yg = mydata.analysis.yg;
    xp = mydata.analysis.xp;
    yp = mydata.analysis.yp;
    try  %if target is TEM grid
        intx = mydata.analysis.intx;
        inty = mydata.analysis.inty;
    catch exc %if target is Pepper pot
        intx = mydata.analysis.int;
        inty = intx;
    end
    
    sigmaxp = mydata.analysis.sigmaxp;
    sigmayp = mydata.analysis.sigmayp;
    
    xs = xg + driftLength*xp;
    ys = yg + driftLength*yp;
    
    
    covxX(1,1)=sum(intx(:).*xg(:).^2)/sum(intx(:))-(sum(intx(:).*xg(:))/sum(intx(:)))^2;
    covxX(2,2)=sum(intx(:).*xs(:).^2)/sum(intx(:))-(sum(intx(:).*xs(:))/sum(intx(:)))^2;
    covxX(1,2)=sum(intx(:).*xg(:).*xs(:))/sum(intx(:))-(sum(intx(:).*xg(:))/sum(intx(:)))*(sum(intx(:).*xs(:))/sum(intx(:)));
    covxX(2,1)=covxX(1,2);
    
%     xs0=sum(intx(:).*xs(:))/sum(intx(:)); xg0=sum(intx(:).*xg(:))/sum(intx(:));
%     ys0=sum(inty(:).*ys(:))/sum(inty(:)); yg0=sum(inty(:).*yg(:))/sum(inty(:));
%     Mx(i) = sum(intx(:).*(xs(:)-xs0)./(xg(:)-xg0))/sum(intx(:));
%     My(i) = sum(inty(:).*(ys(:)-ys0)./(yg(:)-yg0))/sum(inty(:));
%     linfun=@(p,x) p(1)+p(2)*x;
%     pfit = nlinfit(xg(:)-xg0,xs(:)-xs0,linfun,[0 Mx(i)],'Weights',intx(:)); Mx(i)=pfit(2);
%     pfit = nlinfit(yg(:)-yg0,ys(:)-ys0,linfun,[0 My(i)],'Weights',inty(:)); My(i)=pfit(2);
    
    covyY(1,1)=sum(inty(:).*yg(:).^2)/sum(inty(:))-(sum(inty(:).*yg(:))/sum(inty(:)))^2;
    covyY(2,2)=sum(inty(:).*ys(:).^2)/sum(inty(:))-(sum(inty(:).*ys(:))/sum(inty(:)))^2;
    covyY(1,2)=sum(inty(:).*yg(:).*ys(:))/sum(inty(:))-(sum(inty(:).*yg(:))/sum(inty(:)))*(sum(inty(:).*ys(:))/sum(inty(:)));
    covyY(2,1)=covyY(1,2);
    

    Cxx(i)=covxX(1,1); CxX(i)=covxX(1,2); CXX(i)=covxX(2,2); Mx(i)=CxX(i)/Cxx(i);
    Cyy(i)=covyY(1,1); CyY(i)=covyY(1,2); CYY(i)=covyY(2,2); My(i)=CyY(i)/Cyy(i);
    
    ex(i)=mydata.res.exn/gamma/beta;
    ey(i)=mydata.res.eyn/gamma/beta;
    
    sxp_mean(i)=sum(intx(:).*sigmaxp(:))/sum(intx(:));
    syp_mean(i)=sum(inty(:).*sigmayp(:))/sum(inty(:));

    ex_cov(i) = 1/driftLength*sqrt(det(covxX));
    ey_cov(i) = 1/driftLength*sqrt(det(covyY));
    ex_sigma(i) = sqrt( covxX(1)*sum(intx(:).*sigmaxp(:).^2)/sum(intx(:)) );
    ey_sigma(i) = sqrt( covyY(1)*sum(inty(:).*sigmayp(:).^2)/sum(inty(:)) );
    
end
N= length(ex_cov);

%% plots(sigmas and emittance)
figure(650); clf
subplot(4,1,1); plot(1:N,sxp_mean/1e-6,'bo',1:N,syp_mean/1e-6,'ro'); ylabel('\sigma_{xp},\sigma_{yp} [urad]'); legend('x','y')
hold on; plot([0 N],mean(sxp_mean/1e-6)*[1 1],'b--',[0 N],mean(syp_mean/1e-6)*[1 1],'r--'); 
text(N,mean(sxp_mean/1e-6),num2str(mean(sxp_mean/1e-6))); text(N,mean(syp_mean/1e-6),num2str(mean(syp_mean/1e-6)));

subplot(4,1,2); plot(1:N,sxp_mean*driftLength/pxconv,'bo',1:N,syp_mean*driftLength/pxconv,'ro'); ylabel('\sigma_{X},\sigma_{Y} [px]')
hold on; plot([0 N],mean(sxp_mean*driftLength/pxconv)*[1 1],'b--',[0 N],mean(syp_mean*driftLength/pxconv)*[1 1],'r--'); 
text(N,mean(sxp_mean*driftLength/pxconv),num2str(mean(sxp_mean*driftLength/pxconv))); text(N,mean(syp_mean*driftLength/pxconv),num2str(mean(syp_mean*driftLength/pxconv)));

subplot(4,1,3); plot(1:N,ex/1e-9,'bo',1:N,ex_cov/1e-9,'r^',1:N,ex_sigma/1e-9,'mv'); ylabel('\epsilon_x [nm rad]'); legend('\epsilon_x','det(cov(x*X))','\epsilonx_{\sigma^2}')
hold on; plot([0 N],mean(ex/1e-9)*[1 1],'b--',[0 N],mean(ex_cov/1e-9)*[1 1],'r--',[0 N],mean(ex_sigma/1e-9)*[1 1],'m--'); 
text(N,mean(ex/1e-9),num2str(mean(ex/1e-9))); text(N,mean(ex_cov/1e-9),num2str(mean(ex_cov/1e-9))); text(N,mean(ex_sigma/1e-9),num2str(mean(ex_sigma/1e-9)));

subplot(4,1,4); plot(1:N,ey/1e-9,'bo',1:N,ey_cov/1e-9,'r^',1:N,ey_sigma/1e-9,'mv'); ylabel('\epsilon_y [nm rad]')
hold on; plot([0 N],mean(ey/1e-9)*[1 1],'b--',[0 N],mean(ey_cov/1e-9)*[1 1],'r--',[0 N],mean(ey_sigma/1e-9)*[1 1],'m--'); 
text(N,mean(ey/1e-9),num2str(mean(ey/1e-9))); text(N,mean(ey_cov/1e-9),num2str(mean(ey_cov/1e-9))); text(N,mean(ey_sigma/1e-9),num2str(mean(ey_sigma/1e-9)));

subplot(4,1,1); title(keywords{1})

%% plots(xg-xs, yg-ys correlations)

figure(651); clf
subplot(4,1,1); plot(1:N,Cxx/1e-9,'bo',1:N,Cyy/1e-9,'ro'); ylabel('\langlex^2\rangle , \langley^2\rangle'); legend('x','y')
hold on; plot([0 N],mean(Cxx/1e-9)*[1 1],'b--',[0 N],mean(Cyy/1e-9)*[1 1],'r--'); 
text(N,mean(Cxx/1e-9),num2str(mean(Cxx/1e-9))); text(N,mean(Cyy/1e-9),num2str(mean(Cyy/1e-9)));

subplot(4,1,2); plot(1:N,CXX/1e-9,'bo',1:N,CYY/1e-9,'ro'); ylabel('\langleX^2\rangle , \langleY^2\rangle');
hold on; plot([0 N],mean(CXX/1e-9)*[1 1],'b--',[0 N],mean(CYY/1e-9)*[1 1],'r--'); 
text(N,mean(CXX/1e-9),num2str(mean(CXX/1e-9))); text(N,mean(CYY/1e-9),num2str(mean(CYY/1e-9)));

subplot(4,1,3); plot(1:N,CxX/1e-9,'bo',1:N,CyY/1e-9,'ro'); ylabel('\langlexX\rangle , \langleyY\rangle');
hold on; plot([0 N],mean(CxX/1e-9)*[1 1],'b--',[0 N],mean(CyY/1e-9)*[1 1],'r--'); 
text(N,mean(CxX/1e-9),num2str(mean(CxX/1e-9))); text(N,mean(CyY/1e-9),num2str(mean(CyY/1e-9)));

subplot(4,1,3); plot(1:N,CxX/1e-9,'bo',1:N,CyY/1e-9,'ro'); ylabel('\langlexX\rangle , \langleyY\rangle');
hold on; plot([0 N],mean(CxX/1e-9)*[1 1],'b--',[0 N],mean(CyY/1e-9)*[1 1],'r--'); 
text(N,mean(CxX/1e-9),num2str(mean(CxX/1e-9))); text(N,mean(CyY/1e-9),num2str(mean(CyY/1e-9)));

subplot(4,1,4); plot(1:N,Mx,'bo',1:N,My,'ro'); ylabel('Magnification');
hold on; plot([0 N],mean(Mx)*[1 1],'b--',[0 N],mean(My)*[1 1],'r--'); 
text(N,mean(Mx),num2str(mean(Mx))); text(N,mean(My),num2str(mean(My)));

subplot(4,1,1); title(keywords{1})