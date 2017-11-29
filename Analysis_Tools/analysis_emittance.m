%% Image rotation to remove x-y correlation

for i=1:length(results)
    
    xyangle = 1/2*atand(2*results(i).S(1,3)/(results(i).S(1,1)-results(i).S(3,3)));
    xyangle_interp = 1/2*atand(2*results(i).S_interp(1,3)/(results(i).S_interp(1,1)-results(i).S_interp(3,3)));
    
    
    Mrot= [...
        cosd(xyangle) 0 sind(xyangle) 0; ...
        0 cosd(xyangle) 0 sind(xyangle); ...
        -sind(xyangle) 0 cosd(xyangle) 0; ...
        0 -sind(xyangle) 0 cosd(xyangle)];
    Srot = Mrot*results(i).S*Mrot';
    
    
    Mrot_interp= [...
        cosd(xyangle_interp) 0 sind(xyangle_interp) 0; ...
        0 cosd(xyangle_interp) 0 sind(xyangle_interp); ...
        -sind(xyangle_interp) 0 cosd(xyangle_interp) 0; ...
        0 -sind(xyangle_interp) 0 cosd(xyangle_interp)];
    Srot_interp = Mrot_interp*results(i).S_interp*Mrot_interp';
    
    results(i).angle = xyangle;
    results(i).angle_interp = xyangle_interp;
    results(i).Srot = Srot;
    results(i).Srot_interp = Srot_interp;
    results(i).exrot = gamma*beta*sqrt(det(Srot(1:2,1:2)));
    results(i).eyrot = gamma*beta*sqrt(det(Srot(3:4,3:4)));
    results(i).exrot_interp = gamma*beta*sqrt(det(Srot_interp(1:2,1:2)));
    results(i).eyrot_interp = gamma*beta*sqrt(det(Srot_interp(3:4,3:4)));
    
    
end

%% Plot emittance results

figure(152); clf
plotemittance = 3;

clear ex ey

for i=1:length(results)
    
    switch plotemittance
        case 1
            ex(i) = results(i).ex;
            ey(i) = results(i).ey;
            ti='Projected emittance (no interpolation)';
        case 2
            ex(i) = results(i).e1;
            ey(i) = results(i).e2;
            ti='Intrinsic emittance (no interpolation)';
        case 3
            ex(i) = results(i).ex_interp;
            ey(i) = results(i).ey_interp;
            ti='Projected emittance (interpolation)';
        case 4
            ex(i) = results(i).e1_interp;
            ey(i) = results(i).e2_interp;
            ti='Intrinsic emittance (interpolation)';
        case 5
            ex(i) = results(i).exrot;
            ey(i) = results(i).eyrot;
            ti='Rotated beam emittance (no interpolation)';
        case 6
            ex(i) = results(i).exrot_interp;
            ey(i) = results(i).eyrot_interp;    
            ti='Rotated beam emittance (interpolation)';
        otherwise
            break;
    end
end

plot(1:length(ex),ex/1e-9,'b>'); hold on;
plot(1:length(ey),ey/1e-9,'r^');
limy=get(gca,'ylim'); limy(1)=0; ylim(limy);
xlim([0 length(ex)+1])
plot([0 length(ex)+1],mean(ex(real(ex)~=0 & ~isnan(ex)))/1e-9*[1 1],'b--')
plot([0 length(ex)+1],mean(ey(real(ey)~=0 & ~isnan(ex)))/1e-9*[1 1],'r--')
text(0,mean(ex(real(ex)~=0 & ~isnan(ex)))/1e-9,num2str(mean(ex(real(ex)~=0 & ~isnan(ex)))/1e-9,'%.1f'),'VerticalAlignment','top')
text(0,mean(ey(real(ey)~=0 & ~isnan(ex)))/1e-9,num2str(mean(ey(real(ey)~=0 & ~isnan(ex)))/1e-9,'%.1f'),'VerticalAlignment','top')

xlabel('Image #')
ylabel('Normalized emittance [nm rad]')
legend('Horizontal','Vertical','location','southeast')
title(ti)
return

%% Emittance vs Charge

figure(152); clf
plotemittance = 3;

clear ex ey Q

for i=1:length(results)
    Q(i) = results(i).charge;
    switch plotemittance
        case 1
            ex(i) = results(i).ex;
            ey(i) = results(i).ey;
            ti='Projected emittance (no interpolation)';
        case 2
            ex(i) = results(i).e1;
            ey(i) = results(i).e2;
            ti='Intrinsic emittance (no interpolation)';
        case 3
            ex(i) = results(i).ex_interp;
            ey(i) = results(i).ey_interp;
            ti='Projected emittance (interpolation)';
        case 4
            ex(i) = results(i).e1_interp;
            ey(i) = results(i).e2_interp;
            ti='Intrinsic emittance (interpolation)';
        case 5
            ex(i) = results(i).exrot;
            ey(i) = results(i).eyrot;
            ti='Rotated beam emittance (no interpolation)';
        case 6
            ex(i) = results(i).exrot_interp;
            ey(i) = results(i).eyrot_interp;    
            ti='Rotated beam emittance (interpolation)';
        otherwise
            break;
    end
end

plot(Q,ex/1e-9,'b>'); hold on;
plot(Q,ey/1e-9,'r^');
% limy=get(gca,'ylim'); limy(1)=0; ylim(limy);

% xlabel('Image #')
ylabel('Normalized emittance [nm rad]')
legend('Horizontal','Vertical','location','southeast')
title(ti)
return

%% Analysis beam size, divergence and correlation coefficients
clear Sigx Sigxp Sigy Sigyp corrS12 corrS13 corrS14 corrS23 corrS24 corrS34
clear Sigx_interp Sigxp_interp Sigy_interp Sigyp_interp corrS12_interp corrS13_interp corrS14_interp corrS23_interp corrS24_interp corrS34_interp
for i = 1:length(results)
    auxS=results(i).Srot;
    Sigx(i) = sqrt(auxS(1,1));
    Sigxp(i) = sqrt(auxS(2,2));
    Sigy(i) = sqrt(auxS(3,3));
    Sigyp(i) = sqrt(auxS(4,4));
    corrS12(i)= auxS(1,2)/Sigx(i)/Sigxp(i);
    corrS13(i)= auxS(1,3)/Sigx(i)/Sigy(i);
    corrS14(i)= auxS(1,4)/Sigx(i)/Sigyp(i);
    corrS23(i)= auxS(2,3)/Sigxp(i)/Sigy(i);
    corrS24(i)= auxS(2,4)/Sigxp(i)/Sigyp(i);
    corrS34(i)= auxS(3,4)/Sigy(i)/Sigyp(i);

    auxS_interp=results(i).Srot_interp;
    Sigx_interp(i) = sqrt(auxS_interp(1,1));
    Sigxp_interp(i) = sqrt(auxS_interp(2,2));
    Sigy_interp(i) = sqrt(auxS_interp(3,3));
    Sigyp_interp(i) = sqrt(auxS_interp(4,4));
    corrS12_interp(i)= auxS_interp(1,2)/Sigx_interp(i)/Sigxp_interp(i);
    corrS13_interp(i)= auxS_interp(1,3)/Sigx_interp(i)/Sigy_interp(i);
    corrS14_interp(i)= auxS_interp(1,4)/Sigx_interp(i)/Sigyp_interp(i);
    corrS23_interp(i)= auxS_interp(2,3)/Sigxp_interp(i)/Sigy_interp(i);
    corrS24_interp(i)= auxS_interp(2,4)/Sigxp_interp(i)/Sigyp_interp(i);
    corrS34_interp(i)= auxS_interp(3,4)/Sigy_interp(i)/Sigyp_interp(i);
    
end

figure(153); clf

subplot(421);
plot(1:length(results),Sigx/1e-3,1:length(results),Sigy/1e-3); ylabel('\sigma_x \sigma_y [mm]')
subplot(423);
plot(1:length(results),Sigxp/1e-3,1:length(results),Sigyp/1e-3); ylabel('\sigma_{xp} \sigma_{yp} [mrad]')
subplot(425);
plot(1:length(results),corrS12,1:length(results),corrS34);
ylabel('Correlation coeff'); legend('\rho_{x,xp}','\rho_{y,yp}')
subplot(427);
plot(1:length(results),corrS13,1:length(results),corrS14,...
    1:length(results),corrS23,1:length(results),corrS24);
ylabel('Correlation coeff'); legend('\rho_{x,y}','\rho_{x,yp}','\rho_{xp,y}','\rho_{xp,yp}')

subplot(422);
plot(1:length(results),Sigx_interp/1e-3,1:length(results),Sigy_interp/1e-3); ylabel('\sigma_x \sigma_y [mm]')
subplot(424);
plot(1:length(results),Sigxp_interp/1e-3,1:length(results),Sigyp_interp/1e-3); ylabel('\sigma_{xp} \sigma_{yp} [mrad]')
subplot(426);
plot(1:length(results),corrS12_interp,1:length(results),corrS34_interp);
ylabel('Correlation coeff'); legend('\rho_{x,xp}','\rho_{y,yp}')
subplot(428);
plot(1:length(results),corrS13_interp,1:length(results),corrS14_interp,...
    1:length(results),corrS23_interp,1:length(results),corrS24_interp);
ylabel('Correlation coeff'); legend('\rho_{x,y}','\rho_{x,yp}','\rho_{xp,y}','\rho_{xp,yp}')


%% Plot all beam matrix elements (COMPARE PEPPERPOT AND TEM)
load('\\10.0.0.49\pegasus\Pegasusdata\2017_11_03\Lightfield\TEM300_ S1_1.24_S2_0.49_counts_6p7e4.mat')
clear axx axpxp axxp ayy aypyp ayyp axy axyp axpy axpyp aex aey
for i = 1:length(results)
    auxS=results(i).S;
    axx(i) = auxS(1,1);
    axpxp(i) = auxS(2,2);
    axxp(i) = auxS(1,2);
    ayy(i) = auxS(3,3);
    aypyp(i) = auxS(4,4);
    ayyp(i) = auxS(3,4);
    axy(i) = auxS(1,3);
    axyp(i) = auxS(1,4);
    axpy(i) = auxS(2,3);
    axpyp(i) = auxS(2,4);
    aex(i)= results(i).ex_interp;
    aey(i)= results(i).ey_interp;
end

figure(154); clf
subplot(621); hold on; plot(1:length(results),axx,'bo-','linewidth',1); plot([0 length(results)+1],mean(axx)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xx} [m^2]')
subplot(622); hold on; plot(1:length(results),ayy,'bo-','linewidth',1); plot([0 length(results)+1],mean(ayy)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{yy} [m^2]')
subplot(623); hold on; plot(1:length(results),axpxp,'bo-','linewidth',1); plot([0 length(results)+1],mean(axpxp)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xpxp} [rad^2]')
subplot(624); hold on; plot(1:length(results),aypyp,'bo-','linewidth',1); plot([0 length(results)+1],mean(aypyp)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{ypyp} [rad^2]')
subplot(625); hold on; plot(1:length(results),axxp,'bo-','linewidth',1); plot([0 length(results)+1],mean(axxp)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xxp} [m\cdotrad]')
subplot(626); hold on; plot(1:length(results),ayyp,'bo-','linewidth',1); plot([0 length(results)+1],mean(ayyp)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{yyp} [m\cdotrad]')
subplot(627); hold on; plot(1:length(results),axyp,'bo-','linewidth',1); plot([0 length(results)+1],mean(axyp)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xyp} [m\cdotrad]')
subplot(628); hold on; plot(1:length(results),axpy,'bo-','linewidth',1); plot([0 length(results)+1],mean(axpy)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xpy} [m\cdotrad]')
subplot(629); hold on; plot(1:length(results),axy,'bo-','linewidth',1); plot([0 length(results)+1],mean(axy)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xy} [m^2]');
subplot(6,2,10); hold on; plot(1:length(results),axpyp,'bo-','linewidth',1); plot([0 length(results)+1],mean(axpyp)*[1 1],'b:','linewidth',1); ylabel('\Sigma_{xpyp} [rad^2]');
subplot(6,2,11); hold on; plot(1:length(results),aex,'bo-','linewidth',1); plot([0 length(results)+1],mean(aex)*[1 1],'b:','linewidth',1); ylabel('\epsilon_{n,x} [m\cdotrad]'); xlabel('Image #')
subplot(6,2,12); hold on; plot(1:length(results),aey,'bo-','linewidth',1); plot([0 length(results)+1],mean(aey)*[1 1],'b:','linewidth',1); ylabel('\epsilon_{n,y} [m\cdotrad]'); xlabel('Image #')

load('\\10.0.0.49\pegasus\Pegasusdata\2017_11_03\Lightfield\PepperPotCu300_ S1_1.24_S2_0.49_counts_6p7e4.mat')
clear axx axpxp axxp ayy aypyp ayyp axy axyp axpy axpyp aex aey
for i = 1:length(results)
    auxS=results(i).S;
    axx(i) = auxS(1,1);
    axpxp(i) = auxS(2,2);
    axxp(i) = auxS(1,2);
    ayy(i) = auxS(3,3);
    aypyp(i) = auxS(4,4);
    ayyp(i) = auxS(3,4);
    axy(i) = auxS(1,3);
    axyp(i) = auxS(1,4);
    axpy(i) = auxS(2,3);
    axpyp(i) = auxS(2,4);
    aex(i)= results(i).ex_interp;
    aey(i)= results(i).ey_interp;
end

figure(154); %clf
subplot(621); hold on; plot(1:length(results),axx,'r^-','linewidth',1); plot([0 length(results)+1],mean(axx)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xx} [m^2]')
subplot(622); hold on; plot(1:length(results),ayy,'r^-','linewidth',1); plot([0 length(results)+1],mean(ayy)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{yy} [m^2]')
subplot(623); hold on; plot(1:length(results),axpxp,'r^-','linewidth',1); plot([0 length(results)+1],mean(axpxp)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xpxp} [rad^2]')
subplot(624); hold on; plot(1:length(results),aypyp,'r^-','linewidth',1); plot([0 length(results)+1],mean(aypyp)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{ypyp} [rad^2]')
subplot(625); hold on; plot(1:length(results),axxp,'r^-','linewidth',1); plot([0 length(results)+1],mean(axxp)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xxp} [m\cdotrad]')
subplot(626); hold on; plot(1:length(results),ayyp,'r^-','linewidth',1); plot([0 length(results)+1],mean(ayyp)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{yyp} [m\cdotrad]')
subplot(627); hold on; plot(1:length(results),axyp,'r^-','linewidth',1); plot([0 length(results)+1],mean(axyp)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xyp} [m\cdotrad]')
subplot(628); hold on; plot(1:length(results),axpy,'r^-','linewidth',1); plot([0 length(results)+1],mean(axpy)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xpy} [m\cdotrad]')
subplot(629); hold on; plot(1:length(results),axy,'r^-','linewidth',1); plot([0 length(results)+1],mean(axy)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xy} [m^2]');
subplot(6,2,10); hold on; plot(1:length(results),axpyp,'r^-','linewidth',1); plot([0 length(results)+1],mean(axpyp)*[1 1],'r:','linewidth',1); ylabel('\Sigma_{xpyp} [rad^2]'); 
subplot(6,2,11); hold on; plot(1:length(results),aex,'r^-','linewidth',1); plot([0 length(results)+1],mean(aex)*[1 1],'r:','linewidth',1); ylabel('\epsilon_{n,x} [m\cdotrad]'); xlabel('Image #')
subplot(6,2,12); hold on; plot(1:length(results),aey,'r^-','linewidth',1); plot([0 length(results)+1],mean(aey)*[1 1],'r:','linewidth',1); ylabel('\epsilon_{n,y} [m\cdotrad]'); xlabel('Image #')


subplot(621); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(622); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(623); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(624); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(625); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(626); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(627); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(628); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(629); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(6,2,10); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(6,2,11); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);
subplot(6,2,12); limy=get(gca,'ylim'); limy(1)=0; limy(2)=limy(2); ylim(limy);



%%

load('\\10.0.0.49\pegasus\Pegasusdata\2017_11_03\Lightfield\TEM300_ S1_1.24_S2_0.49_counts_6p7e4.mat')
clear axx axpxp axxp ayy aypyp ayyp axy axyp axpy axpyp 
for i=1:length(results)
    auxS=results(i).S_interp;
    axx(i) = auxS(1,1);
    axpxp(i) = auxS(2,2);
    axxp(i) = auxS(1,2);
    ayy(i) = auxS(3,3);
    aypyp(i) = auxS(4,4);
    ayyp(i) = auxS(3,4);
    axy(i) = auxS(1,3);
    axyp(i) = auxS(1,4);
    axpy(i) = auxS(2,3);
    axpyp(i) = auxS(2,4);
end

figure(155); clf
subplot(231); plot(axx,axpxp,'bx'); xlabel('\Sigma_{xx}'); ylabel('\Sigma_{xpxp}');
subplot(232); plot(axx,axxp,'bx'); xlabel('\Sigma_{xx}'); ylabel('\Sigma_{xxp}');
subplot(233); plot(axxp,axpxp,'bx'); xlabel('\Sigma_{xxp}'); ylabel('\Sigma_{xpxp}');
subplot(234); plot(ayy,aypyp,'bx'); xlabel('\Sigma_{yy}'); ylabel('\Sigma_{ypyp}');
subplot(235); plot(ayy,ayyp,'bx'); xlabel('\Sigma_{yy}'); ylabel('\Sigma_{yyp}');
subplot(236); plot(ayyp,aypyp,'bx'); xlabel('\Sigma_{yyp}'); ylabel('\Sigma_{ypyp}');

load('\\10.0.0.49\pegasus\Pegasusdata\2017_11_03\Lightfield\PepperPotCu300_ S1_1.24_S2_0.49_counts_6p7e4.mat')
clear axx axpxp axxp ayy aypyp ayyp axy axyp axpy axpyp 
for i=1:length(results)
    auxS=results(i).S_interp;
    axx(i) = auxS(1,1);
    axpxp(i) = auxS(2,2);
    axxp(i) = auxS(1,2);
    ayy(i) = auxS(3,3);
    aypyp(i) = auxS(4,4);
    ayyp(i) = auxS(3,4);
    axy(i) = auxS(1,3);
    axyp(i) = auxS(1,4);
    axpy(i) = auxS(2,3);
    axpyp(i) = auxS(2,4);
end

figure(155); 
subplot(231); hold on; plot(axx,axpxp,'r+'); xlabel('\Sigma_{xx}'); ylabel('\Sigma_{xpxp}');
subplot(232); hold on; plot(axx,axxp,'r+'); xlabel('\Sigma_{xx}'); ylabel('\Sigma_{xxp}');
subplot(233); hold on; plot(axxp,axpxp,'r+'); xlabel('\Sigma_{xxp}'); ylabel('\Sigma_{xpxp}');
subplot(234); hold on; plot(ayy,aypyp,'r+'); xlabel('\Sigma_{yy}'); ylabel('\Sigma_{ypyp}');
subplot(235); hold on; plot(ayy,ayyp,'r+'); xlabel('\Sigma_{yy}'); ylabel('\Sigma_{yyp}');
subplot(236); hold on; plot(ayyp,aypyp,'r+'); xlabel('\Sigma_{yyp}'); ylabel('\Sigma_{ypyp}');