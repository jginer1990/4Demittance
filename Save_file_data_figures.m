%% Save data
fname=files(ifile).name(1:(end-4));
file_save=[folder 'Data_' files(ifile).name(1:(end-4)) '.mat'];

save(file_save, 'fname'); % Overwrites if already exists

saveparameters.pxconv = pxconv;
saveparameters.gamma = gamma;
saveparameters.beta = beta;
saveparameters.gridSpacing = mask_prop.gridSpacing;
saveparameters.pitch_to_bar_width_ratio = mask_prop.pitch_to_bar_width_ratio;
saveparameters.driftLength = mask_prop.driftLength;
save(file_save, 'saveparameters', '-append');

res.S=S;
res.exn=ex*gamma*beta;
res.eyn=ey*gamma*beta;
res.e1n=e1*gamma*beta;
res.e2n=e2*gamma*beta;
res.S_interp=S_interp;
res.exn_interp=ex_interp*gamma*beta;
res.eyn_interp=ey_interp*gamma*beta;
res.e1n_interp=e1_interp*gamma*beta;
res.e2n_interp=e2_interp*gamma*beta;
if compute_core
    res.exn2Dcore=ex2Dcore*gamma*beta;
    res.eyn2Dcore=ey2Dcore*gamma*beta;
    res.en4Dcore=e4Dcore*(gamma*beta)^2;
end
res.counts = sum(A(:));
save(file_save, 'res', '-append');

if strcmp(target,'TEM')
    analysis.rhon= A;
    analysis.xb= xb;
    analysis.ybc= ybc;
    analysis.sigmaxp= sigmaxp;
    analysis.intx= intx;
    analysis.ybc= ybc;
    analysis.yb= yb;
    analysis.xbc= xbc;
    analysis.sigmayp= sigmayp;
    analysis.inty= inty;
    analysis.xg= xg;
    analysis.yg= yg;
    analysis.xp= xp;
    analysis.yp= yp;
    analysis.xcen= info.xcen;
    analysis.ycen= info.ycen;
    analysis.xpcen= info.xpcen;
    analysis.ypcen= info.ypcen;
    analysis.intcen= info.intcen;
elseif strcmp(target,'PP')
    analysis.rhon= A;
    analysis.xb= info.xb;
    analysis.sigmaxp= info.sigmaxp;
    analysis.yb= info.yb;
    analysis.sigmayp= info.sigmayp;
    analysis.sigmaxpyp= info.sigmaxpyp;
    analysis.int= info.int;
    analysis.xg= info.xg;
    analysis.yg= info.yg;
    analysis.xp= info.xp;
    analysis.yp= info.yp;
end
save(file_save,'analysis','-append');

%% Save figures in .fig and .png
file_save=[folder 'AnalysisLocs_' files(ifile).name(1:(end-4))];
saveas(99,[file_save '.png'])
saveas(99,[file_save '.fig'])

file_save=[folder 'AnalysisInt_' files(ifile).name(1:(end-4))];
saveas(199,[file_save '.png'])
saveas(199,[file_save '.fig'])

if compute_core
    file_save=[folder 'Proj2Dphasespace_' files(ifile).name(1:(end-4))];
    saveas(780,[file_save '.png'])
    saveas(780,[file_save '.fig'])
    
    file_save=[folder 'Core2D_' files(ifile).name(1:(end-4))];
    saveas(890,[file_save '.png'])
    saveas(890,[file_save '.fig'])
    
    file_save=[folder 'Core4D_' files(ifile).name(1:(end-4))];
    saveas(891,[file_save '.png'])
    saveas(891,[file_save '.fig'])
end