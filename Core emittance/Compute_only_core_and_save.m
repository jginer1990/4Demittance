%% Compute ONLY core emittances

searchmatfile = ['Data_' keywords{1} '*.mat'];

filesdata = dir([folder searchmatfile]);

for i=1:1%length(filesdata)
    mydata=load([folder filesdata(i).name]);
    display([num2str(i) '/' num2str(length(filesdata))    mydata.fname '.tif']);
    %     A=load([folder mydata.fname '.mat'],'image'); A=A.image;
    A = imread([folder mydata.fname '.tif']); % Read file
    A = double(A); %convert int to double.
    BG = mean(mean(A(1:30,1:30))); A=A-BG;
    %     [A,Acrop0] = crop_rotate(A,center,rx,ry,angledeg);
    [A,Acrop0] = crop(A,center,rx,ry);
    [Xp,Yp]=meshgrid(1:size(A,2),1:size(A,1));
    
    
    driftLength = mydata.saveparameters.driftLength;
    gamma = mydata.saveparameters.gamma;
    beta = mydata.saveparameters.beta;
    
    res = mydata.res;
    
    try % for TEM
        xb= mydata.analysis.xb;
        yb= mydata.analysis.yb;
        ybc = mydata.analysis.ybc;
        xbc = mydata.analysis.xbc;
        xg = mydata.analysis.xg;
        yg = mydata.analysis.yg;
        xp = mydata.analysis.xp;
        yp = mydata.analysis.yp;
        intx = mydata.analysis.intx;
        inty = mydata.analysis.inty;
        sigmaxp = mydata.analysis.sigmaxp;
        sigmayp = mydata.analysis.sigmayp;
        
        
        clear ybmid ybcen xbmid xbcen
        ybmid(:,2:size(yb,2)+1) = yb(:,:);
        ybmid(:,1) = ybmid(:,2);
        ybmid(:,end+1) = ybmid(:,end);
        ybcen = (ybmid(:,1:end-1)+ybmid(:,2:end)) /2;
        
        xbmid(:,2:size(xb,1)+1) = xb(:,:)';
        xbmid(:,1) = xbmid(:,2);
        xbmid(:,end+1) = xbmid(:,end);
        xbcen = ((xbmid(:,1:end-1)+xbmid(:,2:end))/2)';
        
        xs = xg + driftLength*xp;
        ys = yg + driftLength*yp;
        
        
        [S,info] = beammatrix_TEM(xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp);
        [S_interp,info_interp] = interpolation_TEM(Xp,Yp,info);
        
        [ex2Dcore,ey2Dcore] = core2Demittance(...
            info_interp.xg_interp,info_interp.xp_interp,...
            info_interp.sigmaxp_interp,info_interp.intx_interp,...
            info_interp.yg_interp,info_interp.yp_interp,...
            info_interp.sigmayp_interp,info_interp.inty_interp);
        e4Dcore = core4Demittance(...
            info_interp.xg_interp,info_interp.xp_interp,...
            info_interp.yg_interp,info_interp.yp_interp,...
            info_interp.sigmaxp_interp,info_interp.sigmayp_interp,...
            info_interp.sigmaxpyp_interp,sqrt(info_interp.intx_interp.*info_interp.inty_interp));
        
    catch % for PP
        
        info.target='PP';
        info.xb=mydata.analysis.xb;
        info.yb=mydata.analysis.yb;
        info.sigmaxp=mydata.analysis.sigmaxp;
        info.sigmayp=mydata.analysis.sigmayp;
        info.sigmaxpyp=mydata.analysis.sigmaxpyp;
        info.xg=mydata.analysis.xg;
        info.yg=mydata.analysis.yg;
        info.xs=mydata.analysis.xg+driftLength*mydata.analysis.xp;
        info.ys=mydata.analysis.yg+driftLength*mydata.analysis.yp;
        info.xp=mydata.analysis.xp;
        info.yp=mydata.analysis.yp;
        info.int=mydata.analysis.int;
        info.x0=sum(info.int(:).*info.xg(:))/sum(info.int(:));
        info.y0=sum(info.int(:).*info.yg(:))/sum(info.int(:));
        info.xp0=sum(info.int(:).*info.xp(:))/sum(info.int(:));
        info.yp0=sum(info.int(:).*info.yp(:))/sum(info.int(:));
        info.S=mydata.res.S;
        
        [S_interp,info_interp] = interpolation_PP(Xp,Yp,info);
        
        [ex2Dcore,ey2Dcore] = core2Demittance(...
            info_interp.xg_interp,info_interp.xp_interp,...
            info_interp.sigmaxp_interp,info_interp.int_interp,...
            info_interp.yg_interp,info_interp.yp_interp,...
            info_interp.sigmayp_interp,info_interp.int_interp);
        e4Dcore = core4Demittance(...
            info_interp.xg_interp,info_interp.xp_interp,...
            info_interp.yg_interp,info_interp.yp_interp,...
            info_interp.sigmaxp_interp,info_interp.sigmayp_interp,...
            info_interp.sigmaxpyp_interp,info_interp.int_interp);

    end
    
    res.exn2Dcore=ex2Dcore*gamma*beta;
    res.eyn2Dcore=ey2Dcore*gamma*beta;
    res.en4Dcore=e4Dcore*(gamma*beta)^2;
    
%     save([folder filesdata(i).name], 'res', '-append');
end