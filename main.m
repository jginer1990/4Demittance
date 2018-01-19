
%%This script selects particular files to analyse, loads an image, and sends it to be evaluated 
%%by the phase space routine. 

%% Initialization
clear
close all
% warning off;
% dbstop if all error

addpath('./Analysis_Tools/')
addpath('./Core emittance/')

%% Import parameters and files

Input_data;

searchstr=[];
savestr=[];
for ii=1:length(keywords)
   searchstr = [searchstr '*' keywords{ii}]; 
   savestr = [savestr, '_' keywords{ii}]; 
end

files = dir([folder searchstr]);

%% Routine

for ifile=1:length(files)
    filename = [folder files(ifile).name];
    display([num2str(ifile) '/' num2str(length(files)) ': ' files(ifile).name])
    
    A = imread(filename); % Read file
    A = double(A); %convert int to double.
    figure(1); clf
    imagesc(A);
    BG = mean(mean(A(1:30,1:30))); A=A-BG;
    
    [A,Acrop0] = crop(A,center,rx,ry);

    xvec = pxconv*(1:size(A,2));
    yvec = pxconv*(1:size(A,1));
    
    [X, Y] = meshgrid(xvec, yvec); %build screen coordinate matrices.
    X0 = trapz(trapz(X.*A))/trapz(trapz(A)); %Mean value of x. A is 2D so need trapz(trapz()) for 2D integral. Assume symmetric distribution.
    Y0 = trapz(trapz(Y.*A))/trapz(trapz(A)); %Mean value of y.
    X = X-X0; %center coords.
    Y = Y-Y0;
    xvec = xvec-X0;
    yvec = yvec-Y0;
    
    keep = true;%input('*Analyse this image? (0 or 1): '); %Ask user if they like the image they see.
    compute_core = false; % include computation of 2D/4D core emittance

    if keep
        [Xp,Yp]=meshgrid(1:size(A,2),1:size(A,1));
        X0px=sum(A(:).*Xp(:))/sum(A(:));
        Y0px=sum(A(:).*Yp(:))/sum(A(:));

        try
            if strcmp(target,'TEM')
                [S1,S2] = shear_fourier_v1(A);
            elseif strcmp(target,'PP')
                [S1,S2] = shear_fourier_v2(A);
            end
        catch
            [S1,S2] = shear_manual(A);
        end

        Asheared = shear_image(Xp,Yp,A,S1,S2,X0px,Y0px);

        if strcmp(target,'TEM')
            [locsx_sh,locsy_sh,minsx,minsy] = splitimage_TEM(Asheared,analysis);
            [Locsx_sh,Locsy_sh] = meshgrid(locsx_sh,locsy_sh);
            [Locsx,Locsy] = undo_shear(Locsx_sh,Locsy_sh,X0px,Y0px,S1,S2);
            plot_screen_divided(A,Locsx,Locsy);
            [xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp] = phasespace_TEM(A,X,Y,locsx_sh,locsy_sh,Locsx,Locsy,mask_prop,analysis);
            [S,info] = beammatrix_TEM(xg,yg,xb,yb,xs,ys,xbcen,ybcen,xbc,ybc,xp,yp,intx,inty,sigmaxp,sigmayp);
            [S_interp,info_interp] = interpolation_TEM(Xp,Yp,info);
            if compute_core
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
            end

        elseif strcmp(target,'PP')
            [locsx_sh,locsy_sh] = splitimage_PP(Asheared,analysis);
            [Locsx_sh,Locsy_sh] = meshgrid(locsx_sh,locsy_sh);
            [Locsx,Locsy] = undo_shear(Locsx_sh,Locsy_sh,X0px,Y0px,S1,S2);
            plot_screen_divided(A,Locsx,Locsy);
            [S,info] = phasespace_PP(A,Xp,Yp,X,Y,X0px,Y0px,S1,S2,mask_prop,locsx_sh,locsy_sh);
            [S_interp,info_interp] = interpolation_PP(Xp,Yp,info);
            if compute_core
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
        end
        
        disp(length(locsx_sh))
        disp(length(locsy_sh))
        
        [ex,ey,e1,e2]=Emittance_2D_4D(S);
        [ex_interp,ey_interp,e1_interp,e2_interp]=Emittance_2D_4D(S_interp);

        %% sort e1 & e2
        sorted= sort([e1 e2]);
        if ex<=ey
            e1=sorted(1); e2=sorted(2);
        else
            e1=sorted(2); e2=sorted(1);
        end
        sorted= sort([e1_interp e2_interp]);
        if ex_interp<=ey_interp
            e1_interp=sorted(1); e2_interp=sorted(2);
        else
            e1_interp=sorted(2); e2_interp=sorted(1);
        end

        %% Save results
        savedata=1; %input('*Save results? (0 or 1): ');
        if savedata
            aux.file=files(ifile).name;
            aux.S=S;
            aux.ex=ex*gamma*beta;
            aux.ey=ey*gamma*beta;
            aux.e1=e1*gamma*beta;
            aux.e2=e2*gamma*beta;
            aux.S_interp=S_interp;
            aux.ex_interp=ex_interp*gamma*beta;
            aux.ey_interp=ey_interp*gamma*beta;
            aux.e1_interp=e1_interp*gamma*beta;
            aux.e2_interp=e2_interp*gamma*beta;
            aux.charge = sum(Acrop0(:));
            if compute_core
                aux.ex2Dcore=ex2Dcore*gamma*beta;
                aux.ey2Dcore=ey2Dcore*gamma*beta;
                aux.e4Dcore=e4Dcore*(gamma*beta)^2;
            else
                aux.ex2Dcore=NaN;
                aux.ey2Dcore=NaN;
                aux.e4Dcore=NaN;
            end
            results(ifile)= aux;
            Info(ifile)= info;
        else
            aux.file=files(ifile).name;
            aux.S=NaN(4); aux.ex=NaN; aux.ey=NaN; aux.e1=NaN; aux.e2=NaN;
            aux.S_interp=NaN(4); aux.ex_interp=NaN; aux.ey_interp=NaN; aux.e1_interp=NaN; aux.e2_interp=NaN;
            aux.charge=NaN; aux.ex2Dcore=NaN; aux.ey2Dcore=NaN; aux.e4Dcore=NaN;
            results(ifile)= aux;
            Info(ifile)= info;
        end

        
    end
    
end

analysis_emittance; % Plot summary of shot-to-shot emittance results 
