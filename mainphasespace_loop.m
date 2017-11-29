
%%This script selects particular files to analyse, loads an image, and sends it to be evalauated 
%%by the phase space routine. 

%% Initialization and Parameters
% clc;
clear all;
%close all;
warning off;

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

%%
clear results Info
for ifile=1:length(files)
    filename = [folder files(ifile).name];
    display([num2str(ifile) '/' num2str(length(files)) ': ' files(ifile).name])
    
    A = imread(filename); % Read file
    A = double(A); %convert int to double.
    
    %% pre-processing image : crop image inside the glowing screen edge
    figure(1); clf
    imagesc(A);
    BG= mean(mean(A(1:30,1:30))); A=A-BG;
    
    if true
        t = linspace(0,2*pi,100);
        xel=center(1)+rx*cos(t);
        yel=center(2)+ry*sin(t);
        for i= 1:length(xel)
            BGel(i)=A(round(yel(i)),round(xel(i)));
        end
        
        hold on; plot(xel,yel,'y-')
        Xcrop = 1:size(A,2); Ycrop=1:size(A,1);
        [Xcrop,Ycrop]=meshgrid(Xcrop,Ycrop);
        Xcrop=Xcrop(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx);
        Ycrop=Ycrop(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx);
        Acrop= A(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx); Acrop2=Acrop;
        for i=1:size(Acrop,1)
            for j=1:size(Acrop,2)
                
                if (Xcrop(i,j)-center(1))^2/rx^2+(Ycrop(i,j)-center(2))^2/ry^2>1
                    Acrop(i,j)=mean(BGel);
                    Acrop2(i,j)=0;
                end
            end
        end
                
        A=Acrop-mean(BGel);
        A(A<0)=0; 
%         A=A';
%         A=imrotate(A,90);
        figure(3); clf; imagesc(A); set(gca,'YDir','normal')
    end
    
    %%
    xvec = pxconv*(1:size(A,2));
    yvec = pxconv*(1:size(A,1));
    
    keep = true;%input('*Analyse this image? (0 or 1): '); %Ask user if they like the image they see.
    % keep = 1;
    
    % X = zeros(size(yvec,2),size(xvec,2));
    % Y = zeros(size(yvec,2),size(xvec,2));
    [X, Y] = meshgrid(xvec, yvec); %build screen coordinate matrices .
    At = A; %Create a thresholded version of A for centroid determination.
    X0 = trapz(trapz(X.*At))/trapz(trapz(At)); %Mean value of x. At is 2D so need trapz(trapz()) for 2D integral. Assume symmetric distribution.
    Y0 = trapz(trapz(Y.*At))/trapz(trapz(At)); %Mean value of y.
    X = X-X0; %center coords.
    Y = Y-Y0;
    xvec = xvec-X0;
    yvec = yvec-Y0;
    
    
    
    if keep
        
%         try
            if strcmp(target,'TEM')
                phasespace_shear;
                phasespace_shear_interpolation;
            elseif strcmp(target,'PP')
                phasespace_shear_PP;
                phasespace_shear_interpolation_PP;
            end
            
            [ex,ey,e1,e2]=Emittance_2D_4D(S);
            [ex_interp,ey_interp,e1_interp,e2_interp]=Emittance_2D_4D(S_interp);
            
            %sort e1 & e2
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
                aux.charge = sum(Acrop2(:));
                results(ifile)= aux;
                Info(ifile)= info;
            else
                aux.file=files(ifile).name;
                aux.S=NaN(4); aux.ex=NaN; aux.ey=NaN; aux.e1=NaN; aux.e2=NaN;
                aux.S_interp=NaN(4); aux.ex_interp=NaN; aux.ey_interp=NaN; aux.e1_interp=NaN; aux.e2_interp=NaN;
                results(ifile)= aux;
                Info(ifile)= info;
            end
%         catch
%             display('!Error in routine.')
%             aux.file=files(ifile).name;
%             aux.S=NaN(4); aux.ex=NaN; aux.ey=NaN; aux.e1=NaN; aux.e2=NaN;
%             aux.S_interp=NaN(4); aux.ex_interp=NaN; aux.ey_interp=NaN; aux.e1_interp=NaN; aux.e2_interp=NaN;
%             results(ifile)= aux;
%             Info(ifile).target='Error';
%         end
        
    end
    
end

analysis_emittance; % Plot summary of shot-to-shot emittance results 