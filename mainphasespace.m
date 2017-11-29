
%%This script selects particular files to analyse, loads an image, and sends it to be evalauated 
%%by the phase space routine. 

%% Initialization and Parameters
clc;
clear all;
%close all;
warning off;

addpath('./Analysis_Tools/')

%% Import parameters and files

Input_data;

results=[]; Info=[];

    filename = [folder file1];
    display(file1);

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
    Acrop= A(center(2)-ry:center(2)+ry,center(1)-rx:center(1)+rx);
    
    for i=1:size(Acrop,1)
        for j=1:size(Acrop,2)
            
            if (Xcrop(i,j)-center(1))^2/rx^2+(Ycrop(i,j)-center(2))^2/ry^2>1
                Acrop(i,j)=mean(BGel);
            end
        end
    end
    figure(3); clf; imagesc(Acrop);
    
    A=Acrop-mean(BGel); 
end
A(A<0)=0;

%%
xvec = pxconv*(1:size(A,2));
yvec = pxconv*(1:size(A,1));

% keep = input('*Analyse this image? (0 or 1): '); %Ask user if they like the image they see.
keep = true;

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
    
    phasespace_shear;
    phasespace_shear_interpolation;
    
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
    
    %         savedata= input('*Save results? (0 or 1): ');
    savedata = true;
    if savedata
        aux.file=file1;
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
        results= [results aux];
        Info= [Info info];
    end
    
    display_phase_space;
    
end


