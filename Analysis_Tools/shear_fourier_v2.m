function [S1,S2] = shear_fourier_v2(rhon)
%SHEAR_FOURIER_V2 Use 2D Fourier transform to find shear factors of image

    Afft=abs(fftshift(fft2(rhon))); %Afft = Afft/max(Afft(:));
    figure(91); clf; set(gcf,'Name','2D Fourier transform')
    imagesc(Afft); colormap('parula'); colorbar;
    hold on;
    
    mask = ones(3); mask(5) = 0; % 3x3 mask
    B=imdilate(Afft,mask); 
    mypeaks= Afft>B;
    
    idx=find(mypeaks);
    valpeaks=Afft(idx);
    [valpeaks,ii]=sort(valpeaks,'descend'); idx=idx(ii);
    filt= find(valpeaks>0.1*max(valpeaks));
    
    [row,col]=ind2sub(size(Afft),idx(filt));  plot(col,row,'r.');
    
    rowcen=round(size(Afft,1)/2);
    colcen=round(size(Afft,2)/2);
    
    [~,hor]=sort(abs(row-rowcen)); hor=hor(1:min(7,length(hor))); horx=col(hor); hory=row(hor);
    pp=polyfit(horx,hory,1); plot([min(horx) max(horx)],polyval(pp,[min(horx) max(horx)]),'y-'); S1=pp(1);

    [~,ver]=sort(abs(col-colcen)); ver=ver(1:min(7,length(ver))); verx=col(ver); very=row(ver); 
    pp=polyfit(very,verx,1); plot(polyval(pp,[min(very) max(very)]),[min(very) max(very)],'y-'); S2=pp(1);
    drawnow;


end
