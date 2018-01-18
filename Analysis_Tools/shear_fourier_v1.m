function [S1,S2] = shear_fourier_v1(rhon)
%SHEAR_FOURIER_V1 Use 2D Fourier transform to find shear factors of image

Afft=abs(fftshift(fft2(rhon))); %Afft = Afft/max(Afft(:));
figure(91); clf; set(gcf,'Name','2D Fourier transform')
imagesc(Afft); colormap('parula'); colorbar;
hold on;

mask = ones(3); mask(5) = 0; % 3x3 mask
%     mask = ones(7); mask(4,4) = 0; % 5x5 mask
B=imdilate(Afft,mask); 
mypeaks= Afft>B;

idx=find(mypeaks);
valpeaks=Afft(idx);
[valpeaks,ii]=sort(valpeaks,'descend'); idx=idx(ii);
filt= find(valpeaks>0.1*max(valpeaks));

[row,col]=ind2sub(size(Afft),idx(filt));

rowcen=round(size(Afft,1)/2);
colcen=round(size(Afft,2)/2);

aux=abs(row-rowcen); [~,hor]=sort(aux); hor=hor(2:3);
aux=abs(col-colcen); [~,ver]=sort(aux); ver=ver(2:3);

plot(col,row,'r.');
plot(col(hor),row(hor),'g.'); plot(col(ver),row(ver),'y.'); drawnow;


S2 = (col(ver(1))-col(ver(2)))/(row(ver(1))-row(ver(2)));
S1 = (row(hor(1))-row(hor(2)))/(col(hor(1))-col(hor(2)));

end
