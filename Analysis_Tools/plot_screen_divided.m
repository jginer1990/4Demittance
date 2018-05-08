function [] = plot_screen_divided(A,Locsx,Locsy)
%PLOT_SCREEN_DIVIDED Plots image A with lines denoting Locsx and Locsy

figure(99); clf
imagesc(A); 
set(gcf,'Name','Screen image')
set(gca,'YDir','normal')
hold on
for i=1:size(Locsx,1)
    plot([Locsx(i,1) Locsx(i,end)],[Locsy(i,1) Locsy(i,end)],'y:')
end
for i=1:size(Locsx,2)
    plot([Locsx(1,i) Locsx(end,i)],[Locsy(1,i) Locsy(end,i)],'y:')
end

end

