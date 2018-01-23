function [Int_list,X_4D,Xp_4D,Y_4D,Yp_4D,r0,Scov] = Density_4D_discrete(x,xp,y,yp,sxp,syp,sxpyp,int, Npoints, lb, ub)
% Npoints: Number of points to resolve x' and y'
% lb, ub: Lower and upper bounds of the region of interest for x' and y' 
nb_sigmas = 3;

rho=sxpyp./sxp./syp; 
int(rho.^2>=1)=0; % discard points with bad interpolation of rho

roi=find(int>0);
x=x(roi);
xp=xp(roi);
y=y(roi);
yp=yp(roi);
sxp=sxp(roi);
syp=syp(roi);
sxpyp=sxpyp(roi); 
rho=rho(roi);
int=int(roi);

if nargin<=8 %default limits and number of points
    Npoints = 200*[1 1];
    lb= [min(xp) min(yp)];
    ub= [max(xp) max(yp)];
elseif nargin<=9
    lb= [min(xp) min(yp)];
    ub= [max(xp) max(yp)];
end

for i=1:length(Npoints) %limit number of points to 100
    if Npoints(i)>300
        Npoints(i)=300; display(['Number of points of variable ' num2str(i) ' is too large. Set to maximum 300.'])
    end
end


xp_meshsize = (ub(1)-lb(1))/Npoints(1) ; 
yp_meshsize = (ub(2)-lb(2))/Npoints(2);
min_nb_points=10;
max_nb_points=100;


progress=5;
roi_size = numel(x);

R_list=cell(roi_size,1);
Int_list=cell(roi_size,1);

for i=1:roi_size
    
    npoints_mesh=2*nb_sigmas*sxp(i)/xp_meshsize; 
    if npoints_mesh>max_nb_points; npoints_mesh=max_nb_points; elseif npoints_mesh<min_nb_points; npoints_mesh=min_nb_points; end;
    xp_mesh= linspace(xp(i)-nb_sigmas*sxp(i),xp(i)+nb_sigmas*sxp(i),npoints_mesh);
    npoints_mesh=2*nb_sigmas*syp(i)/yp_meshsize; 
    if npoints_mesh>max_nb_points; npoints_mesh=max_nb_points; elseif npoints_mesh<min_nb_points; npoints_mesh=min_nb_points; end;
    yp_mesh= linspace(yp(i)-nb_sigmas*syp(i),yp(i)+nb_sigmas*syp(i),npoints_mesh);
    [xp_mesh,yp_mesh]=meshgrid(xp_mesh,yp_mesh);
    x_mesh=x(i)*ones(size(xp_mesh));
    y_mesh=y(i)*ones(size(xp_mesh));

    R_list{i}=[x_mesh(:) xp_mesh(:) y_mesh(:) yp_mesh(:)];
    Int_list{i}=int(i)./(2*pi.*sxp(i).*syp(i).*sqrt(1-rho(i).^2)).* ...
                    exp(-1./(2*(1-rho(i).^2)).*( ...
                        (xp_mesh(:)-xp(i)).^2./(sxp(i).^2)+(yp_mesh(:)-yp(i)).^2./(syp(i).^2) ...
                        -2*rho(i).*(xp_mesh(:)-xp(i)).*(yp_mesh(:)-yp(i))./(sxp(i).*syp(i)) ...
                        ) );
    
    if i/roi_size*100>=progress
        display(['4D density: ' num2str(progress) '%'])
        progress=progress+5;
    end
end

Int_list=cell2mat(Int_list);
R_list=cell2mat(R_list);
X_4D = R_list(:,1);
Xp_4D = R_list(:,2);
Y_4D = R_list(:,3);
Yp_4D = R_list(:,4);

Int_list=Int_list/sum(Int_list(:));
r0 = Int_list'*R_list;
Scov = weightedcov(R_list,Int_list);


end