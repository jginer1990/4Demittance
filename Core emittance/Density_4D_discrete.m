function [Int_list,R_list,r0,Scov] = Density_4D_discrete(x,xp,y,yp,sxp,syp,sxpyp,int)

rho=sxpyp./sxp./syp;

nb_sigmas = 3;
npoints_mesh=100;

Int_list=[];
R_list=[];
for i=1:numel(x)
    
    xp_mesh= linspace(xp(i)-nb_sigmas*sxp(i),xp(i)+nb_sigmas*sxp(i),npoints_mesh);
    yp_mesh= linspace(yp(i)-nb_sigmas*syp(i),yp(i)+nb_sigmas*syp(i),npoints_mesh);
    [xp_mesh,yp_mesh]=meshgrid(xp_mesh,yp_mesh);
    x_mesh=x(i)*ones(size(xp_mesh));
    y_mesh=y(i)*ones(size(xp_mesh));
    int_mesh= 1./(2*pi.*sxp(i).*syp(i).*sqrt(1-rho(i).^2)).* ...
                    exp(-1./(2*(1-rho(i).^2)).*( ...
                        (xp_mesh-xp(i)).^2./(sxp(i).^2)+(yp_mesh-yp(i)).^2./(syp(i).^2) ...
                        -2*rho(i).*(xp_mesh-xp(i)).*(yp_mesh-yp(i))./(sxp(i).*syp(i)) ...
                        ) );
    
    R_list=[R_list; [x_mesh(:) xp_mesh(:) y_mesh(:) yp_mesh(:)]];
    Int_list=[Int_list; int_mesh];
end

X=R_list(:,1);
Xp=R_list(:,2);
Y=R_list(:,3);
Yp=R_list(:,4);

Int=Int/sum(Int(:));
r0=[sum(Int(:).*X(:)) , ...
    sum(Int(:).*Xp(:)) , ...
    sum(Int(:).*Y(:)) , ...
    sum(Int(:).*Yp(:)) , ...
    ];
Scov(1,1)=sum(Int(:).*(X(:)-r0(1)).^2) ;
Scov(2,2)=sum(Int(:).*(Xp(:)-r0(2)).^2) ;
Scov(3,3)=sum(Int(:).*(Y(:)-r0(3)).^2) ;
Scov(4,4)=sum(Int(:).*(Yp(:)-r0(4)).^2) ;
Scov(1,2)=sum(Int(:).*(X(:)-r0(1)).*(Xp(:)-r0(2))) ; Scov(2,1)=Scov(1,2);
Scov(1,3)=sum(Int(:).*(X(:)-r0(1)).*(Y(:)-r0(3))) ; Scov(3,1)=Scov(1,3);
Scov(1,4)=sum(Int(:).*(X(:)-r0(1)).*(Yp(:)-r0(4))) ; Scov(4,1)=Scov(1,4);
Scov(2,3)=sum(Int(:).*(Xp(:)-r0(2)).*(Y(:)-r0(3))) ; Scov(3,2)=Scov(2,3);
Scov(2,4)=sum(Int(:).*(Xp(:)-r0(2)).*(Yp(:)-r0(4))) ; Scov(4,2)=Scov(2,4);
Scov(3,4)=sum(Int(:).*(Y(:)-r0(3)).*(Yp(:)-r0(4))) ; Scov(4,3)=Scov(3,4);



end