function em = emittance_fraction(r,Int)
% Caution: Int(all set) must be normalized sum(Int(:))=1

% N=sum(Int); %Normalization only subset
% 
% ndim = size(r,2);
% S=zeros(ndim);
% 
% if ndim>1
%     pairs = nchoosek(1:ndim,2);
%     pairs = [[(1:ndim)' (1:ndim)']; pairs];
% else
%     pairs = [1 1];
% end
%     
% 
% mu=[0 0];
% for i=1:size(pairs,1)
%     mu(1)=sum(Int.*r(:,pairs(i,1)))/N;
%     mu(2)=sum(Int.*r(:,pairs(i,2)))/N;
%     S(pairs(i,1),pairs(i,2))=sum(Int.*(r(:,pairs(i,1))-mu(1)).*(r(:,pairs(i,2))-mu(2))  )/N;
%     S(pairs(i,2),pairs(i,1))=S(pairs(i,1),pairs(i,2));
% end

S = weightedcov(r,Int);

em = sqrt(abs(det(S)));

end