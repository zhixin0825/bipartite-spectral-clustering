function [ ac ] = compute_acc( y,z )
%ACC Summary of this function goes here
%   Detailed explanation goes here
K = size(y,2);
L = size(z,2);
ac=0;
if K==L
    n = size(y,1);
    R = y'*z./n;
    for k = 1:K
        [~,ind]=max(R);
        [v2,ind1]=max(max(R));
        ac = ac+v2;
        if k<K
            i = ind(ind1);
            j = ind1;
            R(i,:)=[];
            R(:,j)=[];
        end
    end
end
end

