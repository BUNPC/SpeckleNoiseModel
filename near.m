function [n,i] = near(x, v,varargin)
%NEAR Find the nearest value in an ordered unique array x to the value v
k=1;
n=zeros(1,length(v));
i=n;
for k=1:length(v)
    if (size(x,2)>1)&(size(x,1)>1)
        xp=reshape(x,1,numel(x));
        [n(1,k),i(1,k)]=near(xp,v(k));
    else
        distancias=abs(x-v(k));
        [i(1,k),i(1,k)]=min(distancias);
        n(1,k)=x(i(k));
        
    end
end
n=n';
i=i';

end

