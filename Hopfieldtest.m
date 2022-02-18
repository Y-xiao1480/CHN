function [x]=Hopfieldtest(x,n,d,p,alpha,u0,u)
B=alpha; 
%u=zeros(n,p);
for it=1:15000
    y=x;
    rd=randperm(n*p);
    for index=1:n*p
        [i,k]=ind2sub(size(x),rd(index));
        dedx=0.5*d(i,:)*x(:,k)+B*(x(i,k)+sum(x(i,:))-x(i,k)-2+p);
        u(i,k)=u(i,k)-dedx*0.0025; 
        y(i,k)=tanh(u(i,k)/u0);
    end
    x=y;
end
