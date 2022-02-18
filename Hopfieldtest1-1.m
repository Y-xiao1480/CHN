function [x]=Hopfieldtest(x,n,d,p,alpha,u0,u)
B=alpha; 
%u=zeros(n,p);
for it=1:2000
    y=x;
    rd=randperm(n*p);
    for index=1:n*p
        [i,k]=ind2sub(size(x),rd(index));
%         dedx=0.5*d(i,:)*x(:,k)+B*(x(i,k)+sum(x(i,:))-x(i,k)-2+p);
        dedx=0.5*d(i,:)*x(:,k)+B*(sum(x(i,:))-x(i,k)-0.5);
        u(i,k)=u(i,k)-dedx*0.001; %iris:0.00002 8000  seed£º0.00002 8000  wine:0.00002 8000   face:10000 0.00001 cosin
%         y(i,k)=tanh(u(i,k)/u0);
        y(i,k)=sigmoid(u(i,k)/u0);
    end
    x=y;
end

cc=1;
%{
    for k =1:p
    for i=1:n
        for j=1:n
            E=(d(i,j)*x(i,k)*x(j,k))+alpha*(sum(x(i,:))-2+p)*(sum(x(i,:))-2+p);
            E=E+E;
        end
    end
    end
    times(it)=it;
    Ep(it)=E;
%}