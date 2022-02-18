clc
clear
close all
load('wine.mat');
if ~exist('label','var')
    label=gnd;
end
p=3;
u0=0.025;
nClass=max(unique(label));
fea=normlizedata(fea,2);
fea=zscore(fea);
[n,m]=size(fea);
options.ReducedDim=fix(0.9*m);
W = PCA(fea,options);
[n,~]=size(fea);
fea=fea*W;

d=gaussinKernel(fea,1.5); %1.5
d=1-d;
d=d-diag(diag(d));

kmean=zeros(10,2);
kmeanspro=zeros(10,2);
isodata=zeros(10,2);
X_QAP=zeros(10,2);
for i =1 :10
u=2*rand(n,p)-1;
U=u0*log(n-1)+u;
s=tanh(U/u0);
s(s>=0)=1;
s(s<0)=-1;
[s]=Hopfieldtest(s,n,d,p,700,u0,U); 
s(s>=0)=1;
s(s<0)=-1;
sum(s,2)
l1=label1(s,n);
[NMI1,AC1]=ACNMI(l1,label);
[ARI1]=RandIndex(label,l1);
X_QAP(i,1)=AC1;
X_QAP(i,2)=NMI1;
X_QAP(i,3)=ARI1;
QAP_labels{i}=l1;
l2=kmeans(d,p);
[NMI2,AC2]=ACNMI(l2,label);
[ARI2]=RandIndex(label,l2);
kmean(i,1)=AC2;
kmean(i,2)=NMI2;
kmean(i,3)=ARI2;
kmean_labels{i}=l2;
[centroid3, l3]=kmeans_pro(d,'kmeans++',p,50);
[ARI3]=RandIndex(label,l3);
[NMI3,AC3]=ACNMI(l3,label);
kmeanspro(i,1)=AC3;
kmeanspro(i,2)=NMI3;
kmeanspro(i,3)=ARI3;
kmeanpro_labels{i}=l3;
[centroid4, l4]=IsoData(d,'ISODATA',p,100,5,5,0.5);
[NMI4,AC4]=ACNMI(l4,label);
[ARI4]=RandIndex(label,l4);
isodata(i,1)=AC4;
isodata(i,2)=NMI4;
isodata(i,3)=ARI4;
isodata_labels{i}=l4;
end




