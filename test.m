clc
clear
close all
load('result.mat')
load('wine.mat')

for i=1:10
    purity1=purity(label,QAP_labels{1,i});
    X_QAP(i,4)=purity1;
    purity2=purity(label,kmean_labels{1,i});
    kmean(i,4)=purity2;
    purity3=purity(label,kmeanpro_labels{1,i});
    kmeanspro(i,4)=purity3;
    purity4=purity(label,isodata_labels{1,i});
    isodata(i,4)=purity4;
end
x_QAP=X_QAP;
purity_QAP_avg=mean(x_QAP(:,4)*100);
purity_QAP_max=max(x_QAP(:,4)*100);
purity_QAP_min=min(x_QAP(:,4)*100);
purity_QAP_std=std(x_QAP(:,4)*100);

purity_kmean_avg=mean(kmean(:,4)*100);
purity_kmean_max=max(kmean(:,4)*100);
purity_kmean_min=min(kmean(:,4)*100);
purity_kmean_std=std(kmean(:,4)*100);

purity_kmeanpro_avg=mean(kmeanspro(:,4)*100);
purity_kmeanpro_max=max(kmeanspro(:,4)*100);
purity_kmeanpro_min=min(kmeanspro(:,4)*100);
purity_kmeanpro_std=std(kmeanspro(:,4)*100);

purity_isodata_avg=mean(isodata(:,4)*100);
purity_isodata_max=max(isodata(:,4)*100);
purity_isodata_min=min(isodata(:,4)*100);
purity_isodata_std=std(isodata(:,4)*100);
