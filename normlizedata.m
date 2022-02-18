function fea=normlizedata(fea,op)
if op==1
    [nSmp,nFea] = size(fea);
    for i = 1:nSmp
        fea(i,:) = fea(i,:) ./ max(1e-12,norm(fea(i,:)));
    end
else
    maxValue = max(max(fea));
    fea = fea/maxValue;
end