%% Batch Correction

load ChronicWoundData.mat %result of ReadData.m

data=log2(1+data);
for i=1:numel(Gene)
    tmp1=mean(data(i,batch_label==1));
    tmp2=mean(data(i,batch_label==2));
    tsd1= std(data(i,batch_label==1));
    tsd2= std(data(i,batch_label==2));

    data(i,batch_label==1)=(data(i,batch_label==1)-tmp1)./tsd1;
    data(i,batch_label==2)=(data(i,batch_label==2)-tmp2)./tsd2;
end
clear i tmp1 tmp2 tsd1 tsd2
save('ChronicWoundDataBatchCorrected.mat')