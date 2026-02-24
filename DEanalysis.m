clear
load ChronicWoundData.mat

%select samples
 k=find(batch_label==1);  filename = 'DE_results_batch1.xlsx'; sheet = 'All_samples';
% k=find(batch_label==2);  filename = 'DE_results_batch2.xlsx'; sheet = 'All_samples';
% k=find(batch_label==1 & treatment_label>1);filename = 'DE_results_batch1.xlsx'; sheet = 'Non_teated';
% k=find(batch_label==2 & treatment_label>1);filename = 'DE_results_batch2.xlsx'; sheet = 'Non_teated';
% k=find(batch_label==1 & treatment_label<2); filename = 'DE_results_batch1.xlsx'; sheet = 'Teated';
% k=find(batch_label==2 & treatment_label<2); filename = 'DE_results_batch2.xlsx'; sheet = 'Teated';

D=data(:,k);
cl=healed_label(k);

%number of zeros ineach group for each gene
for i=1:numel(Gene)
    nz0(i,1)=nnz(D(i,cl==0));
    nz1(i,1)=nnz(D(i,cl==1));
end

% create table for DE analysis
countTable = array2table(D, ...
    'VariableNames', strcat("S", string(1:numel(cl))), ...
    'RowNames', Gene);

cond1 = find(cl == 1);
cond2 = find(cl == 0);

diffTable = rnaseqde(countTable, cond1, cond2);

diffTable = addvars(diffTable, Gene, 'Before', 1, 'NewVariableNames', 'Gene');
diffTable = addvars(diffTable, nz0, 'After', 6, 'NewVariableNames', append('N non-zero non-healers (out of ',string(numel(cl)-sum(cl)),')'));
diffTable = addvars(diffTable, nz1, 'After', 7, 'NewVariableNames', append('N non-zero healers (out of ',string(sum(cl)),')'));

diffTable = renamevars(diffTable, ...
    {'Mean1','Mean2'}, ...
    {'Mean Healers','Mean Non-healers'});
DEupTable   = diffTable(diffTable.Log2FoldChange >= 0, :);
DEdownTable = diffTable(diffTable.Log2FoldChange <  0, :);

DEupTable   = sortrows(DEupTable,   'AdjustedPValue', 'ascend');
DEdownTable = sortrows(DEdownTable, 'AdjustedPValue', 'ascend');

DEupTable   = DEupTable(DEupTable.AdjustedPValue <= 0.01, :);
DEdownTable = DEdownTable(DEdownTable.AdjustedPValue <= 0.01, :);

DEupTable   = DEupTable(DEupTable.Log2FoldChange >= 1, :);
DEdownTable = DEdownTable(DEdownTable.Log2FoldChange <= -1, :);

DEupTable.Log2FoldChange(isinf(DEupTable.Log2FoldChange)) = NaN;
DEdownTable.Log2FoldChange(isinf(DEdownTable.Log2FoldChange)) = NaN;

writetable(DEupTable,   filename, 'WriteRowNames', true, 'Sheet', sheet, 'Range', 'A1');
writetable(DEdownTable, filename, 'WriteRowNames', true, 'Sheet', sheet, 'Range', 'K1');



% %% Some plots
% DEup  =find(diffTable.Log2FoldChange>1 & -log10(diffTable.AdjustedPValue)>2);
% DEdown=find(diffTable.Log2FoldChange<-1 & -log10(diffTable.AdjustedPValue)>2);
% [~,DEany]=mink(diffTable.AdjustedPValue,30);
% 
% figure
% plot(diffTable.Log2FoldChange,-log10(diffTable.AdjustedPValue),'k.')
% hold on
% plot(diffTable.Log2FoldChange(DEup),-log10(diffTable.AdjustedPValue(DEup)),'r.')
% plot(diffTable.Log2FoldChange(DEdown),-log10(diffTable.AdjustedPValue(DEdown)),'r.')
% plot(diffTable.Log2FoldChange(DEany),-log10(diffTable.AdjustedPValue(DEany)),'bo')
% xlabel('Log_2 ( FoldChange )')
% ylabel('-Log_1_0 ( AdjustedPValue )')
% grid on
% 
% 
% n1=sum(cl);
% n2=numel(cl)-n1;
% Gup=Gene(DEup);
% Gdown=Gene(DEdown);
% figure
% for i=1:min(24,numel(Gup))
%     subplot(4,6,i)
%     X=log2(1+D(DEup(i),:));
%     boxplot(X,cl)
%     hold on
%     plot(2*ones(1,n1),X(cl==1),'b.')
%     plot(ones(1,n2),X(cl==0),'b.')
%     title(Gup(i))
%     xticklabels({'Nh','H'})
% end
% 
% figure
% for i=1:min(24,numel(Gup))
%     subplot(4,6,i)
%     X=log2(1+D(DEdown(i),:));
%     boxplot(X,cl)
%     hold on
%     plot(2*ones(n1,1),X(cl==1),'b.')
%     plot(ones(n2,1),X(cl==0),'b.')
%     title(Gup(i))
%     xticklabels({'Nh','H'})
% end
