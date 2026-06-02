load ChronicWoundData.mat
T=readtable('S11_Table_qPCR_results.xlsx');
SampleName=T{:,1};
for i=1:numel(SampleName)
    n(i)=find(ismember(Metadata_filename,SampleName(i)));
end

genes=T.Properties.VariableNames;
genes(1)=[];
%remove discrepancies in gene names between RNAseq and qPCR datsets
genes=upper(genes);
tmp=find(ismember(Gene,'LORICRIN'));   Gene(tmp)={'LOR'};
tmp=find(ismember(genes,'SPRRIB'));     genes(tmp)={'SPRR1B'};

for i=1:numel(genes)
    k(i)=find(ismember(Gene,genes(i)));
end

genexp=T{:,2:10};
DTP=data(k,n)';

%re-arrange samples
tmp=[5 6 10 2 3 4 1 8 9 7];
DTP=DTP(tmp,:);
genexp=genexp(tmp,:);
SampleName=SampleName(tmp);


geneorder={'IL1B' 'KRT6B'  'KRT10'...
           'CCL3' 'KRT6A'  'KRT1' ...
           'OSM'  'SPRR1B' 'LOR'};
fig_qpcr=figure;
for ii=1:numel(geneorder)
    colororder({'k','k'})
    i=find(ismember(genes,geneorder(ii)));
    D2=[(DTP(:,i)) [genexp(:,i).*1]];
    nil = zeros(10,1);
    subplot(3,3,ii)
    p1=bar(SampleName,[D2(:,1),nil]);
    ylabel('RNAseq');
    yyaxis right
    p2=bar(SampleName,[nil,D2(:,2)], 'grouped','r');
    ylabel('qPCR');
    grid on
    if ii==3
        lgd=legend([p1(2),p2(1)],'RNAseq','qPCR');
        lgd.FontSize = 8;
    end
    xtickangle(90)
    title(genes(i))
    set(gca,'fontsize', 8) 
end
exportgraphics(fig_qpcr,'qpcr_val_barplots.pdf','ContentType','vector')
