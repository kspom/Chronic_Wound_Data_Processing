% plotting boxplots for selected genes
clear
load ChronicWoundData.mat
% gene list can be assigned directly:
% GeneList={'KRT16','IL1B','IL6'}; filename='BoxplotsSomeGenes.pdf'
% GeneList={'MMP12','MMP3','IL24'}; 
% filename='downregulated_in_nonhealers4.pdf';

% reading gene list from file
T=readtable('DE_batches_comparison.xlsx','Sheet','Lists_for_reading');
GeneList1=T{:,2};
GeneList=GeneList1(10:18);
filename='down_in_nonhealers2.pdf';

treat2 = treatment_label;
treat2(treat2 == 3) = 2;

group = ...
    (batch_label-1)*4 + ...
    healed_label*2 + ...
    treat2;

labels = { ...
    'NT', 'NC', ... % non-healer treated and control
    'HT', 'HC', ... % healer treated and control
    'NT', 'NC', ... % non-healer treated and control
    'HT', 'HC'};    % healer treated and control

fobj=figure;
for i=1:numel(GeneList)
    subplot(3,3,i)

    k=find(ismember(Gene,GeneList(i)));
    X=log2(1+data(k,:));

    h=boxplot(X, group, 'Labels', labels);
    ylabel('log_2(1+expr)')
    title(GeneList(i))
    xtickangle(45)
    colorboxplot(h)
    xlabel({'Batch 1          Batch 2'})
end

exportgraphics(fobj,filename,'ContentType','vector');

function f=colorboxplot(h)
% Define colors you want for each box
% Each row is RGB [R G B] between 0 and 1
colors = [ 0.8980    0.6431    0.5922;   
          0.7961    0.2863    0.1882;   
          0.6275    0.7137    0.8471;   
          0.2549    0.4314    0.6980;
          0.8980    0.6431    0.5922;   
          0.7961    0.2863    0.1882;   
          0.6275    0.7137    0.8471;   
          0.2549    0.4314    0.6980];

% Get the patch objects
boxes = findobj(h, 'Tag', 'Box');

% Apply colors
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors(j,:), 'FaceAlpha',0.8);
end
end