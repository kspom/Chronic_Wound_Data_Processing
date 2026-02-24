load ChronicWoundData.mat

cl=healed_label;
data=log2(1+data);

% PCA for one batch only:
% cl=healed_label(batch_label==1);
% data=log2(1+data(:,batch_label==1));

[coeff,score,latent,tsquared,explained] = pca(data');
figure
scatter3(score(cl==1,1),score(cl==1,2),score(cl==1,3),'d',...
    'MarkerEdgeColor','#3070B7','MarkerFaceColor','#3070B7')
hold on
scatter3(score(cl==0,1),score(cl==0,2),score(cl==0,3),...
    'MarkerEdgeColor','k','MarkerFaceColor','#DD3B22')
% scatter3(score(batch_label==1,1),score(batch_label==1,2),score(batch_label==1,3),...
%     'MarkerEdgeColor','c','MarkerFaceColor','none')
axis equal
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
legend('Healer','Non-healer')


view([15 30])
exportgraphics(gcf,'scatter.pdf','ContentType','vector')