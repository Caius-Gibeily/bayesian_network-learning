function plotEntropy(H,Hshuf_stats)

lower = Hshuf_stats(:,1);
upper = Hshuf_stats(:,3);
y = (1:length(lower))';

diff_h = H(:,1)-H(:,2);
below = diff_h < lower;           
above = diff_h > upper;          
inside = ~(below | above);


hold on

area(1:size(H,1), upper, 'FaceColor',[0 0 1], 'FaceAlpha',0.3, 'EdgeColor','none')   
plot(Hshuf_stats(:,2),"black", 'LineWidth',2)
area(1:size(H,1), lower, 'FaceColor',[0 1 0], 'FaceAlpha',0.3, 'EdgeColor','none')   

stem(y(below),diff_h(below),   ...
    'filled', 'Color',[0 1 0])

stem( y(above),diff_h(above),   ...
    'filled', 'Color',[0 0 1])

stem(y(inside), diff_h(inside),  ...
    'filled', 'Color',[0 0 0])

xlim([1 length(y)])


end