function plotSilhouettes(silsAll,clips)

dat   = [];
groups = [];
colors = [];

for i = 1:5
    for j = 1:2
        v = silsAll{i,j}(:);
        g = (i-1)*2 + j;

        dat   = [dat; v];
        groups = [groups; repmat(g, numel(v), 1)];

        if j == 1
            colors = [colors; repmat([0 0.4470 0.7410], numel(v), 1)]; % blue
        else
            colors = [colors; repmat([0.4660 0.6740 0.1880], numel(v), 1)]; % green
        end
    end
end

figure; hold on;

colors = {
    [0.1094    0.1562    0.6758],
    [0.0977    0.6875    0.1055]
};


pos   = 1;
width = 0.35;

for i = 1:length(clips)
    for j = 1:2

        y = silsAll{i,j}(:);
        y = y(~isnan(y));

        ygrid = linspace(min(y), max(y), 200);
        f = ksdensity(y, ygrid, 'Function', 'pdf');


        f = f / max(f) * width;


        x = [pos - f, fliplr(pos + f)];
        yv = [ygrid, fliplr(ygrid)];

        patch(x, yv, colors{j}, ...
            'FaceAlpha', 0.4, ...
            'EdgeColor', 'none');


        jitter = (rand(size(y)) - 0.5) * width * 0.7;
        scatter(pos + jitter, y, 12, ...
            'MarkerFaceColor', colors{j}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6);


        q1  = prctile(y, 25);
        med = median(y);
        q3  = prctile(y, 75);


        rectangle('Position', ...
            [pos - width*0.15, q1, width*0.3, q3-q1], ...
            'EdgeColor', 'k', ...
            'LineWidth', 1.2);


        plot([pos - width*0.15, pos + width*0.15], ...
             [med med], 'k-', 'LineWidth', 2);

        pos = pos + 1;
    end
end

xlim([0.5 10.5])
xticks(1.5:2:10.5)
xticklabels(1:5)
xtickangle(45)

xlabel('Movie clip')
ylabel('Median silhouette score')
box on

end

