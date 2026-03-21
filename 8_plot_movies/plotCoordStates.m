function plotCoordStates(flattened,clusterArray)

t = 1:size(flattened,2);

flattened_interp = flattened;  % copy

for i = 1:size(flattened,1)
    y = flattened(i,:);
    valid = y > 0 & ~isnan(y);

    if sum(valid) >= 2
        % interpolate only where y <= 0
        y_interp = interp1(t(valid), y(valid), t, 'linear', nan);
        y(~valid) = y_interp(~valid);
    end

    flattened_interp(i,:) = y;
end

[~,~,clusters] = unique(clusterArray);
uClust = unique(clusters);

cmap = [
    1 0 0;
    0 0 1;
    0 0.5 0;
    1 0.647 0];

hold on

for i = 1:size(flattened_interp,1)
    col = cmap(uClust == clusters(i), :);
    plot(t, flattened_interp(i,:), ...
         'Color', [col 0.2], 'LineWidth', 1)
end

for k = 1:numel(uClust)
    idx = clusters == uClust(k);

    meanTrace = mean(flattened_interp(idx,:), 1, 'omitnan');
    if sum(mode(flattened_interp(idx,:),1)==0) > 0.8*length(t)
        meanTrace = repelem(0,length(t));
    end
    plot(t, meanTrace, ...
         'Color', cmap(k,:), ...
         'LineWidth', 3);
end

for k = 1:numel(uClust)
    idx = clusters == uClust(k);
    mu  = median(flattened_interp(idx,:),1,'omitnan');
    sem = std(flattened_interp(idx,:),0,1,'omitnan');
    
    if sum(mode(flattened_interp(idx,:),1)==0) <= 0.8*length(t)
        fill([t fliplr(t)], ...
             [mu-sem fliplr(mu+sem)], ...
             cmap(k,:), ...
             'FaceAlpha', 0.4, 'EdgeColor', 'none');
    end
end
ylim([0 2000])
axis off
hold off

end