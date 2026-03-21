function  plotMissing(data,c)

shiftedTD = shiftMissing(data{c}.indivArrayTD,1000);
quantsTD = quantile(shiftedTD,[0.025,0.5,0.975],1);
shiftedASD = shiftMissing(data{c}.indivArrayASD,1000);
quantsASD = quantile(shiftedASD,[0.025,0.5,0.975],1);
lw = 1.5;
green = [0.0977    0.6875    0.1055];
blue = [0.1094    0.1562    0.6758];

t = tiledlayout(3,1,'TileSpacing', 'compact','TileIndexing', 'columnmajor');
nexttile(t)
for p = 1:size(data{c}.indivArrayTD,1)
    hold on
    missInds = find(data{c}.indivArrayTD(p,:)==0);
    
    scatter(missInds, ...
        repelem(p,length(missInds)),2,blue,".")
end
ylabel("Participant")
ax = gca;
ax.XAxis.Visible = 'off';

nexttile(t)
for p = 1:size(data{c}.indivArrayASD,1)
    hold on
    missInds = find(data{c}.indivArrayASD(p,:)==0);
    scatter(missInds, ...
        repelem(p,length(missInds)),2,green,".")
end
ylabel("Participant")
ax = gca;
ax.XAxis.Visible = 'off';
nexttile(t)

hold on
time = 1:size(data{c}.indivArrayASD,2);
missTD = mean(data{c}.indivArrayTD==0,1)
missASD = mean(data{c}.indivArrayASD==0,1)


hold on
plot(missTD,"Color",blue,"LineWidth",lw)
plot(missASD,"Color",green,"LineWidth",lw)
plot(median(shiftedTD,1),"Color",blue,"LineStyle","-.","LineWidth",lw)

fill([time fliplr(time)], ...
    [quantsTD(1,:) fliplr(quantsTD(3,:))], ...
    blue, ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none')

plot(median(shiftedASD,1),"Color",green,"LineStyle","-.","LineWidth",lw)
fill([time fliplr(time)], ...
    [quantsASD(1,:) fliplr(quantsASD(3,:))], ...
    green, ...
    'FaceAlpha', 0.4, 'EdgeColor', 'none')

% Show correlation
[corr,~,low,high] = corrcoef(missTD, ...
   missASD)
text(100,0.6,strcat("\rho = ",sprintf("%.2f", corr(1,2)), ...
    " [",sprintf("%.2f",low(1,2)),", ", ...
    sprintf("%.2f",high(1,2)),"]"))

ylim([0 0.7])
xticks(time(1:100:end))
xticklabels(round(time(1:100:end)/30,1))
ylabel({"Proportion of participants", "with missing data"})
xlabel("Time (s)")

end