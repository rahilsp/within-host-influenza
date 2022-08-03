%% Parameter correlation heatmap
% z = 1: Healthy / Pooled
% z = 2: Suppressed
z = 1;

if z == 1
    labels = {'$\beta$','$p$','$\tau$','$s$','$r$','$f$','$k$','$V_0$'};
elseif z == 2
    labels = {'$\beta$','$p$','$\tau$','$s$','$V_0$'};
end

l = size(theta_vec,2);
corr_mat = NaN*ones(l,l);

for i = 1:l
    for j = 1:l
        if i ~= j
         corr_mat(i,j) = corr(theta_vec(:,i),theta_vec(:,j),'Type','Spearman');
%          corr_mat(i,j) = corr(theta_vec(:,i),theta_vec(:,j),'Type','Pearson');
        end
    end
    
    corr_mat = round(corr_mat,3,'decimals');
end


figure
xvalues = 1:l;
h = heatmap(xvalues,xvalues,corr_mat,'Fontsize',14);
h.Colormap = jet;
h.YDisplayLabels = repmat({''}, size(h.YData));  %remove row labels
h.XDisplayLabels = repmat({''}, size(h.XData));  %remove row labels
a2 = axes('Position', h.Position);               %new axis on top
a2.Color = 'none';                               %new axis transparent
a2.YTick = 1:size(h.ColorData,1);                %set y ticks to number of rows
a2.XTick = 1:size(h.ColorData,1);                %set x ticks to number of rows                   
ylim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks
xlim(a2, [0.5, size(h.ColorData,1)+.5])          %center the tick marks
a2.YDir = 'Reverse';                             %flip your y axis to correspond with heatmap's
a2.YTickLabel = labels;     %set you ytick labels here, formatted.
a2.XTickLabel = labels;     %set you xtick labels here, formatted.
a2.FontSize = 22;
h.ColorLimits = [-1,1];
h.MissingDataLabel = {'N/A'};
if z == 1
set(gcf,'Position',[100 100 650 500])
elseif z == 2
set(gcf,'Position',[100 100 550 400])
end
shg

%% 
sav = 0;
i = 5;
j = 6;

figure
scatter(theta_vec(:,i),theta_vec(:,j),6,'filled')
% scatplot(theta_vec(:,i),theta_vec(:,j),[],[],[],[],1,3);
colorbar off
xlabel(labels(i),'Fontsize',30)
ylabel(labels(j),'Fontsize',30)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
box on

% if sav == 1
% filename = [char(label1(i)),'-',char(label1(j))]
% saveas(gcf,filename,'epsc')
% end

