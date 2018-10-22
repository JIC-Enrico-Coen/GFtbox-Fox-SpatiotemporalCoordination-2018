function processModelGraphs()
%This function plots the graphs used in the paper for the model data
%using output from GFtbox.
close all
%% Plot graphs for the early stage spch Epidermis model
plotModelGraphs();
% this function is used to generate the following figures:
% Fig 4E, 4G
% Fig 5D

%% Plot graphs for the early stage spch SubEpidermis model
plotModelGraphs2();
% this function is used to generate the following figures:
% Fig 4F, 4H
end

function plotModelGraphs()
clear all
graphTitles = 'EPIMODEL';
meshDir = 'Z:\work2\Sam\GPT_arLeafCellDiv_V21_20170620\runs\epimodel\meshes\';
fnamelist{7} = 'GPT_arLeafCellDiv_V21_20170620_s000178.mat';
fnamelist{6} = 'GPT_arLeafCellDiv_V21_20170620_s000164.mat';
fnamelist{5} = 'GPT_arLeafCellDiv_V21_20170620_s000156.mat';
fnamelist{4} = 'GPT_arLeafCellDiv_V21_20170620_s000147.mat';
fnamelist{3} = 'GPT_arLeafCellDiv_V21_20170620_s000140.mat';
fnamelist{2} = 'GPT_arLeafCellDiv_V21_20170620_s000132.mat';
fnamelist{1} = 'GPT_arLeafCellDiv_V21_20170620_s000115.mat';

leafWidths = [0.15169 0.22171 0.26626 0.31326 0.38803 0.47258 0.67381];
leafWidthMin = 0.095;
leafWidthMax = 0.68;
colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax);
PSTarray = [123 143 152 159 172 184 203 208];

%% Generate Fig 4E and 4G
scatter_cellArea_lpboundry(meshDir,fnamelist,graphTitles,colourArray)
%% Generate Fig 5D
scatter_growthRate_cellArea(meshDir,fnamelist,graphTitles,colourArray,PSTarray)
end

function plotModelGraphs2()
clear all
graphTitles = 'SUBEPIMODEL';
meshDir = 'Z:\work2\Sam\GPT_arLeafCellDiv_V21_20170620\runs\subepimodel\meshes\';
fnamelist{7} = 'GPT_arLeafCellDiv_V21_20170620_s000178.mat';
fnamelist{6} = 'GPT_arLeafCellDiv_V21_20170620_s000164.mat';
fnamelist{5} = 'GPT_arLeafCellDiv_V21_20170620_s000156.mat';
fnamelist{4} = 'GPT_arLeafCellDiv_V21_20170620_s000147.mat';
fnamelist{3} = 'GPT_arLeafCellDiv_V21_20170620_s000140.mat';
fnamelist{2} = 'GPT_arLeafCellDiv_V21_20170620_s000132.mat';
fnamelist{1} = 'GPT_arLeafCellDiv_V21_20170620_s000115.mat';

leafWidths = [0.15169 0.22171 0.26626 0.31326 0.38803 0.47258 0.67381];
leafWidthMin = 0.095;
leafWidthMax = 0.68;
colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax);
PSTarray = [123 143 152 159 172 184 203 208];

%% Generate Fig 4F and 4H
scatter_cellArea_lpboundry(meshDir,fnamelist,graphTitles,colourArray)
end



function scatter_cellArea_lpboundry(meshDir,fnamelist,graphTitles,colourArray)
%cmap = colormap(jet(numel(fnamelist)));
cmap = colourArray;
ylimMax = 2700; %epi 2700 subepi 900
ylimMin = 0;
xlimMax = 800;
xlimMin = -100;

%plot normal scatter plot Lamina Cells against distance lp-boundry
figure; hold on
for i = numel(fnamelist):-1:1
    load([meshDir,fnamelist{i}]);
    petIdx = find(not(cellfun('isempty', strfind(m.mgenIndexToName,'V_JUSTPET'))));
    somenodes = m.nodes(m.morphogens(:,petIdx)==1,:,:);
    displace = max(somenodes(:,2));    
    [centroids, areas, projections] = getCellAreasAndPositions( m,'morphogen', 'v_justlam', 'threshold', 0.5 , 'mode', 'ave','axis',[0 1 0]);
    areas = areas .* 1e6;
    projections = projections .* 1000;
    displace = displace .* 1000;
    projections = projections - displace;
    plot( projections, areas, 'Marker','o','MarkerFaceColor',cmap(i,:), ...
        'LineStyle','none','MarkerEdgeColor',[0 0 0] )
end
ylabel('Cell Area um2')
xlabel('Distance from lamina-petiole boundry in um')
legend('115h','132h','140h','147h','156h','164h','178h', ...
    'location','northwest');
ylim([ylimMin ylimMax]);
xlim([xlimMin xlimMax]);
title(['Model Data - Lamina Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

%plot normal scatter plot Midvein cells against distance lp-boundry
figure; hold on
for i = numel(fnamelist):-1:1
    load([meshDir,fnamelist{i}]);
    petIdx = find(not(cellfun('isempty', strfind(m.mgenIndexToName,'V_JUSTPET'))));
    somenodes = m.nodes(m.morphogens(:,petIdx)==1,:,:);
    displace = max(somenodes(:,2));    
    [centroids, areas, projections] = getCellAreasAndPositions( m,'morphogen', 'v_justmid', 'threshold', 0.5 , 'mode', 'ave','axis',[0 1 0]);
    areas = areas .* 1e6;
    projections = projections .* 1000;
    displace = displace .* 1000;
    projections = projections - displace;
    plot( projections, areas, 'Marker','o','MarkerFaceColor',cmap(i,:), ...
        'LineStyle','none','MarkerEdgeColor',[0 0 0] )
end
ylabel('Cell Area um2')
xlabel('Distance from lamina-petiole boundry in um')
legend('115h','132h','140h','147h','156h','164h','178h', ...
    'location','northeast');
ylim([ylimMin ylimMax]);
xlim([xlimMin xlimMax]);
title(['Model Data - Midvein Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end

function scatter_growthRate_cellArea(meshDir,fnamelist,graphTitles,colourArray,PSTarray)
cmap = colourArray;
printTxt = 0;
load([meshDir,fnamelist{end}])
m2 = m;

%plot growth rate against cell area for non-distal cells in color and distal cells in gray (PST time)
for i = 1:numel(fnamelist)-1
    load([meshDir,fnamelist{i}]);
    figure; hold on
    distIdx = find(not(cellfun('isempty', strfind( m.secondlayer.valuedict.index2NameMap,'c_justverydist'))));   
    petIdx = find(not(cellfun('isempty', strfind( m.secondlayer.valuedict.index2NameMap,'c_justnotpet'))));

%     cellsOfTintrestPtr = find(m.secondlayer.cellvalues(:,distIdx) > 0);
%     allCellAreasAtStart = m2.userdatastatic.cellareastart{i} .* 1e6;
%     cellAreasAtStart = allCellAreasAtStart(cellsOfTintrestPtr);
%     allCellAreasAtEnd = m2.userdatastatic.cloneareas{i};
%     cellAreasAtEnd = allCellAreasAtEnd(cellsOfTintrestPtr) .* 1e6;    
%     timeInterval = PSTarray(i+1) - PSTarray(i);
%     relativeGrowthRate = (log(cellAreasAtEnd) - log(cellAreasAtStart)) / timeInterval;    
%     plot(log(cellAreasAtStart),relativeGrowthRate, 'Marker','o','MarkerFaceColor',[0.7 0.7 0.7], ...
%          'LineStyle','none','MarkerEdgeColor',[0 0 0]);
     
    %cellsOfTintrestPtr = find((m.secondlayer.cellvalues(:,distIdx) < 1) & (m.secondlayer.cellvalues(:,petIdx) > 0) > 0);
    cellsOfTintrestPtr = find((m.secondlayer.cellvalues(:,petIdx) > 0) > 0);
    
    allCellAreasAtStart = m2.userdatastatic.cellareastart{i} .* 1e6;
    cellAreasAtStart = allCellAreasAtStart(cellsOfTintrestPtr);
    allCellAreasAtEnd = m2.userdatastatic.cloneareas{i};
    cellAreasAtEnd = allCellAreasAtEnd(cellsOfTintrestPtr) .* 1e6;    
    timeInterval = PSTarray(i+1) - PSTarray(i);
    relativeGrowthRate = (log(cellAreasAtEnd) - log(cellAreasAtStart)) / timeInterval;    
    plot(log(cellAreasAtStart),relativeGrowthRate, 'Marker','o','MarkerFaceColor',cmap(i,:), ...
         'LineStyle','none','MarkerEdgeColor',[0 0 0]);     
    ylim([0 0.1]);
    xlim([0 8]);
    
    [fitobject,gof,output] = fit(log(cellAreasAtStart), relativeGrowthRate,'poly1');
    lm = fitlm(log(cellAreasAtStart), relativeGrowthRate,'linear');
    t = lm.Coefficients;
    pval = t{2,4};
    gradient = t{2,1};
    meanCA = mean(log(cellAreasAtStart));
    SE = std(log(cellAreasAtStart)) / (sqrt(numel(cellAreasAtStart)));    
    meankA = mean(relativeGrowthRate);
    SEkA = std(relativeGrowthRate) / (sqrt(numel(relativeGrowthRate)));    
    plot([meanCA,meanCA],[0,0.1],'r','linewidth',3);
    plot([0,8],[meankA,meankA],'r--','linewidth',3);

    p = polyfit(log(cellAreasAtStart), relativeGrowthRate, 1);
    yfit = polyval(p,log(cellAreasAtStart));
    if(i~=1)
        plot(log(cellAreasAtStart),yfit,'Linewidth',3,'Color',[1 0 1])
    end
    if(printTxt)
        text(0.25,0.1,['rsquare=',num2str(gof.rsquare), ' p-val=',num2str(pval), ' m=',num2str(gradient)]);
        text(0.25,0.085,['\mu x=',num2str(meanCA), ' SE = ', num2str(SE), ' SEx1.96 = ', num2str(SE*1.96)]);
        text(0.25,0.09,[' \mu y=', num2str(meankA), ' SE = ', num2str(SEkA), ' SEx1.96 = ', num2str(SEkA*1.96)]);
        ylabel(['Karea h^1 (T interval ',num2str(timeInterval),')'])
        xlabel('Ln Cell area um^2')
    else
        ylim([0 0.1]);
        xlim([0 8]);
    end    
    title(['Model Data - All Cells','- ',graphTitles,' PST']);
    g = gca;
    printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'-',num2str(i),'.png'];
    print(printname,'-dpng','-r300');
        
end
end

function colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax)
numOfIntervals = 100;
cmap = colormap(jet(numOfIntervals));
scaleMap = linspace(leafWidthMin,leafWidthMax,numOfIntervals);

colourArray = zeros(numel(leafWidths),3);
for i = 1:numel(leafWidths)
    [~,idx] = min(abs(scaleMap-leafWidths(i)));
    closest(i) = scaleMap(idx);
    colourArray(i,:) = cmap(idx,:);
end

figure;
colormap(cmap)
tlabels = round(linspace(leafWidthMin,leafWidthMax,11),2);
tlabels = strread(num2str(tlabels),'%s');
colorbar('Location','northoutside','TickLabels',tlabels)
title('leaf width colour bar');
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end