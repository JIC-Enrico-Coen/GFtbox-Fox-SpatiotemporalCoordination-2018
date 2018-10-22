function processGraphs()
%This function plots the graphs used in the paper for the expermental data
%using output from Track n R.
close all

%% Plot graphs for epidermis using data from expID 3002
plotExperimentalData();
%this function will generate the following figures:
% Fig 4A , 4C
% Fig 5B
% Fig Fig 2sup1
% Fig 4I

%% Plot graphs for the subepidermis using data from expID 3002
plotExperimentalData2();
%this function will generate the following figures:
% Fig 4B, 4D
% Fig 4J
% Fig 5C

%% Plot graphs for epidermis using data from exp ID 3289 (finer time resolution)
plotExperimentalData3();
%this function will generate the following figures:
% Fig 3B, 3D
% Fig 3F
% Fig 3C, 3E
% Fig 3G

%% Plot graphs for Wild-Type epidermis using data from exp ID 3278 
plotExperimentalData5();
%this function will generate the following figures:
% Fig 5A
% Fig 8A sup2 and Fig 8B sup2
end

function plotExperimentalData()
%exp 3002 cell areas epidermis
clear all
graphTitles = 'EXP 3002 - Epidermis';
fileAndLocation = '.\FP20170419\ExpID3002\ExpID3002_epidermis\data_ExpID3002_epid_FP20170418.csv';

fnamelist{1} = ['tracking_TL02-TL09_all_divide_new.csv',',','ExpID3002_spch4_TL002_plantA_scalebar.png',',','ExpID3002_spch4_TL004_plantA_scalebar.png'];
fnamelist{2} = ['tracking_TL04-TL09_all_divide_new.csv',',','ExpID3002_spch4_TL004_plantA_scalebar.png',',','ExpID3002_spch4_TL005_plantA_scalebar.png'];
fnamelist{3} = ['tracking_TL05-TL09_all_divide.csv',',','ExpID3002_spch4_TL005_plantA_scalebar.png',',','ExpID3002_spch4_TL006_plantA_scalebar.png'];
fnamelist{4} = ['tracking_TL06-TL09_all_divide.csv',',','ExpID3002_spch4_TL006_plantA_scalebar.png',',','ExpID3002_spch4_TL007_plantA_scalebar.png'];
fnamelist{5} = ['tracking_TL07-TL09_all_divide.csv',',','ExpID3002_spch4_TL007_plantA_scalebar.png',',','ExpID3002_spch4_TL008_plantA_scalebar.png'];
fnamelist{6} = ['tracking_TL08-TL09_all_divide.csv',',','ExpID3002_spch4_TL008_plantA_scalebar.png',',','ExpID3002_spch4_TL009_plantA_scalebar.png'];
fnamelist{7} = ['tracking_TL09-TL10_all_divide.csv',',','ExpID3002_spch4_TL009_plantA_scalebar.png',',','ExpID3002_spch4_TL010_plantA_scalebar.png'];

leafWidths = [0.15 0.22 0.265 0.3 0.385 0.48 0.68];
leafWidthMin = 0.095;
leafWidthMax = 0.68;
colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax);

%% Generate Fig 4A and 4C
scatter_cellArea_lpboundary_FlorentData(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 5B 
scatter_cellArea_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 2sup1
scatter_growthRate_distance(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 4I
%fancyBoxPlot_kAREA_midAndLam(fileAndLocation,fnamelist,graphTitles,colourArray);
normalBoxPlot(fileAndLocation,fnamelist,graphTitles,colourArray);
end

function plotExperimentalData2()
%exp 3002 cell areas sup-epi
clear all
graphTitles = 'EXP 3002 - Sub-epidermis';
fileAndLocation = '.\FP20170419\ExpID3002\ExpID3002_mesophyll\data_ExpID3002_meso_FP20170418.csv';

fnamelist{1} = ['TL02-TL09_divisions.csv', ',','ExpID3002_TL02_mesophyll.png',',','ExpID3002_TL04_mesophyll.png'];
fnamelist{2} = ['TL04-TL09_divisions.csv',',','ExpID3002_TL04_mesophyll.png',',','ExpID3002_TL05_mesophyll.png'];
fnamelist{3} = ['TL05-TL09_divisions.csv',',','ExpID3002_TL05_mesophyll.png',',','ExpID3002_TL06_mesophyll.png'];
fnamelist{4} = ['TL06-TL09_divisions.csv',',','ExpID3002_TL06_mesophyll.png',',','ExpID3002_TL07_mesophyll.png'];
fnamelist{5} = ['TL06-TL09_divisions.csv',',','ExpID3002_TL07_mesophyll.png',',','ExpID3002_TL08_mesophyll.png'];
fnamelist{6} = ['TL06-TL09_divisions.csv',',','ExpID3002_TL08_mesophyll.png',',','ExpID3002_TL09_mesophyll.png'];
fnamelist{7} = ['TL09-TL10_divisions.csv',',','ExpID3002_TL09_mesophyll.png',',','ExpID3002_TL10_mesophyll.png'];

leafWidths = [0.15 0.22 0.265 0.3 0.385 0.48 0.68];
leafWidthMin = 0.095;
leafWidthMax = 0.68;
colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax);
%% Generate Fig 4B and 4D
scatter_cellArea_lpboundary_FlorentData(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 4J
%fancyBoxPlot_kAREA_midAndLam(fileAndLocation,fnamelist,graphTitles,colourArray);
normalBoxPlot(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 5C
scatter_cellArea_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray);
end

function plotExperimentalData3()
%exp 3289 cell area at divison epidermis
clear all
graphTitles = 'EXP 3289 - Epidermis';
fileAndLocation = '.\FP20170419\ExpID3289\data_ExpID3289_FP20170309.csv';

fnamelist{1} = ['tracking_all_TL00-TL05+petiole_edited.csv,','TL00.png',',','TL02_t04_despeckle_median1.png'];
fnamelist{2} = ['tracking_all_TL00-TL05+petiole_edited.csv',',','TL02_t04_despeckle_median1.png',',','TL05_despeckle_median1.png'];
fnamelist{3} = ['tracking_All_TL05-TL12+petiole_edited.csv',',','TL05_despeckle_median1.png',',','TL07_t04.png'];
fnamelist{4} = ['tracking_All_TL05-TL12+petiole_edited.csv',',','TL07_t04.png',',','TL12_despeckle_median1.png'];
fnamelist{5} = ['tracking_All_TL12-TL20_edited.csv',',','TL12_despeckle_median1.png',',','TL16_t04.png'];
fnamelist{6} = ['tracking_All_TL12-TL20_edited.csv',',','TL16_t04.png',',','TL20.png'];

%leafWidths = [0.15 0.22 0.265 0.3 0.385 0.48 0.68]; 
%leafWidths = [0.095 0.122 0.180 0.260 0.356 0.470 0.552];
leafWidths = [0.095 0.095 0.180 0.180 0.360 0.360 0.552]; 
leafWidthMin = 0.095;
leafWidthMax = 0.68;
colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax);
%% Generate Fig 3B and 3D
scatter_cellAreaDiv_distance(fileAndLocation,fnamelist,graphTitles,colourArray); 
%% Generate Fig 3F
hist_cellAreaDiv_frequency(fileAndLocation,fnamelist,graphTitles,colourArray,1);
%% Generate Fig 3C and 3E
scatter_cellCycleDur_distance(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 3G
hist_cellCycleDur_frequency(fileAndLocation,fnamelist,graphTitles,colourArray,1);

end

function plotExperimentalData5()
%exp 3278 WT data
clear all
graphTitles = 'EXP 3278 - WT Epidermis';
fileAndLocation = '.\FP20170419\ExpID3278\data_ExpID3278_FP20170312_EDIT.csv';

fnamelist{1} = ['tracking_samJan17.csv',',','TL01_median.png',',','TL09_T03_median_flip.png'];
fnamelist{2} = ['tracking_samJan17.csv',',','TL09_T03_median_flip.png',',','TL10_median.png'];
fnamelist{3} = ['tracking_samJan17.csv',',','TL10_median.png',',','TL19_T03_median_flip.png'];
fnamelist{4} = ['tracking_samJan17.csv',',','TL19_T03_median_flip.png',',','TL19_T13_median_flip.png'];
fnamelist{5} = ['tracking_samJan17.csv',',','TL19_T13_median_flip.png',',','TL26_median.png'];

leafWidths = [0.170 0.230 0.282 0.392 0.42];
leafWidthMin = 0.095;
leafWidthMax = 0.68;
colourArray = convertColours(leafWidths, leafWidthMin,leafWidthMax);

%% Generate Fig 5A
scatter_WT_cellArea_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray);
%% Generate Fig 8A sup2 and Fig 8B sup2
scatter_WT_cellAreaDiv_distance(fileAndLocation,fnamelist,graphTitles,colourArray);

end

function [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(csvString,startString,endString,fileAndLocation)
t = readtable(fileAndLocation);
%find the idxs in ther table where csvString is present
csvStringIdxs =  strcmp(t.csvFile,csvString);
fistIdx = find(csvStringIdxs,1,'first');
lastIdx = find(csvStringIdxs,1,'last');

%truncate table to include
tt = t(fistIdx:lastIdx,:);

%find the startString Idxs
startStringIdxs =  strcmp(tt.Image1,startString);
startStringFistIdx = find(startStringIdxs,1,'first');
startStringLastIdx = find(startStringIdxs,1,'last');

%find the endString Idxs
endStringIdxs =  strcmp(tt.Image2,endString);
endStringFistIdx = find(endStringIdxs,1,'first');
endStringLastIdx = find(endStringIdxs,1,'last');

%find the union of startStringIdxs & endStringIdxs
unionIdxs = startStringIdxs .* endStringIdxs;
ptr = find(unionIdxs == 1);

%get x and y distances
yDists = tt.CentroY(ptr);
xDists = tt.CentroX(ptr);

%get cellAreas
cellAreas = tt.Area(ptr);

%is cell MID or not
ismid = zeros(size(ptr,1),1);
for i = 1:numel(ptr)
    if(strcmp(tt.IsMID(ptr(i)),'TRUE'))
        ismid(i) = 1;
    else
        ismid(i) = 0;
    end
end
%iscomp
iscomp = zeros(size(ptr,1),1);
for i = 1:numel(ptr)
    if(strcmp(tt.IsCompetent(ptr(i)),'TRUE'))
        iscomp(i) = 1;
    else
        iscomp(i) = 0;
    end
end
%hasDiv
hasDiv = zeros(size(ptr,1),1);
for i = 1:numel(ptr)
    if(strcmp(tt.HasDivided(ptr(i)),'TRUE'))
        hasDiv(i) = 1;
    else
        hasDiv(i) = 0;
    end
end
%cellAreaAtDiv
%cellAreaAtDiv = tt.AreaAtDiv(ptr);
cellAreaAtDiv = zeros(size(ptr,1),1);
for i = 1:numel(ptr)
    if(strcmp(tt.AreaAtDiv(ptr(i)),'NA'))
        cellAreaAtDiv(i) = NaN;
    else
        cellAreaAtDiv(i) = str2double(tt.AreaAtDiv(ptr(i)));
    end
end

cellCycleDuration = zeros(size(ptr,1),1);
for i = 1:numel(ptr)
    if(strcmp(tt.Cell_Cycle_Duration(ptr(i)),'NA'))
        cellCycleDuration(i) = NaN;
    else
        cellCycleDuration(i) = str2double(tt.Cell_Cycle_Duration(ptr(i)));
    end
end

%kArea
%kArea = tt.karea(ptr);
kArea = zeros(size(ptr,1),1);
for i = 1:numel(ptr)
    if(strcmp(tt.karea(ptr(i)),'NA'))
        kArea(i) = NaN;
    else
        if(isType(tt.karea(ptr(i)),'double'))
            kArea(i) = tt.karea(ptr(i));
        else
            kArea(i) = str2double(tt.karea(ptr(i)));
        end
    end
end


%kPAR
kPAR = tt.ky(ptr);

%kPER
kPER = tt.kx(ptr);

%Identity
if isfield(table2struct(tt),'Identity')
    identity = tt.Identity(ptr);
else
    identity = ones(size(ptr,1),1) .* -99;
end

% xDists = xDists ./ 1000;
% yDists = yDists ./ 1000;
% cellAreas = cellAreas ./ 1000;
% cellAreaAtDiv = cellAreaAtDiv ./ 1000;

end

function hist_WT_cellArea_Distance_FlorentData(fileAndLocation,fnamelist,graphTitles,colourArray)
%plot hist of All Cells vs lp boundary 
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv, ...
        cellAreaAtDiv,kArea,kPAR,kPER, identity] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    projections = yDists;
    cellAreasThatMeetSpec = cellAreas;
    bins = -400:50:800;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(cellAreasThatMeetSpec(idx),'omitnan');        
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area um2')
xlabel('Distance from lamina-petiole boundary in um')
legend('TL01-09T03','TL09T03-10','TL10-19T03','TL19T03-19T15',...
    'TL19T15-26', ...
    'location','northwest');
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;

title(['Realworld WT Data - All Cells','- ',graphTitles]);

ylim([0 1400]);
xlim([0 800]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
%print(printname,'-dpng','-r300');

%plot hist of cell area at division vs distance midvein SORT THIS OUT
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find(hasDiv == 1);
    projections = xDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    bins = -150:25:100;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area at Division um2')
xlabel('Distance from midvein in um')
legend(legendTxt,'location','northwest','Interpreter', 'none');

title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 1400]);
xlim([-200 100]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];
end





function scatter_cellArea_lpboundary_FlorentData(fileAndLocation,fnamelist,graphTitles,colourArray)
ylimMax = 2700; %epi 2700 subepi 900
ylimMin = 0;
xlimMax = 800;
xlimMin = -100;
%cmap = colormap(parula(numel(fnamelist)));
%cmap = colormap(jet(numel(fnamelist)));
cmap = (colourArray);

%plot normal scatter plot Lamina Cells
figure; hold on
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    lamPtr = find(ismid == 0);
    plot(yDists(lamPtr),cellAreas(lamPtr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
        'LineStyle','none','MarkerEdgeColor',[0 0 0]);    
end
ylabel('Cell Area um2')
xlabel('Distance from lamina-petiole boundary in um')
% legend('TL02-04 (124h)','TL04-05 (140h)','TL05-06 (147h)','TL06-07(153h)',...
%     'TL07-08(162h)','TL08-09(175h)','TL09-10(182h)', ...
%     'location','northwest');
legend('TL09-10','TL08-09','TL07-08','TL06-07',...
    'TL05-06','TL04-05','TL02-04','location','northwest');
ylim([ylimMin ylimMax]);
xlim([xlimMin xlimMax]);
title(['Realworld Data - Lamina Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');


%plot normal scatter plot Midvein Cells
figure; hold on
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    lamPtr = find(ismid == 1);
    plot(yDists(lamPtr),cellAreas(lamPtr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
        'LineStyle','none','MarkerEdgeColor',[0 0 0]);    
end
ylabel('Cell Area um2')
xlabel('Distance from lamina-petiole boundary in um')
legend('TL09-10','TL08-09','TL07-08','TL06-07',...
    'TL05-06','TL04-05','TL02-04','location','northeast');
ylim([ylimMin ylimMax]);
xlim([xlimMin xlimMax]);
title(['Realworld Data - Midvein Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

end



function scatter_cellCycleDur_distance(fileAndLocation,fnamelist,graphTitles,colourArray)
%cmap = colourArray;
cmap = colormap(lines(numel(fnamelist)));

%plot of duration of cell cycle vs distance from LP-boundary for Midvein
%cells
figure; hold on
allYs = [];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    %ptr = find(ismid==1);
    ptr = find((ismid==1) & (hasDiv==1) & ~isnan(cellCycleDuration));
    allYs = [allYs; cellCycleDuration(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),cellCycleDuration(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
meanAllY = mean(allYs);
SE = std(allYs) / (sqrt(numel(allYs)));
text(-75,5,['\mu=',num2str(meanAllY), ' SE = ', num2str(SE)]);
ylabel('Cell Cycle Duration (hours)')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northoutside');
ylim([0 50]);
xlim([-100 400]);
title(['Realworld Data - Midvein Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');



%plot of duration of cell cycle vs distance from LP-boundary for Lamina
%cells
figure; hold on
trendLineYs = [];
trendLineXs = [];
allYs = [];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((ismid==0) & (hasDiv==1) & ~isnan(cellCycleDuration));
    ptrGreater150 = find((ismid==0) & (hasDiv==1) & (yDists>=150) & ~isnan(cellCycleDuration));
    trendLineYs = [trendLineYs; yDists(ptrGreater150)];
    trendLineXs = [trendLineXs; cellCycleDuration(ptrGreater150)];
    allYs = [allYs; cellCycleDuration(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),cellCycleDuration(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(trendLineYs, trendLineXs, 1);
yfit = polyval(p,trendLineYs);
%plot(trendLineYs,yfit,'Linewidth',3,'Color',[1 0 0])
lm = fitlm(trendLineYs, trendLineXs,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
rsq = lm.Rsquared.Adjusted;
meanAllY = mean(allYs);
SE = std(allYs) / (sqrt(numel(allYs)));
% text(-75,5,['R^2=',num2str(rsq), ' p-val=',num2str(pval), ' m=',num2str(gradient), ' \mu=',num2str(meanAllY)]);
text(-75,5,['\mu=',num2str(meanAllY),' SE = ', num2str(SE)]);
ylabel('Cell Cycle Duration (hours)')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northoutside');
ylim([0 50]);
xlim([-100 400]);
title(['Realworld Data - Lamina Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

end

function scatter_WT_cellAreaDiv_distance(fileAndLocation,fnamelist,graphTitles,colourArray)
cmap = colormap(lines(numel(fnamelist)));

%plot scatter plot Lamina Pavement Cells cell area at divsion vs distance from
%lp-boundary
figure; hold on
trendLineYs = [];
trendLineXs = [];
allYs = [];
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & (ismid==0) & strcmp(identity, 'PC'));
    ptrGreater150 = find((hasDiv==1) & (ismid==0) & (yDists>=150) & strcmp(identity, 'PC'));
    trendLineYs = [trendLineYs; yDists(ptrGreater150)];
    trendLineXs = [trendLineXs; cellAreaAtDiv(ptrGreater150)];
    allYs = [allYs; cellAreaAtDiv(ptr)];
    if isempty(ptr)
        plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),cellAreaAtDiv(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(trendLineYs, trendLineXs, 1);
yfit = polyval(p,trendLineYs);
%plot(trendLineYs,yfit,'Linewidth',3,'Color',[1 0 0])
lm = fitlm(trendLineYs, trendLineXs,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
rsq = lm.Rsquared.Adjusted;
meanAllYs = mean(allYs);
SE = std(allYs) / (sqrt(numel(allYs)));
text(50,450,['mean=',num2str(meanAllYs),' SE=',num2str(SE)]);
ylabel('Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(fliplr(legendTxt),'location','northoutside');
ylim([0 500]);
xlim([0 400]);
title(['Realworld Data - Lamina Pavement Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

%plot scatter plot Lamina Non-Pavement Cells cell area at divsion vs distance from
%lp-boundary
figure; hold on
trendLineYs = [];
trendLineXs = [];
allYs = [];
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & (ismid==0) & ~strcmp(identity, 'PC'));
    ptrGreater150 = find((hasDiv==1) & (ismid==0) & (yDists>=150) & ~strcmp(identity, 'PC'));
    trendLineYs = [trendLineYs; yDists(ptrGreater150)];
    trendLineXs = [trendLineXs; cellAreaAtDiv(ptrGreater150)];
    allYs = [allYs; cellAreaAtDiv(ptr)];
    if isempty(ptr)
        plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),cellAreaAtDiv(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(trendLineYs, trendLineXs, 1);
yfit = polyval(p,trendLineYs);
%plot(trendLineYs,yfit,'Linewidth',3,'Color',[1 0 0])
lm = fitlm(trendLineYs, trendLineXs,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
rsq = lm.Rsquared.Adjusted;
meanAllYs = mean(allYs);
SE = std(allYs) / (sqrt(numel(allYs)));
text(50,450,['mean=',num2str(meanAllYs),' SE=',num2str(SE)]);
ylabel('Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(fliplr(legendTxt),'location','northoutside');
ylim([0 500]);
xlim([0 400]);
title(['Realworld Data - Lamina Non-Pavement Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

%%find mean for all Pavement and Non-Pavement Cells
allYs = [];
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & (ismid==0)); 
    allYs = [allYs; cellAreaAtDiv(ptr)];  
end
meanAllYs = mean(allYs)
SEallYs = std(allYs) / (sqrt(numel(allYs)))
end

function scatter_cellAreaDiv_distance(fileAndLocation,fnamelist,graphTitles,colourArray)

cmap = colormap(lines(numel(fnamelist)));
%cmap = colormap(parula(numel(fnamelist)));
%cmap = colourArray;

%plot scatter plot Lamina Cells cell area at divsion vs distance from
%lp-boundary
figure; hold on
trendLineYs = [];
trendLineXs = [];
lessThan150Ys =[];
allYs = [];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & (ismid==0));
    ptrGreater150 = find((hasDiv==1) & (ismid==0) & (yDists>=150));
    trendLineYs = [trendLineYs; yDists(ptrGreater150)];
    trendLineXs = [trendLineXs; cellAreaAtDiv(ptrGreater150)];
    ptrLessThan150 = find((hasDiv==1) & (ismid==0) & (yDists<150));
    lessThan150Ys = [lessThan150Ys; cellAreaAtDiv(ptrLessThan150)];
    allYs = [allYs; cellAreaAtDiv(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),cellAreaAtDiv(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(trendLineYs, trendLineXs, 1);
yfit = polyval(p,trendLineYs);
%plot(trendLineYs,yfit,'Linewidth',3,'Color',[1 0 0])
lm = fitlm(trendLineYs, trendLineXs,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
rsq = lm.Rsquared.Adjusted;
meanAllY = mean(allYs);
SE = std(allYs) / (sqrt(numel(allYs)));

meanLessThan150y = mean(lessThan150Ys);
SELessThan150y = std(lessThan150Ys) / (sqrt(numel(lessThan150Ys)))

meanGraterThan150y = mean(trendLineYs);
SEGraterThan150y = std(trendLineYs) / (sqrt(numel(trendLineYs)))

% text(-75,600,['R^2=',num2str(rsq), ' p-val=',num2str(pval), ' m=',num2str(gradient), ' \mu all data =',num2str(meanAllY),...
%      ' \mu <150um=',num2str(meanLessThan150y), ' \mu >=150um=',num2str(meanGraterThan150y)]);
text(-75,600,[' \mu all data =',num2str(meanAllY), ' SE = ', num2str(SE) ...
     ' \mu <150um=',num2str(meanLessThan150y), ' \mu >=150um=',num2str(meanGraterThan150y)]);
ylabel('Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northoutside');
ylim([0 800]);
xlim([-100 400]);
title(['Realworld Data - Lamina Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

%plot scatter plot Midvein Cells cell area at divsion vs distance from
%lp-boundary
figure; hold on
allYs = [];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & (ismid==1));
    allYs = [allYs; cellAreaAtDiv(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),cellAreaAtDiv(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
meanAllY = mean(allYs);
SE = std(allYs) / (sqrt(numel(allYs)));
text(-75,600,['\mu=',num2str(meanAllY), ' SE = ', num2str(SE)]);
ylabel('Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northoutside');
ylim([0 800]);
xlim([-100 400]);
title(['Realworld Data - Midvein Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end

function scatter_cellAreaDiv_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray)
cmap = colormap(lines(numel(fnamelist)));
%cmap = colourArray;

%plot scatter plot All Cells cell area at divsion vs growth rate
figure; hold on
cA=[];
kA=[];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find(hasDiv==1);
    cA = [cA;cellAreaAtDiv(ptr)];
    kA = [kA;kArea(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(cellAreaAtDiv(ptr),kArea(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(cA, kA, 1);
yfit = polyval(p,cA);
plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
lm = fitlm(cA, kA,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
text(300,0.09,['R^2=',num2str(lm.Rsquared.Adjusted), ' p-val=',num2str(pval), ' m=',num2str(gradient)]);

xlabel('Cell Area at Division um2')
ylabel('kArea h^-1')
%legend(legendTxt,'location','northoutside');
xlim([0 800]);
ylim([0 0.1]);
title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end

function scatter_WT_cellCycleDur_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray)
cmap = colormap(lines(numel(fnamelist)));
%plot scatter plot All Cells cell-cycle duration vs growth rate
figure; hold on
cA=[];
kA=[];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & ~isnan(cellCycleDuration));
    cA = [cA;cellCycleDuration(ptr)];
    kA = [kA;kArea(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(cellCycleDuration(ptr),kArea(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(cA, kA, 1);
yfit = polyval(p,cA);
plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
lm = fitlm(cA, kA,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
text(20,0.09,['R^2=',num2str(lm.Rsquared.Adjusted), ' p-val=',num2str(pval), ' m=',num2str(gradient)]);

xlabel('Cell Cycle Duration (hours)')
ylabel('kArea h^-1')
%legend(legendTxt,'location','northoutside');
xlim([0 55]);
ylim([0 0.1]);
title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

%plot scatter plot pavement Cells cell-cycle duration vs growth rate
figure; hold on
cA=[];
kA=[];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & ~isnan(cellCycleDuration) & strcmp(identity, 'PC'));
    cA = [cA;cellCycleDuration(ptr)];
    kA = [kA;kArea(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(cellCycleDuration(ptr),kArea(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(cA, kA, 1);
yfit = polyval(p,cA);
plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
lm = fitlm(cA, kA,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
text(20,0.09,['R^2=',num2str(lm.Rsquared.Adjusted), ' p-val=',num2str(pval), ' m=',num2str(gradient)]);

xlabel('Cell Cycle Duration (hours)')
ylabel('kArea h^-1')
%legend(legendTxt,'location','northoutside');
xlim([0 55]);
ylim([0 0.1]);
title(['Realworld Data - Pavement Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');

%plot scatter plot pavement Cells cell-cycle duration vs growth rate
figure; hold on
cA=[];
kA=[];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & ~isnan(cellCycleDuration) & ~strcmp(identity, 'PC'));
    cA = [cA;cellCycleDuration(ptr)];
    kA = [kA;kArea(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(cellCycleDuration(ptr),kArea(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(cA, kA, 1);
yfit = polyval(p,cA);
plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
lm = fitlm(cA, kA,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
text(20,0.09,['R^2=',num2str(lm.Rsquared.Adjusted), ' p-val=',num2str(pval), ' m=',num2str(gradient)]);

xlabel('Cell Cycle Duration (hours)')
ylabel('kArea h^-1')
%legend(legendTxt,'location','northoutside');
xlim([0 55]);
ylim([0 0.1]);
title(['Realworld Data - Non-Pavement Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end

function scatter_cellCycleDur_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray)
cmap = colormap(lines(numel(fnamelist)));

%plot scatter plot All Cells cell cell cycle duration vs growth rate
figure; hold on
cA=[];
kA=[];
for i = numel(fnamelist):-1:1 
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((hasDiv==1) & ~isnan(cellCycleDuration));
    cA = [cA;cellCycleDuration(ptr)];
    kA = [kA;kArea(ptr)];
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(cellCycleDuration(ptr),kArea(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
end
p = polyfit(cA, kA, 1);
yfit = polyval(p,cA);
plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
lm = fitlm(cA, kA,'linear');
t = lm.Coefficients;
pval = t{2,4};
gradient = t{2,1};
text(20,0.09,['R^2=',num2str(lm.Rsquared.Adjusted), ' p-val=',num2str(pval), ' m=',num2str(gradient)]);

xlabel('Cell Cycle Duration (hours)')
ylabel('kArea h^-1')
%legend(legendTxt,'location','northoutside');
xlim([0 55]);
ylim([0 0.1]);
title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end

function scatter_cellArea_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray)
cmap = colourArray;
printTxt = 1;

%plot cell area against LN growth rate of all cells different plots for each timepoint
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    cA = log(cellAreas);
    kA = kArea;
    notnanptr = ~isnan(kA);
    kA = kA(notnanptr);
    cA = cA(notnanptr);
    figure; hold on
    h = plot(cA,kA,'Marker','o','MarkerFaceColor',cmap(i,:), ...
        'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    p = polyfit(cA, kA, 1);
    yfit = polyval(p,cA);
    
    yresid = kA - yfit;
    SSredid = sum(yresid.^2);
    SStotal = (length(kA)-1) * var(kA);
    rsq = 1 - SSredid/SStotal;   
    if rsq >= 0.1
        plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
    end
    lm = fitlm(cA, kA,'linear');
    t = lm.Coefficients;
    pval = t{2,4};
    gradient = t{2,1};
    meanCA = mean(cA);
    SE = std(cA) / (sqrt(numel(cA)));
    meankA = mean(kA);
    SEkA = std(kA) / (sqrt(numel(kA)));
    plot([meanCA,meanCA],[0,0.1],'r','linewidth',3);
    plot([0,8],[meankA,meankA],'r--','linewidth',3);
    if(printTxt)
        text(0.25,0.1,['R^2=',num2str(rsq),' p-val=',num2str(pval),' m=',num2str(gradient)]);
        text(0.25,0.095,[' \mu x=', num2str(meanCA), ' SE = ', num2str(SE), ' SEx1.96 = ', num2str(SE*1.96)]);
        text(0.25,0.09,[' \mu y=', num2str(meankA), ' SE = ', num2str(SEkA), ' SEx1.96 = ', num2str(SEkA*1.96)]);
        xlabel('Ln Cell Area um2')
        ylabel('Karea h-1')        
    else
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    end
    ylim([0 0.1]);
    xlim([0 8]);
    title(['Realworld Data - All Cells','- ',graphTitles,'- ',fnamelist{i}]);
    g = gca;    
    printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'-(T',num2str(i),').png'];
    print(printname,'-dpng','-r300');    
end

end

function scatter_WT_cellArea_growthRate(fileAndLocation,fnamelist,graphTitles,colourArray)
cmap = colourArray;
printTxt = 1;
ylimMax = 0.1; %epi 2700 subepi 900
ylimMin = 0;
xlimMax = 8;
xlimMin = 0;

%plot cell area against LN growth rate of all cells different plots for each timepoint
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    cA = log(cellAreas);
    kA = kArea;
    figure; hold on
    h = plot(cA,kA,'Marker','o','MarkerFaceColor',cmap(i,:), ...
        'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    p = polyfit(cA, kA, 1);
    yfit = polyval(p,cA);
    
    yresid = kA - yfit;
    SSredid = sum(yresid.^2);
    SStotal = (length(kA)-1) * var(kA);
    rsq = 1 - SSredid/SStotal; 
    if rsq >= 0.1
        plot(cA,yfit,'Linewidth',3,'Color',[1 0 1])
    end
    lm = fitlm(cA, kA,'linear');
    t = lm.Coefficients;
    pval = t{2,4};
    gradient = t{2,1};
    meanCA = mean(cA);    
    SE = std(cA) / (sqrt(numel(cA)));
    meankA = mean(kA);
    SEkA = std(kA) / (sqrt(numel(kA)));
    plot([meanCA,meanCA],[0,0.1],'r','linewidth',3);
    plot([0,8],[meankA,meankA],'r--','linewidth',3);
    if(printTxt)
        text(0.25,0.1,['R^2=',num2str(rsq),' p-val=',num2str(pval),' m=',num2str(gradient) ]);
        text(0.25,0.095,['\mu x=', num2str(meanCA), 'SE = ', num2str(SE), ' SEx1.96 = ', num2str(SE*1.96)]);
        text(0.25,0.09,[' \mu y=', num2str(meankA), ' SE = ', num2str(SEkA), ' SEx1.96 = ', num2str(SEkA*1.96)]);
        xlabel('Ln Cell Area um2')
        ylabel('Karea h-1')
    else
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
    end
    ylim([0 0.1]);
    xlim([0 8]);
    title(['Realworld Data - All Cells','- ',graphTitles]);
    g = gca;
    printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'-(T',num2str(i),').png'];
    print(printname,'-dpng','-r300');    
end

end

function hist_cellCycleDur_frequency(fileAndLocation,fnamelist,graphTitles,colourArray,concatTo24h)

%plot hist of cell cycle duration agains frequency for lamina Cells below
%150um from lp-boundry
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    
    
    ptr = find((hasDiv == 1) & (ismid == 0) & ~isnan(cellCycleDuration));
    projections = yDists(ptr);
    cellCycDurOfThatHaveDivided = cellCycleDuration(ptr);
    ptr2 = find(projections <= 150);
    projections = projections(ptr2);
    cellCycDurOfThatHaveDivided = cellCycDurOfThatHaveDivided(ptr2);
    
    
    
    
%     ptr = find((hasDiv == 1) & ~isnan(cellCycleDuration));
%     projections = xDists(ptr);
%     cellCycDurOfThatHaveDivided = cellCycleDuration(ptr);
    
    
    
    mean_cellCycDurThatHaveDivided(i) = mean(cellCycDurOfThatHaveDivided);
    std_cellCycDurThatHaveDivided(i) = std(cellCycDurOfThatHaveDivided);
    store_cellCycDurThatHaveDivided{i} = cellCycDurOfThatHaveDivided;
    bins = 0:2:30;
    for k = 1:numel(bins)-1
        idx = find((cellCycDurOfThatHaveDivided >= bins(k)) &  (cellCycDurOfThatHaveDivided < bins(k+1)));
        meanCellCycDur(k,i) = mean(cellCycDurOfThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
if (concatTo24h)
    idx = 1;
    store_CellCycDurThatHaveDivided_24h = [];
    for j = 1:2:numel(fnamelist)
        numberOfCells_24h(:,idx) = numberOfCells(:,j) + numberOfCells(:,j+1);
        colourArray_24h(idx,:) = colourArray(j,:);
        store_cellCycDurThatHaveDivided_24h{idx} = cat(1,store_cellCycDurThatHaveDivided{j},store_cellCycDurThatHaveDivided{j+1});
        mean_cellCycDurThatHaveDivided_24h(idx) = round(mean(store_cellCycDurThatHaveDivided_24h{idx}),1);
        std_cellCycDurThatHaveDivided_24h(idx) = round(std(store_cellCycDurThatHaveDivided_24h{idx}),1);
        legendTxt24h{idx} = [legendTxt{j},' + ',legendTxt{j+1}];
        idx = idx + 1;
    end
else
    numberOfCells_24h = numberOfCells;
    colourArray_24h = colourArray;
    mean_cellCycDurThatHaveDivided_24h = round(mean_cellCycDurThatHaveDivided,1);
    std_cellCycDurThatHaveDivided_24h = round(std_cellCycDurThatHaveDivided,1);
    legendTxt24h = legendTxt;
end
colormap(colourArray_24h);
b = bar(binCentroids,numberOfCells_24h, 'EdgeColor',[0 0 0]); 

for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end

allMeans = mean(cat(1,store_cellCycDurThatHaveDivided_24h{:}));
allStd = std(cat(1,store_cellCycDurThatHaveDivided_24h{:}));

ylabel('Frequency')
xlabel('Cell Cycle Duration (hours)')
legend(legendTxt24h,'location','northeast','Interpreter', 'none');
title(['Realworld Data - Lamina Cells beloew 150um','- ',graphTitles]);
text(20,3,['mean =',mat2str(mean_cellCycDurThatHaveDivided_24h),' std =',mat2str(std_cellCycDurThatHaveDivided_24h)],'color',[1 0 1]);
text(20,8,['mean of all data =',num2str(allMeans),' std =',num2str(allStd)],'color',[1 0 1]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
%ylim([0 31]);
%xlim([-200 100]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2.5 f.Position(4)];
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
%print(printname,'-dpng','-r300');
end

function hist_cellAreaDiv_frequency(fileAndLocation,fnamelist,graphTitles,colourArray,concatTo24h)

%plot hist of areas at divsion against freq Lam Cells below 150um from
%lp-bound
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find((hasDiv == 1) & (ismid == 0));
    projections = yDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    ptr2 = find(projections <= 150);
    projections = projections(ptr2);
    areaOfCellsThatHaveDivided = areaOfCellsThatHaveDivided(ptr2);
    mean_areaOfCellsThatHaveDivided(i) = mean(areaOfCellsThatHaveDivided);
    std_areaOfCellsThatHaveDivided(i) = std(areaOfCellsThatHaveDivided);
    store_areaOfCellsThatHaveDivided{i} = areaOfCellsThatHaveDivided;
    bins = 0:25:400;
    for k = 1:numel(bins)-1
        idx = find((areaOfCellsThatHaveDivided >= bins(k)) &  (areaOfCellsThatHaveDivided < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
if (concatTo24h)
    idx = 1;
    store_areaOfCellsThatHaveDivided_24h = [];
    for j = 1:2:numel(fnamelist)
        numberOfCells_24h(:,idx) = numberOfCells(:,j) + numberOfCells(:,j+1);
        colourArray_24h(idx,:) = colourArray(j,:);
        store_areaOfCellsThatHaveDivided_24h{idx} = cat(1,store_areaOfCellsThatHaveDivided{j},store_areaOfCellsThatHaveDivided{j+1});
        mean_areaOfCellsThatHaveDivided_24h(idx) = round(mean(store_areaOfCellsThatHaveDivided_24h{idx}),1);
        std_areaOfCellsThatHaveDivided_24h(idx) = round(std(store_areaOfCellsThatHaveDivided_24h{idx}),1);
        idx = idx + 1;
    end
else
    numberOfCells_24h = numberOfCells;
    colourArray_24h = colourArray;
    mean_areaOfCellsThatHaveDivided_24h = round(mean_areaOfCellsThatHaveDivided,1);
    std_areaOfCellsThatHaveDivided_24h = round(std_areaOfCellsThatHaveDivided,1);
end
colormap(colourArray_24h);
b = bar(binCentroids,numberOfCells_24h, 'EdgeColor',[0 0 0]); 

for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end

allMeans = mean(cat(1,store_areaOfCellsThatHaveDivided_24h{:}));
allStd = std(cat(1,store_areaOfCellsThatHaveDivided_24h{:}));

ylabel('Frequency')
xlabel('Cell Area at Division um2')
legend(legendTxt,'location','northeast','Interpreter', 'none');
title(['Realworld Data - Lamina Cells below 150um lp-boundary','- ',graphTitles]);
text(200,6,['mean =',mat2str(mean_areaOfCellsThatHaveDivided_24h),'std =',mat2str(std_areaOfCellsThatHaveDivided_24h)],'color',[1 0 1]);
text(200,8,['mean of all data =',num2str(allMeans),'std =',num2str(allStd)],'color',[1 0 1]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
%ylim([0 21]);
%xlim([-200 100]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2.5 f.Position(4)];


end

function hist_cellAreaDiv_distance(fileAndLocation,fnamelist,graphTitles,colourArray)

%plot hist of cell area at division vs lp-boundary
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find(hasDiv == 1);
    projections = yDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    bins = -100:50:400;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northwest','Interpreter', 'none');
title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 1400]);
xlim([-100 800]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];

%plot hist of lamina cell area vs lp-boundary
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find((hasDiv == 1) & (ismid == 0));
    projections = yDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    bins = -100:50:400;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northwest','Interpreter', 'none');
title(['Realworld Data - Lamina Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 1400]);
xlim([-100 800]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];

%plot hist of midvein cell area vs lp-boundary
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find((hasDiv == 1) & (ismid == 1));
    projections = yDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    bins = -100:50:400;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area at Division um2')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northwest','Interpreter', 'none');
title(['Realworld Data - Midvein Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 1400]);
xlim([-100 800]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];

%plot hist of cell area at division vs distance midvein
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find(hasDiv == 1);
    projections = xDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    bins = -150:25:100;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area at Division um2')
xlabel('Distance from midvein in um')
legend(legendTxt,'location','northwest','Interpreter', 'none');

title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 1400]);
xlim([-200 100]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];

%plot hist of midvein cell area at division vs distance midvein
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find((hasDiv == 1) & (ismid == 1));
    projections = xDists(ptr);
    areaOfCellsThatHaveDivided = cellAreaAtDiv(ptr);
    bins = -150:25:100;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(areaOfCellsThatHaveDivided(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Cell Area at Division um2')
xlabel('Distance from midvein in um')
legend(legendTxt,'location','northwest','Interpreter', 'none');

title(['Realworld Data - Midvein Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 1400]);
xlim([-200 100]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];
end



function scatter_growthRate_distance(fileAndLocation,fnamelist,graphTitles,colourArray)
%cmap = colormap(jet(numel(fnamelist)));
cmap = colourArray;

%plot scatter plot all Cells KPAR vs distance from
%lp-boundary
figure; hold on
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((ismid==1) | (ismid==0));
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(yDists(ptr),kPAR(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
    %uistack(h,'bottom')
end
ylabel('K_{Midline} (% per h)')
xlabel('Distance from lamina-petiole boundary in um')
%legend(fliplr(legendTxt),'location','northeast','Interpreter', 'none');
legend('TL09-10','TL08-09','TL07-08','TL06-07',...
    'TL05-06','TL04-05','TL02-04','location','northeast');
ylim([0 0.04]);
xlim([-100 650]);
title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');


%plot scatter plot all Cells KPER vs distance from
%midvein
figure; hold on
for i = numel(fnamelist):-1:1
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    ptr = find((ismid==1) | (ismid==0));
    if isempty(ptr)
       plot(nan,nan,'Marker','o','MarkerFaceColor',cmap(i,:), ...
           'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    else
        h = plot(xDists(ptr),kPER(ptr),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0]);
    end
    %uistack(h,'bottom')
end
ylabel('K_{PerMidline} (% per h)')
xlabel('Distance from midvein in um')
%legend(fliplr(legendTxt),'location','northeast','Interpreter', 'none');
legend('TL09-10','TL08-09','TL07-08','TL06-07',...
    'TL05-06','TL04-05','TL02-04','location','eastoutside');
ylim([-0.01 0.04]);
xlim([-300 200]);
title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
end

function hist_growthRate_distance(fileAndLocation,fnamelist,graphTitles,colourArray)
%plot hist of KPAR vs lp-boundary for all cells
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find((ismid==0) | (ismid==1));
    projections = yDists(ptr);
    cellAreasThatMeetSpec = kPAR(ptr);
    bins = -100:50:700;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(cellAreasThatMeetSpec(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; colormap(jet(numel(fnamelist)))
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Kpar (% per h)')
xlabel('Distance from lamina-petiole boundary in um')
legend(legendTxt,'location','northeast','Interpreter', 'none');

title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 0.04]);
xlim([-100 650]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];

%plot hist of KPER vs midvein for all cells
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
    legendTxt{i} = [c{2},'-',c{3}];
    [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
        cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);    
    ptr = find((ismid==0) | (ismid==1));
    projections = xDists(ptr);
    cellAreasThatMeetSpec = kPER(ptr);
    bins = -250:25:200;
    for k = 1:numel(bins)-1
        idx = find((projections >= bins(k)) &  (projections < bins(k+1)));
        meanCellArea(k,i) = mean(cellAreasThatMeetSpec(idx),'omitnan');
        numberOfCells(k,i) = numel(idx);
        binCentroids(k) = (bins(k) + bins(k+1)) / 2;
        binRanges{k} = [num2str(bins(k)), '-' ,num2str(bins(k+1)-1)];
    end
end
figure; %colormap(jet(numel(fnamelist)))
colormap(colourArray);
b = bar(binCentroids,meanCellArea, 'EdgeColor',[0 0 0]); 
for k = 1:numel(b)
    b(k).EdgeColor = 'black';
end
ylabel('Average Kper (% per h)')
xlabel('Distance from midvein in um')
legend(legendTxt,'location','northeast','Interpreter', 'none');

title(['Realworld Data - All Cells','- ',graphTitles]);
g = gca;
g.XTick = round(binCentroids,2);
g.XTickLabel = binRanges;
g.FontSize = 7;
ylim([0 0.04]);
xlim([-300 200]);
f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];
end

function normalBoxPlot(fileAndLocation,fnamelist,graphTitles,colourArray)
%cmap = colormap(jet(numel(fnamelist)));
cmap = colourArray;
figure; hold on;
startPlotX = [1 2 3 4];
smallPlotOffset = 0.1;
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');

[xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    
    % comp and lam    
    ptra = ~ismid & iscomp;
    ptrb = ~ismid & hasDiv & ~iscomp;
    ptr = ptra | ptrb;    
    data1{i} = kArea(ptr);
    
       
    % noncomp and lam    
    ptr = ~ismid & ~iscomp & ~hasDiv;      
    data2{i} = kArea(ptr);
            
    % comp and mid    
    ptra = ismid & iscomp;
    ptrb = ismid & hasDiv & ~iscomp;
    ptr = ptra | ptrb;     
    data3{i} = kArea(ptr);
    
    % noncomp and mid
    ptr = ismid & ~iscomp & ~hasDiv;
    data4{i} = kArea(ptr);    
    
end

alldata = [data1 [nan] data2 [nan] data3 [nan] data4];
%remove data with less than 15 sample points
for j = 1:numel(alldata)
    if(~isnan(alldata{j}))
        if(numel(alldata{j}) < 15)
            alldata{j} = double.empty(0);
        end
    end
end

alldataEXP = [];
grp = [];
for kk = 1:numel(alldata)
    if(~isempty(alldata{kk}))
        dataToAdd = alldata{kk}';
    else
        dataToAdd = 999;
    end
    alldataEXP = [alldataEXP dataToAdd];
    if(~isempty(alldata{kk}))
        grp = [grp ones(1,numel(alldata{kk})) .* kk];
    else
        grp = [grp kk];
    end
end
figure;
b = boxplot(alldataEXP,grp,'Notch','on');
h = findobj(gca,'Tag','Box');

allCmap = repmat(cmap,4,1);
kkk = 1;
for i = numel(h):-1:1
    if (~isnan(h(i).YData(1)))
        h(i).Color = allCmap(kkk,:);
        kkk = kkk + 1;
    end    
end

ylim([0 0.08]);

f = gcf;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];
ylabel('Karea (h^{-1})','Interpreter', 'tex');
title(graphTitles);
end

function fancyBoxPlot_kAREA_midAndLam(fileAndLocation,fnamelist,graphTitles,colourArray)
%cmap = colormap(jet(numel(fnamelist)));
cmap = colourArray;
figure; hold on;
startPlotX = [1 2 3 4];
smallPlotOffset = 0.1;
for i = 1:numel(fnamelist)
    c = strsplit(fnamelist{i},',');
%     [xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
%     cellAreaAtDiv,kArea,kPAR,kPER] = Florent_getData(c{1},c{2},c{3},fileAndLocation);


[xDists,yDists,cellAreas,ismid,iscomp,hasDiv,...
    cellAreaAtDiv,kArea,kPAR,kPER, identity, cellCycleDuration] = Florent_getData(c{1},c{2},c{3},fileAndLocation);
    
    % comp and lam    
    ptra = ~ismid & iscomp;
    ptrb = ~ismid & hasDiv & ~iscomp;
    ptr = ptra | ptrb;    
    data1 = kArea(ptr);
    
    err = std(data1) / sqrt(numel(data1));
    if ~isempty(data1)
        errorbar(startPlotX(1)+smallPlotOffset*i,nanmean(data1),err,'CapSize',18,'color',cmap(i,:),'Linewidth',2);
        %plot(startPlotX(1)+smallPlotOffset*i,data1,'.','color',cmap(i,:));    
        plot(startPlotX(1)+smallPlotOffset*i,nanmean(data1),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0],'MarkerSize',5);
        
    end
    
    % noncomp and lam    
    ptr = ~ismid & ~iscomp & ~hasDiv;      
    data2 = kArea(ptr);
    
    err = std(data2) / sqrt(numel(data2));
    if ~isempty(data2)
        errorbar(startPlotX(2)+smallPlotOffset*i,nanmean(data2),err,'color',cmap(i,:),'Linewidth',2);
        %plot(startPlotX(2)+smallPlotOffset*i,data2,'.','color',cmap(i,:));    
        plot(startPlotX(2)+smallPlotOffset*i,nanmean(data2),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0],'MarkerSize',5);        
    end 
        
    % comp and mid    
    ptra = ismid & iscomp;
    ptrb = ismid & hasDiv & ~iscomp;
    ptr = ptra | ptrb;     
    data3 = kArea(ptr);
    
    err = std(data3) / sqrt(numel(data3));
    if ~isempty(data3)
        errorbar(startPlotX(3)+smallPlotOffset*i,nanmean(data3),err,'color',cmap(i,:),'Linewidth',2);
        %plot(startPlotX(3)+smallPlotOffset*i,data3,'.','color',cmap(i,:));    
        plot(startPlotX(3)+smallPlotOffset*i,nanmean(data3),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0],'MarkerSize',5);        
    end
        
    % noncomp and mid
    ptr = ismid & ~iscomp & ~hasDiv;
    data4 = kArea(ptr);
    
    err = std(data4) / sqrt(numel(data4));
    if ~isempty(data4)
        errorbar(startPlotX(4)+smallPlotOffset*i,nanmean(data4),err,'color',cmap(i,:),'Linewidth',2);
        %plot(startPlotX(4)+smallPlotOffset*i,data4,'.','color',cmap(i,:));    
        plot(startPlotX(4)+smallPlotOffset*i,nanmean(data4),'Marker','o','MarkerFaceColor',cmap(i,:), ...
            'LineStyle','none','MarkerEdgeColor',[0 0 0],'MarkerSize',5);        
    end    
end
% 
% [h,icons,plots,str] = legend('TL02-04', 'TL04-05', 'TL05-06', ...
%     'TL06-07', 'TL07-08', 'TL08-09', 'TL09-10', ...
%     'location','eastoutside');
for i = 1:numel(fnamelist)
    icons(numel(fnamelist)+i*2).Color = cmap(i,:);
end
ylabel('Karea (% per h)');
title(graphTitles);
g = gca;
% xlim([100 500]);
% ylim([0 2700]);
g.XTick = 1:0.1:5;
for i = 1:numel(g.XTickLabel)
    g.XTickLabel{i} = '';
end
g.XTickLabel{4} = 'Lam Comp';
g.XTickLabel{14} = 'Lam Arrest';
g.XTickLabel{24} = 'Mid Comp';
g.XTickLabel{34} = 'Mid Arrest';

printname = [g.Title.String,'-',g.XLabel.String,'-',g.YLabel.String,'.png'];
print(printname,'-dpng','-r300');
%g.YScale = 'log';
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

end





















