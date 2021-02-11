% This script will read in two ROIs and calculate the DICE coefficient. The
% Dice coefficient is a measurement of spatial overlap. Here is a briefing:
% http://sve.loni.ucla.edu/instructions/metrics/dice/

% DICE = 2 * (A intersection B) / (numel(A) + numel(B))

% It is very important that the two ROIs are co-registered and in the same
% space as one another. Note that this script does not check for that.

clear all; close all; clc;

addpath(genpath('/Applications/jsonlab-2.0/'));

save_figures = 'yes';

rootDir = '/Volumes/240/ashs-validation/';
data = 'upenn';

if strcmp(data, 'upenn')
    
    % Identify where parcellations are.
    blprojectid = 'proj-601324fdc081031288ba4792'; %upenn
    % Identify where repeat parcellations are.
    blprojectid_repeat = 'proj-601324fdc081031288ba4792-repeat'; %upenn
    
elseif strcmp(data, 'devti')
    
    % Identify where parcellations are.
    blprojectid = 'proj-5fbd00ce7e8ecbb59eabfe3b'; %devti project
    % Identify where REPEAT parcellations are.
    blprojectid_repeat = 'proj-5fbd00ce7e8ecbb59eabfe3b-repeat'; %devti project
    
    % Import age data.
    beh = readtable([rootDir 'supportFiles/participant_information.csv'], 'TreatAsEmpty', {'.', 'na'});
    
end

subfield_list = {'CA1', 'CA2', 'DG', 'CA3', 'SUB', 'ERC', 'BA35', 'BA36'};

% Get contents of the directory where the tract measures for this subject are stored.
sub = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
sub = sub(arrayfun(@(x) x.name(1), sub) ~= '.');

% Keep only names that are subject folders.
sub = sub(arrayfun(@(x) x.name(1), sub) == 's');

% Get contents of the directory where the tract measures for this subject are stored.
sub_repeat = dir(fullfile(rootDir, blprojectid_repeat));

% Remove the '.' and '..' files.
sub_repeat = sub_repeat(arrayfun(@(x) x.name(1), sub_repeat) ~= '.');

% Keep only names that are subject folders.
sub_repeat = sub_repeat(arrayfun(@(x) x.name(1), sub_repeat) == 's');

for s = 1:size(sub, 1)
    
    %% Get contents of the directory where the volume data for this subject are stored.
    parc = dir(fullfile(sub(s).folder, sub(s).name, 'dt-neuro-parc-stats.tag-volume*/parc_nodes.csv'));
    
    % Remove the '.' and '..' files.
    parc = parc(arrayfun(@(x) x.name(1), parc) ~= '.');
    
    % Read in nifti containing parcellation.
    A = readtable(fullfile(parc(1).folder, parc(1).name));
    
    %% Get contents of the directory where the volume REPEAT data for this subject are stored.
    parc_repeat = dir(fullfile(sub_repeat(s).folder, sub_repeat(s).name, 'dt-neuro-parc-stats.tag-volume*/parc_nodes.csv'));
    
    % Remove the '.' and '..' files.
    parc_repeat = parc_repeat(arrayfun(@(x) x.name(1), parc_repeat) ~= '.');
    
    % Read in nifti containing REPEAT parcellation.
    B = readtable(fullfile(parc(1).folder, parc(1).name));
    
    if strcmp(data, 'devti')
        
        % Get the age-group for this participant: child = 1, adolescent = 2, adult = 3.
        subID = str2double(sub(s).name(end-2:end));
        disp(subID)
        if beh.child(find(beh.participant == subID)) == 1
            age(s) = 1;
        elseif beh.adolescent(find(beh.participant == subID)) == 1
            age(s) = 2;
        elseif beh.adult(find(beh.participant == subID)) == 1
            age(s) = 3;
        end
        
    end
    
    %% Get volume ofeach subfield.
    
    for subfield = 1:length(subfield_list)
        
        % Find the row index for this subfield in the ashs segmentation..
        idx = find(strcmp(A.structureID, subfield_list{subfield}));
        
        % Grab the volume of this subfield for this subject in the ashs segmentation.
        vol(s, subfield) = A.volume(idx);
        
        % Find the row index for this subfield in the ashs repeat segmentation..
        idx = find(strcmp(B.structureID, subfield_list{subfield}));
        
        % Grab the volume of this subfield for this subject in the ashs repeat segmentation.
        repeat_vol(s, subfield) = B.volume(idx);
        
        clear A2 B2
    end
    
    clear A B A1 B2
    
end

figure(1)
hold on;
% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 300;
xtickvalues = 1:length(subfield_list)+1;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
labelrotation = 45;
yticklength = 0;
xticklength = 0.05;

%old colors
% ca1_color = [0 0.4470 0.7410]; %blue
% ca2_color = [0.8500 0.3250 0.0980]; %orange
% dg_color = [0.4940 0.1840 0.5560]; %purple
% ca3_color = [0.9290 0.6940 0.1250]; %yellow
% sub_color = [162 129 0]/255; %gold
% erc_color = [0.4660 0.6740 0.1880]; %green
% ba35_color = [0.3010 0.7450 0.9330]; %cyan
% ba36_color = [0.6350 0.0780 0.1840]; %red
% gdsc_color = [0 0 0]; %black

%ashs segmentation colors
ca1_color = [255 0 0]/255;
ca2_color = [0 255 0]/255;
dg_color = [128 64 255]/255;
ca3_color = [255 255 0]/255;
sub_color = [240 86 224]/255;
erc_color = [210 180 140]/255;
ba35_color = [55 180 255]/255;
ba36_color = [0 0 255]/255;
gdsc_color = [0 0 0]/255; %black

coloralpha = .4;
marker_adult = 'o';
marker_adolescent = 'square';
marker_child = '^';

% Individual data points for subjects and subfields.
if strcmp(data, 'upenn')
    scatter(vol(:, 1), repeat_vol(:, 1), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 2), repeat_vol(:, 2), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 3), repeat_vol(:, 3), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 4), repeat_vol(:, 4), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 5), repeat_vol(:, 5), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 6), repeat_vol(:, 6), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 7), repeat_vol(:, 7), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(:, 8), repeat_vol(:, 8), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
elseif strcmp(data, 'devti')
    markersize = 200;
    %adult = 3
    scatter(vol(find(age==3), 1), repeat_vol(find(age==3), 1), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 2), repeat_vol(find(age==3), 2), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 3), repeat_vol(find(age==3), 3), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 4), repeat_vol(find(age==3), 4), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 5), repeat_vol(find(age==3), 5), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 6), repeat_vol(find(age==3), 6), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 7), repeat_vol(find(age==3), 7), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==3), 8), repeat_vol(find(age==3), 8), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    %adolescent = 2
    scatter(vol(find(age==2), 1), repeat_vol(find(age==2), 1), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 2), repeat_vol(find(age==2), 2), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 3), repeat_vol(find(age==2), 3), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 4), repeat_vol(find(age==2), 4), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 5), repeat_vol(find(age==2), 5), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 6), repeat_vol(find(age==2), 6), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 7), repeat_vol(find(age==2), 7), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==2), 8), repeat_vol(find(age==2), 8), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    %child = 1
    scatter(vol(find(age==1), 1), repeat_vol(find(age==1), 1), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 2), repeat_vol(find(age==1), 2), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 3), repeat_vol(find(age==1), 3), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 4), repeat_vol(find(age==1), 4), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 5), repeat_vol(find(age==1), 5), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 6), repeat_vol(find(age==1), 6), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 7), repeat_vol(find(age==1), 7), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(vol(find(age==1), 8), repeat_vol(find(age==1), 8), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    
    
end


% % Means and standard deviations for each subfield.
% x_mean = cat(2, mean(dsc, 1), mean(gdsc));
% x_sd = cat(2, std(dsc, 1), std(gdsc));
% xval = linspace(1, length(x_mean), length(x_mean));
% scatter(xval, x_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
% errorbar(xval, x_mean, x_sd, 'Color', 'k', 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xlimlo = 0; xlimhi = 5000;
xax = get(gca, 'xaxis');
xax.Limits = [xlimlo xlimhi];
xax.TickValues = [xlimlo (xlimlo+xlimhi)/2 xlimhi];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.TickLabels = {num2str(xlimlo, '%2.0f'), num2str((xlimlo+xlimhi)/2, '%2.0f'), num2str(xlimhi, '%2.0f')};
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
ylimlo = 0; ylimhi = 5000;
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.0f'), num2str((ylimlo+ylimhi)/2, '%2.0f'), num2str(ylimhi, '%2.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off
% set(a,'color',[220 220 220]/255);
% set(gcf,'color',[220 220 220]/255);
% set(gcf, 'InvertHardCopy', 'off');
a.YLabel.String = {'Volume,'; 'Repeat ASHS-segmentation (mm^{3})'};
a.YLabel.FontSize = 19;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Volume,'; 'ASHS-segmentation (mm^{3})'};
a.XLabel.FontSize = 19;
a.XLabel.FontAngle = 'normal';
pbaspect([1 1 1])

legend(subfield_list, 'Location', 'southeast')
legend box off
%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', ['plot_scatter_testretest_volume_' data]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_scatter_testretest_volume_' data]), '-depsc')
    
end

hold off;

% % Set up output table of mean/stdev/min/max of dsc.
% t = [mean(dsc(:, 1)) std(dsc(:, 1)) min(dsc(:, 1)) max(dsc(:, 1));
%     mean(dsc(:, 2)) std(dsc(:, 2)) min(dsc(:, 2)) max(dsc(:, 2));
%     mean(dsc(:, 3)) std(dsc(:, 3)) min(dsc(:, 3)) max(dsc(:, 3));
%     mean(dsc(:, 4)) std(dsc(:, 4)) min(dsc(:, 4)) max(dsc(:, 4));
%     mean(dsc(:, 5)) std(dsc(:, 5)) min(dsc(:, 5)) max(dsc(:, 5));
%     mean(dsc(:, 6)) std(dsc(:, 6)) min(dsc(:, 6)) max(dsc(:, 6));
%     mean(dsc(:, 7)) std(dsc(:, 7)) min(dsc(:, 7)) max(dsc(:, 7));
%     mean(dsc(:, 8)) std(dsc(:, 8)) min(dsc(:, 8)) max(dsc(:, 8));
%     mean(gdsc) std(gdsc) min(gdsc) max(gdsc)];
%
% t_out = cat(2, [subfield_list 'gdsc']', array2table(t));
% t_out.Properties.VariableNames = {'subfield', 'mean', 'stdev', 'min', 'max'};
% writetable(t_out, fullfile(rootDir, 'supportFiles', 'dsc_validation_upenn.csv'));



