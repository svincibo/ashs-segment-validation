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

% Identify where parcellations are.
blprojectid_ashs = 'proj-601324fdc081031288ba4792'; %upenn
% Identify where Hand-drawn parcellations are.
blprojectid_hand = 'proj-601324fdc081031288ba4792-handdrawn'; %upenn
subfield_list = {'CA1', 'CA2', 'DG', 'CA3', 'SUB', 'ERC', 'BA35', 'BA36'};

% Get contents of the directory where the tract measures for this subject are stored.
sub_ashs = dir(fullfile(rootDir, blprojectid_ashs));

% Remove the '.' and '..' files.
sub_ashs = sub_ashs(arrayfun(@(x) x.name(1), sub_ashs) ~= '.');

% Keep only names that are subject folders.
sub_ashs = sub_ashs(arrayfun(@(x) x.name(1), sub_ashs) == 's');

% Get contents of the directory where the tract measures for this subject are stored.
sub_hand = dir(fullfile(rootDir, blprojectid_hand));

% Remove the '.' and '..' files.
sub_hand = sub_hand(arrayfun(@(x) x.name(1), sub_hand) ~= '.');

% Keep only names that are subject folders.
sub_hand = sub_hand(arrayfun(@(x) x.name(1), sub_hand) == 's');

for s = 1:size(sub_ashs, 1)
    
%    %% Get contents of the directory where the intracranial volume (ICV) data for this subject is stored.
%     icvdir = dir(fullfile(sub_ashs(s).folder, sub_ashs(s).name, '*parc-stats-freesurfer*/whole_brain.csv'));
%     
%     % Remove the '.' and '..' files.
%     icvdir = icvdir(arrayfun(@(x) x.name(1), icvdir) ~= '.');
%     
%     % Read in nifti containing parcellation.
%     icvtable = readtable(fullfile(icvdir(1).folder, icvdir(1).name));
%     
%     icv = icvtable.Total_Intracranial_volume;
%     
    %% Get contents of the directory where the volume data for this subject are stored.
    ashs = dir(fullfile(sub_ashs(s).folder, sub_ashs(s).name, 'dt-neuro-parc-stats.tag-volume*/parc_nodes.csv'));
    
    % Remove the '.' and '..' files.
    ashs = ashs(arrayfun(@(x) x.name(1), ashs) ~= '.');
    
    % Read in nifti containing parcellation.
    A = readtable(fullfile(ashs(1).folder, ashs(1).name));
    
    %% Get contents of the directory where the volume REPEAT data for this subject are stored.
    hand = dir(fullfile(sub_hand(s).folder, sub_hand(s).name, 'dt-neuro-parc-stats.tag-volume*/parc_nodes.csv'));
    
    % Remove the '.' and '..' files.
    hand = hand(arrayfun(@(x) x.name(1), hand) ~= '.');
    
    % Read in nifti containing REPEAT parcellation.
    B = readtable(fullfile(ashs(1).folder, ashs(1).name));
        
    %% Calculate "dice similarity coefficient" (dsc) for each subfield.
    
    for subfield = 1:length(subfield_list)
        
        % Find the row index for this subfield in the ashs segmentation..
        idx = find(strcmp(A.structureID, subfield_list{subfield}));
        
        % Grab the volume of this subfield for this subject in the ashs segmentation.
        ashs_vol(s, subfield) = A.volume(idx);
        
        % Find the row index for this subfield in the hand segmentation..
        idx = find(strcmp(B.structureID, subfield_list{subfield}));
        
        % Grab the volume of this subfield for this subject in the hand segmentation.
        hand_vol(s, subfield) = B.volume(idx);        
        
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

% Individual data points for subjects and subfields.
idx_gdsc = length(subfield_list) + 1;
scatter(ashs_vol(:, 1), hand_vol(:, 1), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 2), hand_vol(:, 2), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 3), hand_vol(:, 3), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 4), hand_vol(:, 4), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 5), hand_vol(:, 5), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 6), hand_vol(:, 6), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 7), hand_vol(:, 7), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(ashs_vol(:, 8), hand_vol(:, 8), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

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
a.YLabel.String = {'Volume,'; 'Hand-segmentation (mm^{3})'};
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
    
    print(fullfile(rootDir, 'plots', 'plot_scatter_validation_volume'), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', 'plot_scatter_validation_volume'), '-depsc')
    
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



