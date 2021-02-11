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
    
    %% Get contents of the directory where the parc data for this subject are stored.
    ashs = dir(fullfile(sub_ashs(s).folder, sub_ashs(s).name, 'dt-neuro-parcellation-volume*/parc.nii.gz'));
    
    % Remove the '.' and '..' files.
    ashs = ashs(arrayfun(@(x) x.name(1), ashs) ~= '.');
    
    % Read in nifti containing parcellation.
    A = niftiRead(fullfile(ashs(1).folder, ashs(1).name));
    
    % Get labels for different codes in the ASHS nifti.
    A_labels = struct2table(jsondecode(fileread(fullfile(ashs(1).folder, 'label.json'))));
    
    %% Get contents of the directory where the parc REPEAT data for this subject are stored.
    hand = dir(fullfile(sub_hand(s).folder, sub_hand(s).name, 'dt-neuro-parcellation-volume*/parc.nii.gz'));
    
    % Remove the '.' and '..' files.
    hand = hand(arrayfun(@(x) x.name(1), hand) ~= '.');
    
    % Read in nifti containing REPEAT parcellation.
    B = niftiRead(fullfile(hand(1).folder, hand(1).name));
    
    % Get labels for different codes in the ASHS nifti.
    B_labels = struct2table(jsondecode(fileread(fullfile(hand(1).folder, 'label.json'))));
    
    %% Calculate the generalized dice similarity coefficient (gdsc) for this subject.
    
    % Ensure that the segmentation is binary. Select only those that should
    % be included in the gdsc calculation (Yushkevic et al., Human Brain
    % Mapping, 2015).
    % Find the row indices for indicated subfields.
    idx = find(contains(A_labels.name, subfield_list));
    subfield_codes = str2double(A_labels.label(idx));
    
    % Binarize.
    A1.data = ismember(A.data, subfield_codes);
    B1.data = ismember(B.data, subfield_codes);
    
    % Calculate DICE coefficient.
    common = (A1.data & B1.data);
    a = sum(common(:));
    b = sum(A1.data(:));
    c = sum(B1.data(:));
    gdsc(s) = 2*(a/(b+c));
    
    %% Calculate "dice similarity coefficient" (dsc) for each subfield.
    
    for subfield = 1:length(subfield_list)
        
        % Find the row index for this subfield.
        idx = find(strcmp(A_labels.name, subfield_list{subfield}));
        A_subfield_code = str2num(A_labels.label{idx});
        
        % Find the row index for this subfield.
        idx = find(strcmp(B_labels.name, subfield_list{subfield}));
        B_subfield_code = str2num(B_labels.label{idx});
        
        % Select only the voxels that correspond to the current ROI and make binary.
        A2.data = (A.data == A_subfield_code);
        B2.data = (B.data == B_subfield_code);
        
        % Calculate DICE coefficient for this ROI for this subject.
        common = (A2.data & B2.data);
        a = sum(common(:));
        b = sum(A2.data(:));
        c = sum(B2.data(:));
        dsc(s, subfield) = 2*(a/(b+c));
        
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
scatter(ones([s 1]), dsc(:, 1), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(2*ones([s 1]), dsc(:, 2), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(3*ones([s 1]), dsc(:, 3), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(4*ones([s 1]), dsc(:, 4), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(5*ones([s 1]), dsc(:, 5), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(6*ones([s 1]), dsc(:, 6), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(7*ones([s 1]), dsc(:, 7), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(8*ones([s 1]), dsc(:, 8), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

scatter(idx_gdsc*ones([s 1]), gdsc', 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', gdsc_color, 'MarkerEdgeColor', gdsc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% % Means and standard deviations for each subfield.
% x_mean = cat(2, mean(dsc, 1), mean(gdsc));
% x_sd = cat(2, std(dsc, 1), std(gdsc));
% xval = linspace(1, length(x_mean), length(x_mean));
% scatter(xval, x_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
% errorbar(xval, x_mean, x_sd, 'Color', 'k', 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xlim_lo = 0.5; xlim_hi = length(subfield_list) + 1 + 0.5;
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = [subfield_list 'GDSC']
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;
xax.TickLabelRotation = labelrotation;

% yaxis
ylimlo = 0.40; ylimhi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off
% set(a,'color',[220 220 220]/255);
% set(gcf,'color',[220 220 220]/255);
% set(gcf, 'InvertHardCopy', 'off');
a.YLabel.String = 'Dice Similarity Coefficient';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', 'plot_scatter_validation_dsc'), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', 'plot_scatter_validation_dsc'), '-depsc')
    
end

hold off;

% Set up output table of mean/stdev/min/max of dsc.
t = [mean(dsc(:, 1)) std(dsc(:, 1)) min(dsc(:, 1)) max(dsc(:, 1));
    mean(dsc(:, 2)) std(dsc(:, 2)) min(dsc(:, 2)) max(dsc(:, 2));
    mean(dsc(:, 3)) std(dsc(:, 3)) min(dsc(:, 3)) max(dsc(:, 3));
    mean(dsc(:, 4)) std(dsc(:, 4)) min(dsc(:, 4)) max(dsc(:, 4));
    mean(dsc(:, 5)) std(dsc(:, 5)) min(dsc(:, 5)) max(dsc(:, 5));
    mean(dsc(:, 6)) std(dsc(:, 6)) min(dsc(:, 6)) max(dsc(:, 6));
    mean(dsc(:, 7)) std(dsc(:, 7)) min(dsc(:, 7)) max(dsc(:, 7));
    mean(dsc(:, 8)) std(dsc(:, 8)) min(dsc(:, 8)) max(dsc(:, 8));
    mean(gdsc) std(gdsc) min(gdsc) max(gdsc)];

t_out = cat(2, [subfield_list 'gdsc']', array2table(t));
t_out.Properties.VariableNames = {'subfield', 'mean', 'stdev', 'min', 'max'};
writetable(t_out, fullfile(rootDir, 'supportFiles', 'dsc_validation_upenn.csv'));



