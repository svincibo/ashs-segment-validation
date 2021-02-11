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
    % Identify where REPEAT parcellations are.
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
    
    %% Get contents of the directory where the parc data for this subject are stored.
    parc = dir(fullfile(sub(s).folder, sub(s).name, 'dt-neuro-parcellation-volume*/parc.nii.gz'));
    
    % Remove the '.' and '..' files.
    parc = parc(arrayfun(@(x) x.name(1), parc) ~= '.');
    
    % Read in nifti containing parcellation.
    A = niftiRead(fullfile(parc(1).folder, parc(1).name));
    
    % Get labels for different codes in the ASHS nifti.
    A_labels = struct2table(jsondecode(fileread(fullfile(parc(1).folder, 'label.json'))));
    
    %% Get contents of the directory where the parc REPEAT data for this subject are stored.
    parc_repeat = dir(fullfile(sub_repeat(s).folder, sub_repeat(s).name, 'dt-neuro-parcellation-volume*/parc.nii.gz'));
    
    % Remove the '.' and '..' files.
    parc_repeat = parc_repeat(arrayfun(@(x) x.name(1), parc_repeat) ~= '.');
    
    % Read in nifti containing REPEAT parcellation.
    B = niftiRead(fullfile(parc_repeat(1).folder, parc_repeat(1).name));
    
    % Get labels for different codes in the ASHS nifti.
    B_labels = struct2table(jsondecode(fileread(fullfile(parc_repeat(1).folder, 'label.json'))));
    
    %% Calculate the generalized dice similarity coefficient (gdsc) for this subject.
    
    % Ensure that the segmentation is binary. Select only those that should
    % be included in the gdsc calculation (Yushkevic et al., Human Brain
    % Mapping, 2015).
    if strcmp(data, 'upenn')
        % Find the row indices for indicated subfields.
        idx = find(contains(A_labels.name, subfield_list));
        subfield_codes = str2double(A_labels.label(idx));
        
        % Binarize.
        A1.data = ismember(A.data, subfield_codes);
        B1.data = ismember(B.data, subfield_codes);
        
    elseif strcmp(data, 'devti')
        
        % Binarize.
        A1.data = (A.data ~= 0);
        B1.data = (B.data ~= 0);
        
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
marker_adolescent = 'square';
marker_child = '^';

% Individual data points for subjects and subfields.
idx_gdsc = length(subfield_list) + 1;

if strcmp(data, 'upenn')
    scatter(ones([s 1]), dsc(:, 1), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(2*ones([s 1]), dsc(:, 2), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(3*ones([s 1]), dsc(:, 3), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(4*ones([s 1]), dsc(:, 4), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(5*ones([s 1]), dsc(:, 5), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(6*ones([s 1]), dsc(:, 6), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(7*ones([s 1]), dsc(:, 7), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(8*ones([s 1]), dsc(:, 8), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(idx_gdsc*ones([s 1]), gdsc', 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', gdsc_color, 'MarkerEdgeColor', gdsc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
elseif strcmp(data, 'devti')
    gdsc=gdsc';
    markersize = 200;
    % adults = 3
    count = length(find(age==3));
    scatter(ones([count 1])+.25, dsc(find(age==3), 1), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(2*ones([count 1])+.25, dsc(find(age==3), 2), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(3*ones([count 1])+.25, dsc(find(age==3), 3), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(4*ones([count 1])+.25, dsc(find(age==3), 4), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(5*ones([count 1])+.25, dsc(find(age==3), 5), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(6*ones([count 1])+.25, dsc(find(age==3), 6), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(7*ones([count 1])+.25, dsc(find(age==3), 7), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(8*ones([count 1])+.25, dsc(find(age==3), 8), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(idx_gdsc*ones([count 1])+.25, gdsc(find(age==3)), 'Marker', marker_adult, 'SizeData', markersize/2, 'MarkerFaceColor', gdsc_color, 'MarkerEdgeColor', gdsc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    % adolescents = 2
    count = length(find(age==2));
    scatter(ones([count 1]), dsc(find(age==2), 1), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(2*ones([count 1]), dsc(find(age==2), 2), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(3*ones([count 1]), dsc(find(age==2), 3), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(4*ones([count 1]), dsc(find(age==2), 4), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(5*ones([count 1]), dsc(find(age==2), 5), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(6*ones([count 1]), dsc(find(age==2), 6), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(7*ones([count 1]), dsc(find(age==2), 7), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(8*ones([count 1]), dsc(find(age==2), 8), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(idx_gdsc*ones([count 1]), gdsc(find(age==2)), 'Marker', marker_adolescent, 'SizeData', markersize/2, 'MarkerFaceColor', gdsc_color, 'MarkerEdgeColor', gdsc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    % children = 1
    count = length(find(age==1));
    scatter(ones([count 1])-.25, dsc(find(age==1), 1), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ca1_color, 'MarkerEdgeColor', ca1_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(2*ones([count 1])-.25, dsc(find(age==1), 2), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ca2_color, 'MarkerEdgeColor', ca2_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(3*ones([count 1])-.25, dsc(find(age==1), 3), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', dg_color, 'MarkerEdgeColor', dg_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(4*ones([count 1])-.25, dsc(find(age==1), 4), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ca3_color, 'MarkerEdgeColor', ca3_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(5*ones([count 1])-.25, dsc(find(age==1), 5), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', sub_color, 'MarkerEdgeColor', sub_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(6*ones([count 1])-.25, dsc(find(age==1), 6), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', erc_color, 'MarkerEdgeColor', erc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(7*ones([count 1])-.25, dsc(find(age==1), 7), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ba35_color, 'MarkerEdgeColor', ba35_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(8*ones([count 1])-.25, dsc(find(age==1), 8), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', ba36_color, 'MarkerEdgeColor', ba36_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
    scatter(idx_gdsc*ones([count 1])-.25, gdsc(find(age==1)), 'Marker', marker_child, 'SizeData', markersize/2, 'MarkerFaceColor', gdsc_color, 'MarkerEdgeColor', gdsc_color, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
end

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
ylimlo = 0.80; ylimhi = 1;
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
    
    print(fullfile(rootDir, 'plots', ['plot_scatter_dsc_' data]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_scatter_dsc_' data]), '-depsc')
    
end

hold off;

% Set up output table of mean/stdev/min/max of dsc.
% if strcmp(data, 'upenn')
    t = [mean(dsc(:, 1)) std(dsc(:, 1)) min(dsc(:, 1)) max(dsc(:, 1));
        mean(dsc(:, 2)) std(dsc(:, 2)) min(dsc(:, 2)) max(dsc(:, 2));
        mean(dsc(:, 3)) std(dsc(:, 3)) min(dsc(:, 3)) max(dsc(:, 3));
        mean(dsc(:, 4)) std(dsc(:, 4)) min(dsc(:, 4)) max(dsc(:, 4));
        mean(dsc(:, 5)) std(dsc(:, 5)) min(dsc(:, 5)) max(dsc(:, 5));
        mean(dsc(:, 6)) std(dsc(:, 6)) min(dsc(:, 6)) max(dsc(:, 6));
        mean(dsc(:, 7)) std(dsc(:, 7)) min(dsc(:, 7)) max(dsc(:, 7));
        mean(dsc(:, 8)) std(dsc(:, 8)) min(dsc(:, 8)) max(dsc(:, 8));
        mean(gdsc) std(gdsc) min(gdsc) max(gdsc)];
% elseif strcmp(data, 'devti')
%     t = [mean(dsc(:, 1)) std(dsc(:, 1)) min(dsc(:, 1)) max(dsc(:, 1));
%         mean(dsc(:, 2)) std(dsc(:, 2)) min(dsc(:, 2)) max(dsc(:, 2));
%         mean(dsc(:, 3)) std(dsc(:, 3)) min(dsc(:, 3)) max(dsc(:, 3));
%         mean(dsc(:, 4)) std(dsc(:, 4)) min(dsc(:, 4)) max(dsc(:, 4));
%         mean(dsc(:, 5)) std(dsc(:, 5)) min(dsc(:, 5)) max(dsc(:, 5));
%         mean(gdsc) std(gdsc) min(gdsc) max(gdsc)];
% end
t_out = cat(2, [subfield_list 'gdsc']', array2table(t));
t_out.Properties.VariableNames = {'subfield', 'mean', 'stdev', 'min', 'max'};
writetable(t_out, fullfile(rootDir, 'supportFiles', ['dsc_' data '.csv']));

m = mean(beh.age(find(beh.child==1)));
sd = std(beh.age(find(beh.child==1)));
mi = min(beh.age(find(beh.child==1)));
ma = max(beh.age(find(beh.child==1)));
disp(['Children, M = ' num2str(m) ', SD = ' num2str(sd) ', [' num2str(mi) ',' num2str(ma) '].'])

m = mean(beh.age(find(beh.adolescent==1)));
sd = std(beh.age(find(beh.adolescent==1)));
mi = min(beh.age(find(beh.adolescent==1)));
ma = max(beh.age(find(beh.adolescent==1)));
disp(['Adolescents, M = ' num2str(m) ', SD = ' num2str(sd) ', [' num2str(mi) ',' num2str(ma) '].'])

m = mean(beh.age(find(beh.adult==1)));
sd = std(beh.age(find(beh.adult==1)));
mi = min(beh.age(find(beh.adult==1)));
ma = max(beh.age(find(beh.adult==1)));
disp(['Adults, M = ' num2str(m) ', SD = ' num2str(sd) ', [' num2str(mi) ',' num2str(ma) '].'])
