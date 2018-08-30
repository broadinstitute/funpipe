%---------------------------------------------------------------------------------
% PLOT ALIGNMENT DENSITY ACROSS GENOME FOR MULTIPLE STRAINS
% Files required to run the function - check the Load section
% matlab 2016b
%---------------------------------------------------------------------------------

% parameter to be prespecified:
% 1. coverage density matrix: *.den
% 2. chromosome coordinates and RGB Color
% 3. centromere coordinates for specific genomes

n_samples = 64  % change this number
prefix = 'batch1_AD_progeny'
n_sc = 28
n_bins = 7560   % # sliding windows for the genome

n_samples = 15  % change this number
prefix = 'batch1_15AD_cov_fc'
n_sc = 28
n_bins = 7560   % # sliding windows for the genome



figure('rend','painters','pos',[100 10 1000 n_samples * 50])   # automate height
%---------------------------------------------------------------------------------
% load the data files
%---------------------------------------------------------------------------------
alnDen = load(strcat(prefix,'.den'));
chrColor = load('/cil/shed/sandboxes/xiaoli/fungal-pipeline/analysis/crypto/cneo_h99_jec21_chr_coord_5k.dat');  % chromosome coordinates and RGB color
% load plots/cneo_centro_coords.dat;   % centromere coordinate

%---------------------------------------------------------------------------------
% assign the data files
%---------------------------------------------------------------------------------

colNum = 7; % Initialized to the starting column number.

%---------------------------------------------------------------------------------
% plot the data
%---------------------------------------------------------------------------------

for idx = 1:n_samples
    h(idx) = subplot(n_samples, 1, idx);

    % density data for each strain plotted scaffold-wise
    for sc = 1:n_sc
        i1 = chrColor(sc, 2); % Start of scaffold
        i2 = chrColor(sc, 3); % End of scaffold

        % colours for the scaffolds assigned here.
        P = chrColor(sc, 4); % R
        Q = chrColor(sc, 5); % G
        R = chrColor(sc, 6); % B

        % plot the density data -- Uncomment stem or line
        stem(alnDen(i1:i2, 2), alnDen(i1:i2, colNum), 'MarkerSize', 0.025, 'Color', [P Q R]); % Uncomment for a stem plot.
        % plot(alnDen(i1:i2, 2), alnDen(i1:i2, colNum), 'Color', [P Q R], 'square'); % Uncomment for a line plot.
        box on; hold all;
    end

    % calibrate the normal and cnv levels
    plot(alnDen(:, 2), alnDen(:, 3), 'Color', [0.411765 0.411765 0.411765]); % Normal coverage, i.e. no gain or loss in copy number
    plot(alnDen(:, 2), alnDen(:, 4), 'Color', [0.411765 0.411765 0.411765]); % Gain of 1 copy (for diploids)
    plot(alnDen(:, 2), alnDen(:, 5), 'Color', [0.545098 0.537255 0.537255]); % Gain of 1 copy (for haploids)
    plot(alnDen(:, 2), alnDen(:, 6), 'Color', [0.545098 0.545098 0.478431]); % Gain of 2 copies (for haploids)


    set(h(idx), ...
        'XLim', [0 n_bins], ...   % need to automate
        'XTick', [], ...
        'XTickLabel', [], ...
        'YLim', [0 3], ...
        'YTickLabel', [ ], ...
        'TickDir', 'out', ...
        'TickLen', [0 0]);

    colNum = colNum+1; % col value for next colNumain
end
print(strcat(prefix, '.png'), '-dpng')

% %---------------------------------------------------------------------------------
% % plot the centromere data
% %---------------------------------------------------------------------------------
%
%for sc = 1:14
%   a = cneo_centro_coords(sc, 7);%
%   b = 0;
%   plot(a, b, 'o', 'MarkerEdgeColor', [0.129412 0.129412 0.129412],
%        'MarkerFaceColor', [0.933333 0.933333 0], 'MarkerSize', 8);
%   box on; hold all;
%end

print('pilot.pdf', '-dpdf', '-bestfit')
print('pilot.png', '-dpng')
