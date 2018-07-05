%---------------------------------------------------------------------------------
% PLOT ALIGNMENT DENSITY ACROSS GENOME FOR MULTIPLE STRAINS
% Files required to run the function - check the Load section
% matlab 2016b
%---------------------------------------------------------------------------------
figure('rend','painters','pos',[100 10 1200 4000])
%---------------------------------------------------------------------------------
% load the data files
%---------------------------------------------------------------------------------
load batch1_75_AD_A.den;
% load pilot_cov_fc.den; % Normalized Alignment Cover Data
load ../JEC21_NCBI/ploidy/plots/cneo_sc_color.dat;
load plots/cneo_centro_coords.dat;

%---------------------------------------------------------------------------------
% assign the data files
%---------------------------------------------------------------------------------

alnDen = batch1_75_AD_A; % Assign the plot subset.den to the 'alnDen' variable
colNum = 7; % Initialized to the starting strain column number.

%---------------------------------------------------------------------------------
% plot the data
%---------------------------------------------------------------------------------

n_samples = 75
for idx = 1:n_samples    % change this number
    disp(idx)
    h(idx) = subplot(n_samples, 1, idx);

    % density data for each strain plotted scaffold-wise
    for sc = 1:14
        i1 = cneo_sc_color(sc, 2); % Start of scaffold
        i2 = cneo_sc_color(sc, 3); % End of scaffold

        % colours for the scaffolds assigned here.
        P = cneo_sc_color(sc, 4); % R
        Q = cneo_sc_color(sc, 5); % G
        R = cneo_sc_color(sc, 6); % B

        % plot the density data -- Uncomment stem or line
        % stem(alnDen(i1:i2, 2), alnDen(i1:i2, colNum), 'MarkerSize', 0.025, 'Color', [P Q R]); % Uncomment for a stem plot.
        plot(alnDen(i1:i2, 2), alnDen(i1:i2, colNum), 'Color', [P Q R], 'square'); % Uncomment for a line plot.
        box on; hold all;
    end

    % calibrate the normal and cnv levels
    plot (alnDen(:, 2), alnDen(:, 3), 'Color', [0.411765 0.411765 0.411765]); % Normal coverage, i.e. no gain or loss in copy number
    plot (alnDen(:, 2), alnDen(:, 4), 'Color', [0.411765 0.411765 0.411765]); % Gain of 1 copy (for diploids)
    plot (alnDen(:, 2), alnDen(:, 5), 'Color', [0.545098 0.537255 0.537255]); % Gain of 1 copy (for haploids)
    plot (alnDen(:, 2), alnDen(:, 6), 'Color', [0.545098 0.545098 0.478431]); % Gain of 2 copies (for haploids)


    set(h(idx), ...
        'XLim', [0 3781], ...
        'XTick', [], ...
        'XTickLabel', [], ...
        'YLim', [0 3], ...
        'YTickLabel', [ ], ...
        'TickDir', 'out', ...
        'TickLen', [0 0]);

    colNum = colNum+1; % col value for next colNumain
end
print('batch1_75_AD_A.png', '-dpng')

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
