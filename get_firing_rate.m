basepaths = {
    "D:\Dissertacao\Set119 - IP3R2KO\WT\A0\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\WT\A1\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\WT\D0\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\WT\D1\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\WT\D2\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\WT\D3\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\WT\D4\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\KO\B0\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\KO\B1\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\KO\C0\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\KO\C1\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\KO\C2\phy_ms4",
    "D:\Dissertacao\Set119 - IP3R2KO\KO\C4\phy_ms4"
};
%   2.5 %%
%start_time = 3 * 15 * 60;
%end_time = 4 * 15 * 60;

%start_time = (4 * 15 * 60); %vai buscar o início da gravação de 2%
%end_time = (5 * 15 * 60);
%   4%%
start_time = 0;
end_time = 15 * 60;

all_results = table();
channel_order = [19 18 24 17 23 16 21 31 20 30 25 7 22 12 29 11 28 10 27 8 26 2 3 14 6 13 9 1 5 0 4 15 51 50 56 49 55 48 53 63 52 62 57 39 54 44 61 43 60 42 59 40 58 34 35 46 38 45 41 33 37 32 36 47];
%channel_order=[19 18 24 17 23 16 21 31 20 30 25 7 22 12 29 11 28 10 27 8 26 2 3 14 6 13 9 1 5 0 4 15 38 32 37 36 42 35 41 34 40 33 39] %% ORDEM PARA PFC SÓ COM 1 SHANK


for i = 1:length(basepaths)
    basepath = basepaths{i};

    % Extrair ID do animal
    tokens = regexp(basepath, 'IP3R2KO\\(WT|KO)\\([A-E]\d+)', 'tokens');
    if isempty(tokens)
        continue;
    end
    animal_id = join(tokens{1}, ' ');

    % Carregar cell_metrics e spikes
    cell_metrics_path = fullfile(basepath, 'recording.cell_metrics.cellinfo.mat');
    if exist(cell_metrics_path, 'file') ~= 2
        continue;
    end
    load(cell_metrics_path);  % carrega 'cell_metrics'

    spikes_path = fullfile(basepath, 'recording.spikes.cellinfo.mat');
    if exist(spikes_path, 'file') ~= 2
        continue;
    end
    load(spikes_path);  % carrega 'spikes'

    num_cells = length(spikes.times);
    firing_rates_interval = zeros(num_cells, 1);
    cell_types_column = cell(num_cells, 1);
    region_column = cell(num_cells, 1);

    % Verificar se temos tipos de células
    if isfield(cell_metrics, 'putativeCellType')
        putative_cell_types = cell_metrics.putativeCellType;
    else
        putative_cell_types = repmat({'desconhecido'}, num_cells, 1);
    end

    % Calcular firing rate e região
    for j = 1:num_cells
        spike_times = spikes.times{j};
        spikes_in_interval = spike_times(spike_times >= start_time & spike_times <= end_time);
        firing_rates_interval(j) = length(spikes_in_interval) / (end_time - start_time);
        cell_types_column{j} = putative_cell_types{j};

        % Determinar região com base no maxWaveformCh
        ch = spikes.maxWaveformCh(j);
        pos_in_order = find(channel_order == ch);
        if isempty(pos_in_order)
            region_column{j} = 'unknown';
        elseif pos_in_order <= 32
            region_column{j} = 'hip';
        else
            region_column{j} = 'pfc';
        end
    end

    % Criar tabela parcial
    T = table( ...
        repmat(animal_id, num_cells, 1), ...
        cell_types_column, ...
        region_column, ...
        (1:num_cells)', ...
        firing_rates_interval, ...
        'VariableNames', {'Animal', 'CellType', 'Region', 'CellID', 'FiringRate_Hz'} ...
    );

    all_results = [all_results; T];
end



writetable(all_results, "C:\Users\patri\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Prism Graphs\Set 119 - IP3R2KO\firing_rates_4.xlsx");




