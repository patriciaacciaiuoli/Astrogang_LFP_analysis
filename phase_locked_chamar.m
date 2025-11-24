clear all
close all

% Adiciona ao MATLAB o caminho da toolbox Buzcode (conjunto de funções para análise de neurofisiologia)
addpath(genpath('D:\Astrogang Ephys_Patricia\buzcode-master\buzcode-master'))
% Adiciona o caminho da toolbox CellExplorer (para explorar e carregar métricas das células)
addpath(genpath('C:\Users\Medicina\Downloads\CellExplorer-main'))

% Adiciona a pasta onde estão scripts ou funções de análise personalizadas
addpath('C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_BrunaMatos\Análise de Dados');
addpath(genpath('C:\Users\Medicina\Desktop'))


% Definicao de animais para analisar
animais = {'KO_C0','WT_D1','KO_C4','WT_A3','WT_D0','WT_D1','KO_C4','WT_A3'}; 

% Diretório base onde estão guardados os resultados do spike sorting
base_dir = 'D:\Outputs spikesorting';


% ciclo for para analisar cada animal
for a = 1:length(animais)
    animal = animais{a};  % animal atual no loop
    
    % Define o caminho completo até à pasta da sessão (onde está o ficheiro phy_ms4)
    basepath = fullfile(base_dir, animal, 'phy_ms4'); 
    basename = {'recording'}; %nome base da sessao

    
    % definir ordem dos canais (mapa)
    channel_order = [19 18 24 17 23 16 21 31 20 30 25 7 22 12 29 11 28 10 27 8 26 2 3 14 6 13 9 1 5 0 4 15 51 50 56 49 55 48 53 63 52 62 57 39 54 44 61 43 60 42 59 40 58 34 35 46 38 45 41 33 37 32 36 47];

    cd(basepath)

    %Definir intervalos temporais (gravação com 15 min por percentagem :
    %protocolo 1
    c1=2 %concentracao para analisar
    c2=4 %concentracao para analisar

    % Cálculo dos tempos de início e fim de cada condição (em segundos)
    tempoinicioc1=(4-c1)/0.5*15
    tempoinicioc2=(4-c2)/0.5*15
    tempofinalc1=(4-c1)/0.5*15+15
    tempofinalc2=(4-c2)/0.5*15+15

    
    % Carrega os spikes processados pelo Buzcode
    spikes = bz_GetSpikes('basepath', basepath);

    % Carrega o LFP a partir do ficheiro .dat reamostrado a 50 Hz
    lfp = bz_GetLFP(channel_order, 'basepath', basepath, 'fromDat', true, 'downsample', 50);

    % Obtém a duração total da gravação em segundos
    duracao_segundos = lfp.timestamps(end);
    display(duracao_segundos)
    numcells = spikes.numcells; % Número de células detectadas

    % Cria um campo de classe genérica para cada célula (todas inicialmente "all")
    cellclass = repmat({'all'}, length(spikes.times), 1);

    spikes.sessionName = 'phy_ms4';
    spikes.basename = 'recording';

    % Se o campo 'region' não existir, define como 'unknown'
    if ~isfield(spikes, 'region')
        spikes.region = 'unknown';
    end

    % Se o campo 'numcells' não existir, calcula-o
    if ~isfield(spikes, 'numcells')
        spikes.numcells = length(spikes.times);
    end
    save('recording.spikes.cellinfo.mat', 'spikes'); % Guarda a estrutura spikes num ficheiro .mat

    band1 = [4,12];   % Theta
    band2 = [12,20];  % Beta
    band3 = [20,40];  % Low Gamma
    band4 = [40,90];  % High Gamma

    cell_metrics = loadCellMetricsBatch('basepaths',{basepath},'basenames',basename); % Carrega as métricas de todas as células (produzidas pelo CellExplorer)
    %cell_metrics = CellExplorer('metrics',cell_metrics);

    
    % Extrai as regiões cerebrais e tipos celulares de cada unidade
    regions_cells = cell_metrics.brainRegion;
    cell_type = cell_metrics.putativeConnections;
    cell_type_2 = cell_metrics.putativeCellType;

    tipos = {'Pyramidal Cell', 'Narrow Interneuron', 'Wide Interneuron'};
    zonas = {'CTX','HIP'}; % Regiões cérebro

    
    % ciclo for para analisar por tipo de célula e região
    for j =1:length(zonas)

        for i =1:length(tipos)
            tipo = tipos{i};
            zona = zonas{j};
            display(zona)
            idx = strcmp(cell_type_2, tipo) & strcmp(regions_cells, zona); % Índice lógico das células que pertencem a este tipo e região

            % Se não houver células deste tipo, passa para o próximo
            if ~any(idx)
                fprintf('Nenhuma célula do tipo: %s\n', tipo);
                continue;
            end

            % selecionar apenas as células do tipo atual
            spikes_tipo = spikes;
            campos = fieldnames(spikes);

            % Filtra todos os campos da estrutura spikes, mantendo apenas as células correspondentes
            for c = 1:length(campos)
                campo = campos{c};
                val = spikes.(campo);
                if iscell(val) && length(val) == length(cell_type_2)
                    spikes_tipo.(campo) = val(idx);
                elseif isnumeric(val) && size(val,1) == length(cell_type_2)
                    spikes_tipo.(campo) = val(idx,:);
                end
            end

            % Atualiza o número de células e os IDs
            spikes_tipo.numcells = sum(idx);
            spikes_tipo.UID = 1:spikes_tipo.numcells;

            % Cria uma variável com nome correspondente ao tipo celular (ex.: spikes_Pyramidal_Cell)
            nome_var = ['spikes_' strrep(tipo, ' ', '_')];
            assignin('base', nome_var, spikes_tipo);
            spikes_tipo = evalin('base', nome_var);

            % Calcular phase locking para cada banda de frequencia e para
            % as duas percentagens c1 e c2
            PhaseLockingData_c1_1 = bz_PhaseModulatio(spikes_tipo, lfp, band1,'intervals', [tempoinicioc1 tempofinalc1],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c1,'region',zona);
            PhaseLockingData_c1_2 = bz_PhaseModulatio(spikes_tipo, lfp, band2,'intervals', [tempoinicioc1 tempofinalc1],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c1,'region',zona);
            PhaseLockingData_c1_3 = bz_PhaseModulatio(spikes_tipo, lfp, band3,'intervals', [tempoinicioc1 tempofinalc1],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c1,'region',zona);
            PhaseLockingData_c1_4 = bz_PhaseModulatio(spikes_tipo, lfp, band4,'intervals', [tempoinicioc1 tempofinalc1],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c1,'region',zona);

            PhaseLockingData_c2_1 = bz_PhaseModulatio(spikes_tipo, lfp, band1,'intervals', [tempoinicioc2 tempofinalc2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c2,'region',zona);
            PhaseLockingData_c2_2 = bz_PhaseModulatio(spikes_tipo, lfp, band2,'intervals', [tempoinicioc2 tempofinalc2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c2,'region',zona);
            PhaseLockingData_c2_3 = bz_PhaseModulatio(spikes_tipo, lfp, band3,'intervals', [tempoinicioc2 tempofinalc2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c2,'region',zona);
            PhaseLockingData_c2_4 = bz_PhaseModulatio(spikes_tipo, lfp, band4,'intervals', [tempoinicioc2 tempofinalc2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',c2,'region',zona);


        end
    end

end
