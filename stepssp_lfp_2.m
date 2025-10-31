clc;
clear all;
close all;

addpath(genpath("C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\buzcode-master"));
addpath(genpath("C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\CellExplorer-main"));

basepath = 'E:\Set 115 - FOXO1-\WT\E0\phy_ms4'; 
basename = 'recording';

channel_order = [19 18 24 17 23 16 21 31 20 30 25 7 22 12 29 11 28 10 27 8 26 2 3 14 6 13 9 1 5 0 4 15 51 50 56 49 55 48 53 63 52 62 57 39 54 44 61 43 60 42 59 40 58 34 35 46 38 45 41 33 37 32 36 47];
%channel_order=[19 18 24 17 23 16 21 31 20 30 25 7 22 12 29 11 28 10 27 8 26 2 3 14 6 13 9 1 5 0 4 15 38 32 37 36 42 35 41 34 40 33 39] %% ORDEM PARA PFC SÓ COM 1 SHANK
cd(basepath);

tempo_inicio_2 = (4 * 15 * 60); %vai buscar o início da gravação de 2%
tempo_final_2 = (5 * 15 * 60);  %vai buscar o fim da gravação de 2%

%tempo_inicio_2_5 = (3 * 15 * 60);  
%tempo_final_2_5 = (4 * 15 * 60);

tempo_inicio_4=0;     %vai buscar o início da gravação de 4%
tempo_final_4= 15*60; %vai buscar o fim da gravação de 4%

spikes = bz_GetSpikes('basepath', basepath); %spikes do animal 
lfp = bz_GetLFP(channel_order, 'basepath', basepath, 'fromDat', true, 'downsample', 50); %lfp do animal correspondente

numcells = spikes.numcells; 
cellclass = repmat({'all'}, length(spikes.times), 1);

% Adicionar os campos necessários à estrutura spikes
spikes.sessionName = 'phy_ms4'; 
spikes.basename='recording';
if ~isfield(spikes, 'region')
    spikes.region = 'unknown';  % Nome da região
end
if ~isfield(spikes, 'numcells')
    spikes.numcells = length(spikes.times);  % Número de células baseado no número de spikes
end
save('recording.spikes.cellinfo.mat', 'spikes');
band1 = [4,12];  % Faixa de frequências de 4 a 12 Hz
band2=[12,20];
band3=[20,40];
band4=[40,90];

%mkdir('PhaseModulationFig')

% if ~isnumeric(band1) || numel(band1) ~= 2
%     error('O parâmetro "passband" deve ser um vetor numérico com 2 elementos.');
% end

% Last 15 min 
cell_metrics = loadCellMetricsBatch('basepaths',{basepath},'basenames',{basename});
%cell_metrics = CellExplorer('metrics',cell_metrics);
cell_type=cell_metrics.putativeConnections;
cell_type_2=cell_metrics.putativeCellType;

% tipos = unique(cell_metrics.putativeCellType);
tipos = {'Pyramidal Cell', 'Narrow Interneuron', 'Wide Interneuron'};
fprintf('Tipos de célula encontrados:\n');
disp(tipos);

for i = 1:length(tipos)
    tipo = tipos{i};
    idx = strcmp(cell_type_2, tipo);
    
    fprintf('Tipo: %s, células selecionadas: %d\n', tipo, sum(idx));
    
    % if sum(idx) == 0
    %     fprintf('Sem células para %s, salto o processamento\n', tipo);
    %     continue;
    % end

    spikes_tipo = spikes;
    campos = fieldnames(spikes);
    for c = 1:length(campos)
        campo = campos{c};
        val = spikes.(campo);
        if iscell(val) && length(val) == length(cell_type_2)
            spikes_tipo.(campo) = val(idx);
        elseif isnumeric(val) && size(val,1) == length(cell_type_2)
            spikes_tipo.(campo) = val(idx,:);
        end
    end

    spikes_tipo.numcells = sum(idx);
    spikes_tipo.UID = 1:spikes_tipo.numcells;
    nome_var = ['spikes_' strrep(tipo, ' ', '_')];
    assignin('base', nome_var, spikes_tipo);
    spikes_tipo = evalin('base', nome_var);

    

    PhaseLockingData_2_5_1 = bz_PhaseModulation(spikes_tipo, lfp, band1,'intervals', [tempo_inicio_2 tempo_final_2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',2);
    PhaseLockingData_2_5_2 = bz_PhaseModulation(spikes_tipo, lfp, band2,'intervals', [tempo_inicio_2 tempo_final_2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',2);
    PhaseLockingData_2_5_3 = bz_PhaseModulation(spikes_tipo, lfp, band3,'intervals', [tempo_inicio_2 tempo_final_2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',2);
    PhaseLockingData_2_5_4 = bz_PhaseModulation(spikes_tipo, lfp, band4,'intervals', [tempo_inicio_2 tempo_final_2],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',2);

    PhaseLockingData_4_1 = bz_PhaseModulation(spikes_tipo, lfp, band1,'intervals', [tempo_inicio_4 tempo_final_4],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',4);
    PhaseLockingData_4_2 = bz_PhaseModulation(spikes_tipo, lfp, band2,'intervals', [tempo_inicio_4 tempo_final_4],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',4);
    PhaseLockingData_4_3 = bz_PhaseModulation(spikes_tipo, lfp, band3,'intervals', [tempo_inicio_4 tempo_final_4],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',4);
    PhaseLockingData_4_4 = bz_PhaseModulation(spikes_tipo, lfp, band4,'intervals', [tempo_inicio_4 tempo_final_4],'samplingRate',600,'powerThresh',0,'plotting',true,'saveMat',true,'celltype2',tipo,'concentration',4);
        
end


cellIndex_list = [];
cellType_list = {};
region_list = {};
concentration_list = [];

% Determinar a região de cada célula
cell_regions = cell(spikes.numcells, 1);
for i = 1:spikes.numcells
    ch = spikes.maxWaveformCh(i);
    pos_in_order = find(channel_order == ch);
    if isempty(pos_in_order)
        cell_regions{i} = 'unknown';
    elseif pos_in_order <= 32
        cell_regions{i} = 'hip';
    else
        cell_regions{i} = 'pfc';
    end
end

% Tipos de células a incluir
%tipos = {'Pyramidal Cell', 'Narrow Interneuron', 'Wide Interneuron'};

% Criar linhas com a ordem exata de processamento
for i = 1:length(tipos)
    tipo = tipos{i};
    idx = strcmp(cell_type_2, tipo);  % células deste tipo
    indices = find(idx);  % índices reais das células desse tipo
    for c = indices(:)'  % manter ordem exata
        for conc = [2.5, 4]
            cellIndex_list(end+1,1) = c;
            cellType_list{end+1,1} = tipo;
            region_list{end+1,1} = cell_regions{c};
            concentration_list(end+1,1) = conc;
        end
    end
end


% Criar tabela
cell_table = table(cellIndex_list, cellType_list, region_list, concentration_list, ...
    'VariableNames', {'CellIndex', 'CellType', 'Region', 'Concentration'});

% Guardar tabela como ficheiro Excel
excel_filename = fullfile(basepath, [basename '_CellRegions.xlsx']);
writetable(cell_table, excel_filename);





% nome_arquivo_geral = fullfile(basepath, [basename '_cell_type.txt']);
% fileID = fopen(nome_arquivo_geral, 'w');
% 
% for i = 1:length(cellType_list)
%     fprintf(fileID, '%s\t%.1f\n', cellType_list{i}, concentration_list(i));
% end
% 
% fclose(fileID);

%%

sessionInfo = bz_getSessionInfo(basepath);

fileinfo = dir(fullfile(basepath, [basename '.dat']));
nSamples = fileinfo.bytes / 2;  % 2 bytes por amostra (int16)
approx_total_samples = round(nSamples / sessionInfo.nChannels);

fprintf('Tamanho estimado do .dat para %d canais: %.2f segundos\n', ...
    sessionInfo.nChannels, approx_total_samples / 30000);