clearvars
close all
clc

addpath(genpath('C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\buzcode-master\buzcode-master'));
addpath(genpath('C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\CellExplorer-main\CellExplorer-main'));

lista_animais = ["WT_D0","WT_D1","WT_D2","WT_D3","WT_D4","KO_C0","KO_C1","KO_C2","KO_C4"];
concentracao = 4;
filename_temp = "C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\Canais escolhidos\Set 119 - IP3R2KO\canais_escolhidos.xlsx";
filename_mapa  = "C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\Canais removidos\Set 119 - IP3R2KO\mapa_final.xlsx";
baseDataFolder = "D:\Astrogang Ephys_Patricia\Recordings\Set 119 - IP3R2KO";
outputBase ='C:\Users\Medicina\MEOCloud\Astrogang_Ephys\ag_PatriciaAzenha\Análise dados\Outputs\Set 119 - IP3R2KO\Outputs_PAC';

% bandas = { [4 12;1 4],[4 12; 4 12],[4 12;12 20],[4 12;20 90], ...
%            [12 20;1 4],[12 20;4 12], [12 20;12 20],[12 20;20 90]};
bandas = {[12 20;20 40],[12 20;40 90]};

% banda_name = {"phase_theta_ampl_delta","phase_theta_ampl_theta","phase_theta_ampl_beta", ...
%               "phase_theta_ampl_gama","phase_beta_ampl_delta", ...
%               "phase_beta_ampl_theta","phase_beta_ampl_beta","phase_beta_ampl_gama"};
banda_name = {"phase_beta_ampl_lowgama","phase_beta_ampl_highgama"};

T_temp = readtable(filename_temp,'ReadVariableNames',true);
T_temp.Genotype = strtrim(string(T_temp.Genotype));
T_temp.Animal   = strtrim(string(T_temp.Animal));
if iscell(T_temp.Anesthesia_Level)
    T_temp.Anesthesia_Level = str2double(strtrim(T_temp.Anesthesia_Level));
end

T_mapa = readtable(filename_mapa,'ReadVariableNames',true);
T_mapa.Genotype = strtrim(string(T_mapa.Genotype));
T_mapa.Animal   = strtrim(string(T_mapa.Animal));
if iscell(T_mapa.Anesthesia_Level)
    T_mapa.Anesthesia_Level = str2double(strtrim(T_mapa.Anesthesia_Level));
end

if concentracao == 4
    start_time = 0;
    end_time = 15*60;
elseif concentracao == 2.5
    start_time = 3*15*60;
    end_time = 4*15*60;
else
    error('Concentração não reconhecida: %.2f', concentracao);
end
intervalo = [start_time end_time];

for j = 1:length(bandas)
    bandafreq1 = bandas{j}(1,:);
    bandafreq2 = bandas{j}(2,:);
    band_label = banda_name{j};

    for ai = 1:length(lista_animais)
        entry = char(lista_animais(ai));
        tok = strsplit(entry,'_');
        if numel(tok)~=2, continue; end
        genotype = tok{1};
        animal   = tok{2};

        mask_temp = (T_temp.Genotype==genotype) & (T_temp.Animal==animal) & (T_temp.Anesthesia_Level==concentracao);
        if ~any(mask_temp), continue; end

        raw_ch = string(T_temp.Selected_Channels{find(mask_temp,1)});
        if ismissing(raw_ch) || strlength(raw_ch) == 0, continue; end
        ch_list = str2double(strtrim(strsplit(raw_ch,',')));
        if any(isnan(ch_list)), continue; end

        mask_mapa = (T_mapa.Genotype==genotype) & (T_mapa.Animal==animal) & (T_mapa.Anesthesia_Level==concentracao);
        if ~any(mask_mapa), continue; end

        mapa_str = string(T_mapa{find(mask_mapa,1),4});
        mapa_vec = str2double(strtrim(strsplit(mapa_str,',')));
        if any(isnan(mapa_vec)), continue; end
        canais_efetivos = mapa_vec(ch_list+1);

        path = fullfile(baseDataFolder, genotype, animal, 'phy_ms4');
        if ~exist(path,'dir'), continue; end
        addpath(genpath(path));

        try
            cd(path);
            lfp = bz_GetLFP(canais_efetivos, 'basepath', char(path), 'fromDat', true, 'downsample', 50,'intervals', intervalo);
        catch
            continue;
        end

        outFolder = fullfile(outputBase, sprintf('conc_%g',concentracao), band_label, sprintf('%s_%s',genotype,animal));
        if ~exist(outFolder,'dir'), mkdir(outFolder); end

        pairs = [canais_efetivos(5) canais_efetivos(1); ...
                 canais_efetivos(3) canais_efetivos(9); canais_efetivos(7) canais_efetivos(11); canais_efetivos(3) canais_efetivos(3); canais_efetivos(9) canais_efetivos(9);canais_efetivos(9) canais_efetivos(3)];

        % pairs = [canais_efetivos(11) canais_efetivos(7)];

        for p = 1:size(pairs,1)
            c_amp   = pairs(p,1);
            c_phase = pairs(p,2);
        
            try
                [comod, fig] = bz_PhaseAmpCouplingByAmp(lfp, bandafreq1, bandafreq2, ...
                                    'ampCh', c_amp, 'phaseCh', c_phase, ...
                                    'intervals', intervalo, 'makePlot', true);
        
                save(fullfile(outFolder, sprintf('comod_ch%d_ch%d.mat', c_amp, c_phase)), 'comod');
        
                if ~isempty(fig)
                    savepath = fullfile(outFolder, sprintf('PAC_Channel_signal%d_modBy%d.png', c_amp, c_phase));
                    exportgraphics(fig, savepath, 'Resolution', 300);
                    close(fig);
                end
        
            catch
                continue;
            end
        end
    end
end