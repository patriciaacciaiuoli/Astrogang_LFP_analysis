clear all
close all

 
% definicao de animais a analisar
lista_animais=["WT_D1"]

%Delta (1-4Hz); Theta (4-12Hz); Beta (12-20Hz); Low Gamma (20-40Hz).
% bandas = { [4 12; 20 90], [12 20; 30 90], [20 40; 4 30] };
bandas = { [4 12;1 4],[4 12; 4 12],[4 12;12 20],[4 12;20 40],[4 12;40 90]};
% Exemplo: [4 12; 1 4] significa fase entre 4–12 Hz e amplitude entre 1–4 Hz

% Nome associado a cada combinação (usado para guardar os resultados)
banda_name={"phase_theta_ampl_delta","phase_theta_ampl_theta","phase_theta_ampl_beta","phase_theta_ampl_lowgama","phase_theta_ampl_highgama"};


% -- ADICIONAR CAMINHO DA FUNÇÃO DE CROSS-FREQUENCY COUPLING (Buzcode) --
addpath('D:\Astrogang Ephys_Patricia\buzcode-master\buzcode-master\analysis\CrossFrequencyCoupling');

% Ciclo for para cada banda e animal
for j = 1:length(bandas)
    banda_name_i=banda_name{j};  % Nome da banda atual
    banda = bandas{j};           % Frequências atuais
    bandafreq1 = banda(1,:);     % Faixa de fase (eixo x)
    bandafreq2 = banda(2,:);     % Faixa de amplitude (eixo y)
    concentracao=2;              % Concentração

    for animal_number=(1:length(lista_animais))
        row = animal_number;
        animal=lista_animais(animal_number);  % Nome do animal atual
        path=sprintf('D:\\Outputs spikesorting\\%s\\phy_ms4',animal);
        addpath(genpath(path));   % Adiciona o caminho da gravação

        
        % Lê o ficheiro Excel com os canais pretendidos para análise
        filename=sprintf('C:\\Users\\Medicina\\MEOCloud\\Astrogang_Ephys\\ag_BrunaMatos\\Outputs Dados Automático(desvio padrão)\\%d\\canais_pretendidos.xlsx',concentracao);
        

        T = readtable(filename);                    % Lê a tabela do Excel
        raw_string = T.CanaisPretendidos{row};      % extrai a string da célula
        lista_str = strsplit(raw_string, ',');      % separa pelos vírgulas

        channel_order_1 = str2double(lista_str);    % ocnverte para números
        combinacoes = nchoosek(channel_order_1, 2); % Cria todas as combinações possíveis entre canais

        
        % Lista com a ordem dos canais
        channel_order=[19 18 24 17 23 16 21 31 20 30 25 7 22 12 29 11 28 10 27 8 26 2 3 14 6 13 9 1 5 0 4 15 51 50 56 49 55 48 53 63 52 62 57 39 54 44 61 43 60 42 59 40 58 34 35 46 38 45 41 33 37 32 36 47];
        %channel_order_1=[]
        %channel_order_2=[]
        % for i = 1:size(combinacoes, 1)
        %     idx1 = combinacoes(i,1);          % índice a buscar
        %     idx2 = combinacoes(i,2);
        %     channel_order_1(end+1) = channel_order(idx1 + 1);
        %     channel_order_2(end+1) = channel_order(idx2 + 1);
        % 
        % end

        cd(path)
        lfp = bz_GetLFP(channel_order, 'basepath', path, 'fromDat', true, 'downsample', 50);
        % all_sameProbe_hip_list = [24 22;24 27;24 14;24 0;31 22;31 27;31 14;31 0;22 14;22 0;27 14; 27 0];
        % all_sameProbe_pfc_list = [56 54;56 59;56 46;56 36;63 54;63 59;63 46;63 36;54 46;54 36;59 46;59 36];
        % all_sameShank_hip_list = [24 31;22 27;14 0];
        % all_sameShank_pfc_list = [56 63;54 59;46 36];
        % all_diff_list      = [];
        % all_sameShank_hip = zeros(50,50,0);
        % all_sameShank_pfc = zeros(50,50,0);
        % all_sameProbe_hip = zeros(50,50,0);
        % all_sameProbe_pfc = zeros(50,50,0);
        % all_diff          = zeros(50,50,0);
        % for i=(1:length(combinacoes))
        % 

        % Definição do intervalo de tempo para análise (concentração)
        % tempo_inicio_2 = 0;
        % 
        % tempo_final_2 = (15 * 60);
        tempo_inicio_2 = 4*15*60;
        tempo_final_2 = 5*15*60;
        % 
        interval=[tempo_inicio_2 tempo_final_2]
         % c1=channel_order_1(i)
         % c2=channel_order_2(i)
        c1=54;
        c2=22;
        
        % Cálculo do PAC
        comod=bz_PhaseAmpCouplingByAmp(lfp,bandafreq1,bandafreq2,'intervals', interval,'ampCh',c1,'phaseCh',c2);
        % M = comod.phaseamphist;
        % 
        % 
        % 
        % if ismember([c1 c2], all_sameProbe_hip_list, 'rows') || ismember([c2 c1], all_sameProbe_hip_list, 'rows')
        %     all_sameProbe_hip(:,:,end+1) = M;
        % elseif ismember([c1 c2], all_sameProbe_pfc_list, 'rows') || ismember([c2 c1], all_sameProbe_pfc_list, 'rows')
        %     all_sameProbe_pfc(:,:,end+1) = M;
        % elseif ismember([c1 c2], all_sameShank_hip_list, 'rows') || ismember([c2 c1], all_sameShank_hip_list, 'rows')
        %     all_sameShank_hip(:,:,end+1) = M;
        % elseif ismember([c1 c2], all_sameShank_pfc_list, 'rows') || ismember([c2 c1], all_sameShank_pfc_list, 'rows')
        %     all_sameShank_pfc(:,:,end+1) = M;
        % else 
        %     all_diff(:,:,end+1) = M;
        % end
        

        % Cria pasta de saída caso nao exista e guarda os resultados
        outputFolder = sprintf('C:\\Users\\Medicina\\MEOCloud\\Astrogang_Ephys\\ag_BrunaMatos\\Outputs Dados Automático(desvio padrão)\\%d\\pac\\%s\\%s',concentracao,banda_name_i, animal);
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end

        % Define o caminho completo para guardar a imagem
        savepath = fullfile(outputFolder, sprintf('PAC(Channel_%d_%d).png', c1,c2));
        saveas(gcf, savepath)  % ou use outro nome/formato
        % end
        % mean_sameShank_hip = mean(all_sameShank_hip, 3, 'omitnan');
        % mean_sameProbe_pfc = mean(all_sameProbe_pfc, 3, 'omitnan');
        % mean_sameShank_pfc = mean(all_sameShank_pfc, 3, 'omitnan');
        % mean_sameProbe_hip = mean(all_sameProbe_hip, 3, 'omitnan');
        % mean_diff      = mean(all_diff, 3, 'omitnan');
        % 
        % figure;
        % imagesc(comod.phasebins, comod.ampbins, mean_sameShank_pfc);
        % axis xy; colormap jet; colorbar;
        % xlabel('Phase (rad)'); ylabel('Amp (Z)');
        % title(sprintf('Média PAC - Intra-Shank - %s - %s', animal, banda_name_i));
        % saveas(gcf, fullfile(outputFolder, 'PAC_mean_sameShank_hip.png'));
        % close(gcf);
        % 
        % figure;
        % imagesc(comod.phasebins, comod.ampbins, mean_sameShank_hip);
        % axis xy; colormap jet; colorbar;
        % xlabel('Phase (rad)'); ylabel('Amp (Z)');
        % title(sprintf('Média PAC - Intra-Shank - %s - %s', animal, banda_name_i));
        % saveas(gcf, fullfile(outputFolder, 'PAC_mean_sameShank_pfc.png'));
        % close(gcf);
        % 
        % % --- Mesma Probe ---
        % figure;
        % imagesc(comod.phasebins, comod.ampbins, mean_sameProbe_hip);
        % axis xy; colormap jet; colorbar;
        % xlabel('Phase (rad)'); ylabel('Amp (Z)');
        % title(sprintf('Mean PAC - Intra-Probe - %s - %s', animal, banda_name_i));
        % saveas(gcf, fullfile(outputFolder, 'PAC_mean_sameProbe_hip.png'));
        % close(gcf);
        % 
        % figure;
        % imagesc(comod.phasebins, comod.ampbins, mean_sameProbe_pfc);
        % axis xy; colormap jet; colorbar;
        % xlabel('Phase (rad)'); ylabel('Amp (Z)');
        % title(sprintf('Mean PAC - Intra-Probe - %s - %s', animal, banda_name_i));
        % saveas(gcf, fullfile(outputFolder, 'PAC_mean_sameProbe_pfc.png'));
        % close(gcf);
        % % --- Diferentes ---
        % figure;
        % imagesc(comod.phasebins, comod.ampbins, mean_diff);
        % axis xy; colormap jet; colorbar;
        % xlabel('Phase (rad)'); ylabel('Amp (Z)');
        % title(sprintf('Mean PAC - Inter-Porbe - %s - %s', animal, banda_name_i));
        % saveas(gcf, fullfile(outputFolder, 'PAC_mean_diff.png'));
        % close(gcf);

    end
end