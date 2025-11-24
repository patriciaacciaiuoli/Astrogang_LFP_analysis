clear all
close all

lista_animais=["KO_C0" "KO_C1" "KO_C2" "KO_C5" "WT_A0" "WT_D0" "WT_D1" "WT_D2" "WT_D9"]
animal_names=["KO C0" "KO C1" "KO C2" "KO C5" "WT A0" "WT D0" "WT D1" "WT D2" "WT D9"] % nomes para titulos dos gráficos

all_sumphase = []; % Matriz para guardar PAC médio de cada animal
phasebins_total_2_global = []; 
c=2;

% canais a analisar
c1=22;
c2=54;


for animal=1:length(lista_animais)
    animal_name=lista_animais(animal)
    nome=animal_names(animal)

    % Caminho da pasta onde estão os PACs individuais salvos
    input_folder=sprintf("C:\\Users\\Medicina\\MEOCloud\\Astrogang_Ephys\\ag_BrunaMatos\\Outputs Dados Automático(desvio padrão)\\%d\\Outputs_Pac\\%s",c,animal_name)
    addpath(input_folder)  % Adiciona pasta ao MATLAB path

    % Nome do arquivo .mat com o PAC do animal
    data_set_folder=sprintf('%s_phase_theta_ampl_gama_channel1%d_channel2%d.mat',animal_name,c1,c2)
    dataset = load(data_set_folder); % Carrega o arquivo

    % extração de bins de amplitude e de fase
    phasebins_total=dataset.comod.ampbins;
    phasebins_total_2 = dataset.comod.phasebins(1:2:end);
     if isempty(phasebins_total_2_global)
        phasebins_total_2_global = phasebins_total_2; % salva apenas uma vez
    end
    sumphase_list=[] % Vetor para guardar PAC médio do animal por bin de fase

    % Cilo for para calcular pac medio por bin de fase
    for phasebin=1:2:(length(phasebins_total))

        % Seleciona os valores do PAC para os dois bins de amplitude correspondentes
        lista=(dataset.comod.phaseamphist(:, phasebin:phasebin+1));
        lista=lista(:);
        soma=mean(lista,'omitnan'); % calcular a média
        sumphase_list(end+1)=soma;
    end

    % Preenche matriz all_sumphase garantindo mesma dimensão para todos os animais
    maxlen = length(phasebins_total_2_global);
    maxlen = length(phasebins_total_2_global);
    temp = nan(1, maxlen);
    temp(1:length(sumphase_list)) = sumphase_list;
    all_sumphase = [all_sumphase; temp];
    


    % PLOT do pac individual(polar))
    figure;
    pax = polaraxes;
    
    hold on;
    
    
    % Loop barras radiais do pac
    for k = 1:length(sumphase_list)
    
    theta = [phasebins_total_2(k)-0.125, phasebins_total_2(k)+0.125]; % width of bar
    r = [0 sumphase_list(k)];


    if contains(animal_name, 'KO')
        color='r' % cor vermelha
    else
        color='k'  %cor preta
    end
    polarplot([theta(1) theta(1)], [0 r(2)], color, 'LineWidth', 2); % ado esquerdo
    polarplot([theta(2) theta(2)], [0 r(2)], color, 'LineWidth', 2); % lado direito
    polarplot([theta(1) theta(2)], [r(2) r(2)],color, 'LineWidth', 2); % topo
    end
    polarplot(phasebins_total_2, sumphase_list, '-o','Color', color); % linha que une ponntos médios

    title(nome);
    hold off;
    fig=sprintf("output%d%d_outliers.png",c1,c2)
    saveas(gcf, fullfile(input_folder, fig)); % guarda figuras individuais
    
end

% Cálculo do pac medio entre todos os animais
mean_sumphase = mean(all_sumphase, 1, 'omitnan');

% Calcula média separada de KO e WT:
is_KO = contains(lista_animais, "KO");
mean_KO = mean(all_sumphase(is_KO,:), 1, 'omitnan');
mean_WT = mean(all_sumphase(~is_KO,:), 1, 'omitnan');

figure;
pax = polaraxes;

hold on;

% --- KO (vermelho) ---
for k = 1:length(mean_KO)
    theta = [phasebins_total_2_global(k)-0.125, phasebins_total_2_global(k)+0.125];
    r = [0 mean_KO(k)];
    hKO = polarplot([theta(1) theta(2)], [r(2) r(2)], 'r', 'LineWidth', 2); % handle de exemplo
    polarplot([theta(1) theta(1)], [0 r(2)], 'r', 'LineWidth', 2);
    polarplot([theta(2) theta(2)], [0 r(2)], 'r', 'LineWidth', 2);
end

% --- WT (preto) ---
for k = 1:length(mean_WT)
    theta = [phasebins_total_2_global(k)-0.125, phasebins_total_2_global(k)+0.125];
    r = [0 mean_WT(k)];
    polarplot([theta(1) theta(1)], [0 r(2)], 'k', 'LineWidth', 2);
    polarplot([theta(2) theta(2)], [0 r(2)], 'k', 'LineWidth', 2);
    hWT = polarplot([theta(1) theta(2)], [r(2) r(2)], 'k', 'LineWidth', 2); % handle de exemplo
end

hWTline = polarplot(phasebins_total_2_global, mean_WT, '-o', 'Color', 'k', 'LineWidth', 2);
hKOline = polarplot(phasebins_total_2_global, mean_KO, '-o', 'Color', 'r', 'LineWidth', 2);
title('Mean KO vs WT');

% --- Legenda com cores corretas ---
legend([hKOline, hWTline], {'\color{red}KO', '\color{black}WT'}, 'TextColor', 'black', 'Location', 'best');

hold off;

% gardar figura média e em excel
figmed=sprintf("outputmed%d%d.png",c1,c2)
input_folder2=sprintf("C:\\Users\\Medicina\\MEOCloud\\Astrogang_Ephys\\ag_BrunaMatos\\Outputs Dados Automático(desvio padrão)\\%d\\Outputs_Pac_outliers",c)
saveas(gcf, fullfile(input_folder2, figmed));

output_csv = fullfile(input_folder2, sprintf("dados_PAC_%d%d.csv", c1, c2));

T = table(phasebins_total_2_global(:), mean_WT(:), mean_KO(:), ...
    'VariableNames', {'Phase', 'PAC_WT', 'PAC_KO'});

writetable(T, output_csv);


% Caminho para sguardar
input_folder2 = sprintf("C:\\Users\\Medicina\\MEOCloud\\Astrogang_Ephys\\ag_BrunaMatos\\Outputs Dados Automático(desvio padrão)\\%d\\Outputs_Pac", c);
output_csv = fullfile(input_folder2, sprintf("dados_PAC_KO_WT_%d%d_outliers.xlsx", c1, c2));

% Identifica quais animais são KO e WT
is_KO = contains(lista_animais, "KO");
is_WT = ~is_KO;

% Separa os valores
valores_KO = all_sumphase(is_KO, :);   % linhas = animais KO
valores_WT = all_sumphase(is_WT, :);   % linhas = animais WT

% Transpõe para colunas por fase
T = table(phasebins_total_2_global(:), 'VariableNames', {'Phase'});

% Adiciona colunas KO
for i = 1:size(valores_KO, 1)
    varname = sprintf('KO_%d', i);
    T.(varname) = valores_KO(i, :)';
end

% Adiciona colunas WT
for i = 1:size(valores_WT, 1)
    varname = sprintf('WT_%d', i);
    T.(varname) = valores_WT(i, :)';
end

% Guarda no Excel
writetable(T, output_csv);

disp("Arquivo Excel salvo (com todas as fases e valores individuais de KO e WT):");
disp(output_csv);
