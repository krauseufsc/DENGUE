% Definições das funções
inp = zeros(1, 13);
function out = field(inp, t, n)
    s = inp(1);
    i1 = inp(2);
    i2 = inp(3);
    r1 = inp(4);
    r2 = inp(5);
    s1 = inp(6);
    s2 = inp(7);
    i12 = inp(8);
    i21 = inp(9);
    r = inp(10);
    sv = inp(11);
    v1 = inp(12);
    v2 = inp(13);

    m = sv + v1 + v2;
    phi = 0.8;
    mu = 1 / 65;
    alfa = 2;
    gama = 52;
    nu = 36.5;
    teta = 2 * nu;
    omega = 2 * pi * 6;
    xi = nu * (1 + 0.4 * cos(omega * t));
    beta = 2 * gama;

    ds = -(beta / m) * s * (v1 + v2) + mu * (n - s);
    di1 = (beta / m) * s * v1 - (gama + mu) * i1;
    di2 = (beta / m) * s * v2 - (gama + mu) * i2;
    dr1 = gama * i1 - (alfa + mu) * r1;
    dr2 = gama * i2 - (alfa + mu) * r2;
    ds1 = -(beta / m) * s1 * v2 + alfa * r1 - mu * s1;
    ds2 = -(beta / m) * s2 * v1 + alfa * r2 - mu * s2;
    di12 = (beta / m) * s1 * v2 - (gama + mu) * i12;
    di21 = (beta / m) * s2 * v1 - (gama + mu) * i21;
    dr = gama * (i12 + i21) - mu * r;
    dsv = -(teta / n) * sv * (i1 + i2 + phi * (i12 + i21)) + xi * m - nu * sv;
    dv1 = (teta / n) * sv * (i1 + phi * i21) - nu * v1;
    dv2 = (teta / n) * sv * (i2 + phi * i12) - nu * v2;

    out = [ds, di1, di2, dr1, dr2, ds1, ds2, di12, di21, dr, dsv, dv1, dv2];
endfunction

function armazenar_resultados(sim, k, in, S_todas, I1_todas, I2_todas, R1_todas, R2_todas, ...
                             S1_todas, S2_todas, I12_todas, I21_todas, R_todas, SV_todas, V1_todas, V2_todas)
    % Armazena os resultados de todas as variáveis para a simulação atual
    idx = k - 18250;
    S_todas(sim, idx) = in(1);
    I1_todas(sim, idx) = in(2);
    I2_todas(sim, idx) = in(3);
    R1_todas(sim, idx) = in(4);
    R2_todas(sim, idx) = in(5);
    S1_todas(sim, idx) = in(6);
    S2_todas(sim, idx) = in(7);
    I12_todas(sim, idx) = in(8);
    I21_todas(sim, idx) = in(9);
    R_todas(sim, idx) = in(10);
    SV_todas(sim, idx) = in(11);
    V1_todas(sim, idx) = in(12);
    V2_todas(sim, idx) = in(13);
endfunction

function out = rk(inp, t, dt, n)
    k1 = field(inp, t, n);
    k2 = field(inp + (dt / 2) * k1, t + dt / 2, n);
    k3 = field(inp + (dt / 2) * k2, t + dt / 2, n);
    k4 = field(inp + dt * k3, t + dt, n);

    temp = k1 + 2 * k2 + 2 * k3 + k4;
    out = max(0, inp + (dt / 6) * temp);
endfunction

function in = aplicar_erro_aleatorio(in, sim)
    % Aplica erro aleatório de ±5% em todas as variáveis
    fatores_erro = 0.95 + 0.10 * rand(1, 13);
    in = max(0, in .* fatores_erro);
    fprintf('  Simulação %d: Erro aleatório de ±5%% aplicado no passo 18250\n', sim);
endfunction

function [S_todas, I1_todas, I2_todas, R1_todas, R2_todas, S1_todas, S2_todas, ...
          I12_todas, I21_todas, R_todas, SV_todas, V1_todas, V2_todas] = ...
          executar_simulacao(sim, condicoes_iniciais, dt, kk, S_todas, I1_todas, I2_todas, ...
                            R1_todas, R2_todas, S1_todas, S2_todas, I12_todas, I21_todas, ...
                            R_todas, SV_todas, V1_todas, V2_todas)
    % Executa uma simulação completa
    fprintf('Executando simulação %d...\n', sim);

    in = condicoes_iniciais(1, :);
    tempo = 0;
    n0 = sum(in(1:10));

    for k = 1:kk
        in = rk(in, tempo, dt, n0);
        tempo = tempo + dt;

        if k == 18250
            in = aplicar_erro_aleatorio(in, sim);
        end

        if k > 18250
            idx = k - 18250;
            S_todas(sim, idx) = in(1);
            I1_todas(sim, idx) = in(2);
            I2_todas(sim, idx) = in(3);
            R1_todas(sim, idx) = in(4);
            R2_todas(sim, idx) = in(5);
            S1_todas(sim, idx) = in(6);
            S2_todas(sim, idx) = in(7);
            I12_todas(sim, idx) = in(8);
            I21_todas(sim, idx) = in(9);
            R_todas(sim, idx) = in(10);
            SV_todas(sim, idx) = in(11);
            V1_todas(sim, idx) = in(12);
            V2_todas(sim, idx) = in(13);
        end
    end
endfunction

% Main script - 5 diferentes cenários
qq = 5;

% Todas as simulações começam com as MESMAS condições iniciais
condicoes_iniciais = zeros(qq, 13);

% Condições iniciais idênticas para todas as simulações
for sim = 1:qq
    condicoes_iniciais(sim, :) = [700, 200, 100, 0, 0, 0, 0, 0, 0, 0, 9000, 500, 500];
end

dt = 1 / 365;
kk = 100 * 365;

% Arrays para armazenar resultados de todas as simulações
% Dimensão: (número de simulações, número de passos)
S_todas = zeros(qq, kk/2);
I1_todas = zeros(qq, kk/2);
I2_todas = zeros(qq, kk/2);
R1_todas = zeros(qq, kk/2);
R2_todas = zeros(qq, kk/2);
S1_todas = zeros(qq, kk/2);
S2_todas = zeros(qq, kk/2);
I12_todas = zeros(qq, kk/2);
I21_todas = zeros(qq, kk/2);
R_todas = zeros(qq, kk/2);
SV_todas = zeros(qq, kk/2);
V1_todas = zeros(qq, kk/2);
V2_todas = zeros(qq, kk/2);

% Executar cada simulação
for sim = 1:qq
    [S_todas, I1_todas, I2_todas, R1_todas, R2_todas, S1_todas, S2_todas, ...
     I12_todas, I21_todas, R_todas, SV_todas, V1_todas, V2_todas] = ...
     executar_simulacao(sim, condicoes_iniciais, dt, kk, S_todas, I1_todas, I2_todas, ...
                       R1_todas, R2_todas, S1_todas, S2_todas, I12_todas, I21_todas, ...
                       R_todas, SV_todas, V1_todas, V2_todas);
end

fprintf('Todas as simulações concluídas!\n');

% Salvar resultados em arquivo
save('resultados_multiplas_simulacoes.mat', 'S_todas', 'I1_todas', 'I2_todas', ...
     'R1_todas', 'R2_todas', 'S1_todas', 'S2_todas', 'I12_todas', 'I21_todas', ...
     'R_todas', 'SV_todas', 'V1_todas', 'V2_todas');

fprintf('Resultados salvos em resultados_multiplas_simulacoes.mat\n');

% Plots comparativos
passos = 1:kk/2;
cores = {'b', 'r', 'g', 'm', 'c'};
legendas = {'Sim 1', 'Sim 2', 'Sim 3', 'Sim 4', 'Sim 5'};

figure('Position', [100, 100, 1400, 1000]);

subplot(4,4,1)
hold on;
for sim = 1:qq
    plot(passos, S_todas(sim, :), cores{sim});
end
title('S (Naive Humans)');
xlabel('passos');
ylabel('S');
legend(legendas, 'Location', 'best', 'FontSize', 6);
hold off;

subplot(4,4,2)
hold on;
for sim = 1:qq
    plot(passos, I1_todas(sim, :), cores{sim});
end
title('I1 (Primary Inf Strain 1)');
xlabel('passos');
ylabel('I1');
hold off;

subplot(4,4,3)
hold on;
for sim = 1:qq
    plot(passos, I2_todas(sim, :), cores{sim});
end
title('I2 (Primary Inf Strain 2)');
xlabel('passos');
ylabel('I2');
hold off;

subplot(4,4,4)
hold on;
for sim = 1:qq
    plot(passos, R1_todas(sim, :), cores{sim});
end
title('R1 (Cross Immune from 1)');
xlabel('passos');
ylabel('R1');
hold off;

subplot(4,4,5)
hold on;
for sim = 1:qq
    plot(passos, R2_todas(sim, :), cores{sim});
end
title('R2 (Cross Immune from 2)');
xlabel('passos');
ylabel('R2');
hold off;

subplot(4,4,6)
hold on;
for sim = 1:qq
    plot(passos, S1_todas(sim, :), cores{sim});
end
title('S1 (Susceptible to 2)');
xlabel('passos');
ylabel('S1');
hold off;

subplot(4,4,7)
hold on;
for sim = 1:qq
    plot(passos, S2_todas(sim, :), cores{sim});
end
title('S2 (Susceptible to 1)');
xlabel('passos');
ylabel('S2');
hold off;

subplot(4,4,8)
hold on;
for sim = 1:qq
    plot(passos, I12_todas(sim, :), cores{sim});
end
title('I12 (Secondary Inf Strain 2)');
xlabel('passos');
ylabel('I12');
hold off;

subplot(4,4,9)
hold on;
for sim = 1:qq
    plot(passos, I21_todas(sim, :), cores{sim});
end
title('I21 (Secondary Inf Strain 1)');
xlabel('passos');
ylabel('I21');
hold off;

subplot(4,4,10)
hold on;
for sim = 1:qq
    plot(passos, R_todas(sim, :), cores{sim});
end
title('R (Totally recovered)');
xlabel('passos');
ylabel('R');
hold off;

subplot(4,4,11)
hold on;
for sim = 1:qq
    plot(passos, SV_todas(sim, :), cores{sim});
end
title('SV (Susceptible vectors)');
xlabel('passos');
ylabel('SV');
hold off;

subplot(4,4,12)
hold on;
for sim = 1:qq
    plot(passos, V1_todas(sim, :), cores{sim});
end
title('V1 (Vectors Strain 1)');
xlabel('passos');
ylabel('V1');
hold off;

subplot(4,4,13)
hold on;
for sim = 1:qq
    plot(passos, V2_todas(sim, :), cores{sim});
end
title('V2 (Vectors Strain 2)');
xlabel('passos');
ylabel('V2');
hold off;

print -dpng 'comparacao_5_simulacoes.png';
fprintf('Gráfico salvo em comparacao_5_simulacoes.png\n');

% Exibir estatísticas finais
fprintf('\n=== ESTATÍSTICAS FINAIS (último passo) ===\n');
variaveis = {'S', 'I1', 'I2', 'R1', 'R2', 'S1', 'S2', 'I12', 'I21', 'R', 'SV', 'V1', 'V2'};
dados_finais = {S_todas, I1_todas, I2_todas, R1_todas, R2_todas, S1_todas, S2_todas, ...
                I12_todas, I21_todas, R_todas, SV_todas, V1_todas, V2_todas};

for i = 1:length(variaveis)
    fprintf('\n%s:\n', variaveis{i});
    for sim = 1:qq
        fprintf('  Sim %d: %.2f\n', sim, dados_finais{i}(sim, end));
    end
end
