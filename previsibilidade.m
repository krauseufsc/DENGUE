% ===================================================================
% Modelo epidemiológico com dois sorotipos - versão PARALELIZADA
%
% Estratégia de paralelização:
%   - Sim 1 roda SEQUENCIALMENTE (as demais dependem de seus I12/I21)
%   - Sims 2..qq rodam em PARALELO via parfor (parallel package)
%     recebendo I12_ref e I21_ref da sim 1 como entrada somente-leitura
%
% Requer: pkg install -forge parallel  (se ainda não instalado)
% ===================================================================

% ---------------------------------------------------------------------------
% Função do campo vetorial do modelo
% ---------------------------------------------------------------------------
inp = zeros(1, 13);
function out = field(inp, t, n)
  s   = inp(1);  i1  = inp(2);  i2  = inp(3);
  r1  = inp(4);  r2  = inp(5);  s1  = inp(6);
  s2  = inp(7);  i12 = inp(8);  i21 = inp(9);
  r   = inp(10); sv  = inp(11); v1  = inp(12); v2 = inp(13);

  m     = sv + v1 + v2;
  phi   = 0.8;
  mu    = 1 / 65;
  alfa  = 2;
  gama  = 52;
  nu    = 36.5;
  teta  = 2 * nu;
  omega = 2 * pi * 6;
  xi    = nu * (1 + 0.4 * cos(omega * t));
  beta  = 2 * gama;

  ds   = -(beta/m)*s*(v1+v2)              + mu*(n-s);
  di1  =  (beta/m)*s*v1                   - (gama+mu)*i1;
  di2  =  (beta/m)*s*v2                   - (gama+mu)*i2;
  dr1  =  gama*i1                         - (alfa+mu)*r1;
  dr2  =  gama*i2                         - (alfa+mu)*r2;
  ds1  = -(beta/m)*s1*v2  + alfa*r1       - mu*s1;
  ds2  = -(beta/m)*s2*v1  + alfa*r2       - mu*s2;
  di12 =  (beta/m)*s1*v2                  - (gama+mu)*i12;
  di21 =  (beta/m)*s2*v1                  - (gama+mu)*i21;
  dr   =  gama*(i12+i21)                  - mu*r;
  dsv  = -(teta/n)*sv*(i1+i2+phi*(i12+i21)) + xi*m - nu*sv;
  dv1  =  (teta/n)*sv*(i1+phi*i21)        - nu*v1;
  dv2  =  (teta/n)*sv*(i2+phi*i12)        - nu*v2;

  out = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2];
endfunction

% ---------------------------------------------------------------------------
% Integrador Runge-Kutta de 4ª ordem
% ---------------------------------------------------------------------------
function out = rk(inp, t, dt, n)
  k1 = field(inp,                t,        n);
  k2 = field(inp + (dt/2)*k1,   t + dt/2, n);
  k3 = field(inp + (dt/2)*k2,   t + dt/2, n);
  k4 = field(inp + dt*k3,        t + dt,   n);
  out = max(0, inp + (dt/6)*(k1 + 2*k2 + 2*k3 + k4));
endfunction

% ---------------------------------------------------------------------------
% Aplica perturbação aleatória de ±5% no instante de bifurcação
% ---------------------------------------------------------------------------
function in = aplicar_erro_aleatorio(in, sim)
  fatores_erro = 0.95 + 0.10 * rand(1, 13);
  in = max(0, in .* fatores_erro);
  fprintf('  Simulação %d: perturbação ±5%% aplicada no passo 18250\n', sim);
endfunction

% ---------------------------------------------------------------------------
% Executa uma simulação completa e devolve struct com séries temporais
%
% Parâmetros:
%   sim       - índice da simulação (usado apenas para log/perturbação)
%   cond_ini  - vetor 1×13 de condições iniciais
%   dt, kk    - passo de tempo e total de passos
%   I12_ref   - série I12 da sim 1 ([] quando a própria sim 1 estiver rodando)
%   I21_ref   - série I21 da sim 1 ([] quando a própria sim 1 estiver rodando)
%
% As sims dependentes recebem I12_ref/I21_ref prontos (somente leitura),
% eliminando qualquer condição de corrida e permitindo parfor.
% ---------------------------------------------------------------------------
function resultado = executar_simulacao(sim, cond_ini, dt, kk, I12_ref, I21_ref)
  usar_ref = ~isempty(I12_ref);   % false somente para a simulação 1
  num_post = kk - 18250;          % passos gravados após a perturbação

  in    = cond_ini;
  tempo = 0;
  n0    = sum(in(1:10));

  % Pré-aloca arrays locais (evita realocação dentro do loop)
  S_loc   = zeros(1, num_post);
  I1_loc  = zeros(1, num_post);
  I2_loc  = zeros(1, num_post);
  R1_loc  = zeros(1, num_post);
  R2_loc  = zeros(1, num_post);
  S1_loc  = zeros(1, num_post);
  S2_loc  = zeros(1, num_post);
  I12_loc = zeros(1, num_post);
  I21_loc = zeros(1, num_post);
  R_loc   = zeros(1, num_post);
  SV_loc  = zeros(1, num_post);
  V1_loc  = zeros(1, num_post);
  V2_loc  = zeros(1, num_post);

  for k = 1:kk
    in    = rk(in, tempo, dt, n0);
    tempo = tempo + dt;

    % Aplica perturbação no passo exato da bifurcação
    if k == 18250
      in = aplicar_erro_aleatorio(in, sim);
    end

    % Grava e, para sims dependentes, força I12/I21 da referência
    if k > 18250
      idx = k - 18250;
      if usar_ref
        in(8) = I12_ref(idx);   % I12 travado na trajetória da sim 1
        in(9) = I21_ref(idx);   % I21 travado na trajetória da sim 1
      end
      S_loc(idx)   = in(1);
      I1_loc(idx)  = in(2);
      I2_loc(idx)  = in(3);
      R1_loc(idx)  = in(4);
      R2_loc(idx)  = in(5);
      S1_loc(idx)  = in(6);
      S2_loc(idx)  = in(7);
      I12_loc(idx) = in(8);
      I21_loc(idx) = in(9);
      R_loc(idx)   = in(10);
      SV_loc(idx)  = in(11);
      V1_loc(idx)  = in(12);
      V2_loc(idx)  = in(13);
    end
  end

  % Empacota em struct para retorno limpo (sem passar arrays gigantes por ref)
  resultado = struct( ...
    'S',   S_loc,   'I1',  I1_loc,  'I2',  I2_loc,  ...
    'R1',  R1_loc,  'R2',  R2_loc,  'S1',  S1_loc,  ...
    'S2',  S2_loc,  'I12', I12_loc, 'I21', I21_loc, ...
    'R',   R_loc,   'SV',  SV_loc,  'V1',  V1_loc,  'V2',  V2_loc);
endfunction

% ===========================================================================
%  SCRIPT PRINCIPAL
% ===========================================================================

pkg load parallel;   % carrega suporte a parfor / workers

qq  = 5;             % número de cenários
dt  = 1 / 365;
kk  = 100 * 365;     % 36 500 passos totais
num_post = kk - 18250;   % 18 250 passos gravados

cond_ini_base = [700, 200, 100, 0, 0, 0, 0, 0, 0, 0, 9000, 500, 500];

% ---------------------------------------------------------------------------
% FASE 1 — Simulação 1 (sequencial, obrigatório: as demais dependem dela)
% ---------------------------------------------------------------------------
fprintf('=== FASE 1: Simulação 1 (referência, sequencial) ===\n');
tic;
res1 = executar_simulacao(1, cond_ini_base, dt, kk, [], []);
fprintf('Simulação 1 concluída em %.1f s\n\n', toc);

% Extrai as séries de referência que as demais simulações irão ler
I12_ref = res1.I12;
I21_ref = res1.I21;

% ---------------------------------------------------------------------------
% FASE 2 — Simulações 2..qq em paralelo
%   Cada worker recebe I12_ref e I21_ref como cópias somente-leitura;
%   não há escrita compartilhada → sem condição de corrida.
% ---------------------------------------------------------------------------
fprintf('=== FASE 2: Simulações 2..%d em paralelo ===\n', qq);
resultados = cell(qq, 1);
resultados{1} = res1;   % já temos o resultado da sim 1

tic;
parfor p = 2:qq
  fprintf('Worker iniciando simulação %d...\n', p);
  resultados{p} = executar_simulacao(p, cond_ini_base, dt, kk, I12_ref, I21_ref);
end
fprintf('Fase paralela concluída em %.1f s\n\n', toc);

fprintf('Todas as simulações concluídas!\n');

% ---------------------------------------------------------------------------
% Monta matrizes de resultado a partir dos structs coletados
% ---------------------------------------------------------------------------
S_todas   = zeros(qq, num_post);
I1_todas  = zeros(qq, num_post);
I2_todas  = zeros(qq, num_post);
R1_todas  = zeros(qq, num_post);
R2_todas  = zeros(qq, num_post);
S1_todas  = zeros(qq, num_post);
S2_todas  = zeros(qq, num_post);
I12_todas = zeros(qq, num_post);
I21_todas = zeros(qq, num_post);
R_todas   = zeros(qq, num_post);
SV_todas  = zeros(qq, num_post);
V1_todas  = zeros(qq, num_post);
V2_todas  = zeros(qq, num_post);

for sim = 1:qq
  r = resultados{sim};
  S_todas(sim,:)   = r.S;
  I1_todas(sim,:)  = r.I1;
  I2_todas(sim,:)  = r.I2;
  R1_todas(sim,:)  = r.R1;
  R2_todas(sim,:)  = r.R2;
  S1_todas(sim,:)  = r.S1;
  S2_todas(sim,:)  = r.S2;
  I12_todas(sim,:) = r.I12;
  I21_todas(sim,:) = r.I21;
  R_todas(sim,:)   = r.R;
  SV_todas(sim,:)  = r.SV;
  V1_todas(sim,:)  = r.V1;
  V2_todas(sim,:)  = r.V2;
end

% ---------------------------------------------------------------------------
% Salva resultados
% ---------------------------------------------------------------------------
save('resultados_multiplas_simulacoes.mat', ...
     'S_todas', 'I1_todas', 'I2_todas', ...
     'R1_todas', 'R2_todas', 'S1_todas', 'S2_todas', ...
     'I12_todas', 'I21_todas', 'R_todas', ...
     'SV_todas', 'V1_todas', 'V2_todas');
fprintf('Resultados salvos em resultados_multiplas_simulacoes.mat\n');

% ---------------------------------------------------------------------------
% Plots comparativos
% ---------------------------------------------------------------------------
passos  = 1:num_post;
cores   = {'b', 'r', 'g', 'm', 'c'};
legendas = {'Sim 1', 'Sim 2', 'Sim 3', 'Sim 4', 'Sim 5'};

figure('Position', [100, 100, 1400, 1000]);

subplot(5,5,1);  hold on;
for sim = 1:qq; plot(passos, S_todas(sim,:),   cores{sim}); end
title('S (Naive Humans)');        xlabel('passos'); ylabel('S');
legend(legendas(1:qq), 'Location', 'best', 'FontSize', 6); hold off;

subplot(5,5,2);  hold on;
for sim = 1:qq; plot(passos, I1_todas(sim,:),  cores{sim}); end
title('I1 (Primary Inf Strain 1)');  xlabel('passos'); ylabel('I1');  hold off;

subplot(5,5,3);  hold on;
for sim = 1:qq; plot(passos, I2_todas(sim,:),  cores{sim}); end
title('I2 (Primary Inf Strain 2)');  xlabel('passos'); ylabel('I2');  hold off;

subplot(5,5,4);  hold on;
for sim = 1:qq; plot(passos, R1_todas(sim,:),  cores{sim}); end
title('R1 (Cross Immune from 1)');   xlabel('passos'); ylabel('R1'); hold off;

subplot(5,5,5);  hold on;
for sim = 1:qq; plot(passos, R2_todas(sim,:),  cores{sim}); end
title('R2 (Cross Immune from 2)');   xlabel('passos'); ylabel('R2'); hold off;

subplot(5,5,6);  hold on;
for sim = 1:qq; plot(passos, S1_todas(sim,:),  cores{sim}); end
title('S1 (Susceptible to 2)');      xlabel('passos'); ylabel('S1'); hold off;

subplot(5,5,7);  hold on;
for sim = 1:qq; plot(passos, S2_todas(sim,:),  cores{sim}); end
title('S2 (Susceptible to 1)');      xlabel('passos'); ylabel('S2'); hold off;

subplot(5,5,8);  hold on;
for sim = 1:qq; plot(passos, I12_todas(sim,:), cores{sim}); end
title('I12 (Secondary Inf Strain 2)'); xlabel('passos'); ylabel('I12'); hold off;

subplot(5,5,9);  hold on;
for sim = 1:qq; plot(passos, I21_todas(sim,:), cores{sim}); end
title('I21 (Secondary Inf Strain 1)'); xlabel('passos'); ylabel('I21'); hold off;

subplot(5,5,10); hold on;
for sim = 1:qq; plot(passos, R_todas(sim,:),   cores{sim}); end
title('R (Totally recovered)');      xlabel('passos'); ylabel('R');  hold off;

subplot(5,5,11); hold on;
for sim = 1:qq; plot(passos, SV_todas(sim,:),  cores{sim}); end
title('SV (Susceptible vectors)');   xlabel('passos'); ylabel('SV'); hold off;

subplot(5,5,12); hold on;
for sim = 1:qq; plot(passos, V1_todas(sim,:),  cores{sim}); end
title('V1 (Vectors Strain 1)');      xlabel('passos'); ylabel('V1'); hold off;

subplot(5,5,13); hold on;
for sim = 1:qq; plot(passos, V2_todas(sim,:),  cores{sim}); end
title('V2 (Vectors Strain 2)');      xlabel('passos'); ylabel('V2'); hold off;

% Subplots de erro (diferença entre sim 1 e sim 2)
subplot(5,5,14); plot(passos, S_todas(1,:)   - S_todas(2,:),   'k');
title('error S');   xlabel('passos'); ylabel('diff S');

subplot(5,5,15); plot(passos, I1_todas(1,:)  - I1_todas(2,:),  'k');
title('error I1');  xlabel('passos'); ylabel('diff I1');

subplot(5,5,16); plot(passos, I2_todas(1,:)  - I2_todas(2,:),  'k');
title('error I2');  xlabel('passos'); ylabel('diff I2');

subplot(5,5,17); plot(passos, R1_todas(1,:)  - R1_todas(2,:),  'k');
title('error R1');  xlabel('passos'); ylabel('diff R1');

subplot(5,5,18); plot(passos, R2_todas(1,:)  - R2_todas(2,:),  'k');
title('error R2');  xlabel('passos'); ylabel('diff R2');

subplot(5,5,19); plot(passos, S1_todas(1,:)  - S1_todas(2,:),  'k');
title('error S1');  xlabel('passos'); ylabel('diff S1');

subplot(5,5,20); plot(passos, S2_todas(1,:)  - S2_todas(2,:),  'k');
title('error S2');  xlabel('passos'); ylabel('diff S2');

subplot(5,5,21); plot(passos, R_todas(1,:)   - R_todas(2,:),   'k');
title('error R');   xlabel('passos'); ylabel('diff R');

subplot(5,5,22); plot(passos, SV_todas(1,:)  - SV_todas(2,:),  'k');
title('error SV');  xlabel('passos'); ylabel('diff SV');

subplot(5,5,23); plot(passos, V1_todas(1,:)  - V1_todas(2,:),  'k');
title('error V1');  xlabel('passos'); ylabel('diff V1');

subplot(5,5,24); plot(passos, V2_todas(1,:)  - V2_todas(2,:),  'k');
title('error V2');  xlabel('passos'); ylabel('diff V2');

print -dpng 'comparacao_simulacoes.png';
fprintf('Gráfico salvo em comparacao_simulacoes.png\n');

% ---------------------------------------------------------------------------
% Estatísticas finais
% ---------------------------------------------------------------------------
fprintf('\n=== ESTATÍSTICAS FINAIS (último passo) ===\n');
variaveis  = {'S','I1','I2','R1','R2','S1','S2','I12','I21','R','SV','V1','V2'};
dados_fin  = {S_todas, I1_todas, I2_todas, R1_todas, R2_todas, S1_todas, S2_todas, ...
              I12_todas, I21_todas, R_todas, SV_todas, V1_todas, V2_todas};
for i = 1:length(variaveis)
  fprintf('\n%s:\n', variaveis{i});
  for sim = 1:qq
    fprintf('  Sim %d: %.2f\n', sim, dados_fin{i}(sim, end));
  end
end
