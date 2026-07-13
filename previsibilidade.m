% ===================================================================
% Modelo epidemiológico com dois sorotipos - versão VETORIZADA (batch)
%
% Análise de sensibilidade: phi é fixo (constante) em todas as amostras
%   Sim 1 (referência): phi = 0.8
%   Sims 2..qq:         phi = 0.8 (mesmo valor, sem variação)
% ===================================================================

% ---------------------------------------------------------------------------
% Campo vetorial escalar (usado na Fase 1) — phi recebido como parâmetro
% ---------------------------------------------------------------------------
inp = zeros(1, 13);
function out = field(inp, t, n, phi)
  s=inp(1); i1=inp(2); i2=inp(3); r1=inp(4); r2=inp(5);
  s1=inp(6); s2=inp(7); i12=inp(8); i21=inp(9);
  r=inp(10); sv=inp(11); v1=inp(12); v2=inp(13);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; omega=2*pi*6;
  xi=nu*(1+0.4*cos(omega*t)); beta=2*gama;

  ds   = -(beta/m)*s*(v1+v2)                   + mu*(n-s);
  di1  =  (beta/m)*s*v1                         - (gama+mu)*i1;
  di2  =  (beta/m)*s*v2                         - (gama+mu)*i2;
  dr1  =  gama*i1                               - (alfa+mu)*r1;
  dr2  =  gama*i2                               - (alfa+mu)*r2;
  ds1  = -(beta/m)*s1*v2  + alfa*r1             - mu*s1;
  ds2  = -(beta/m)*s2*v1  + alfa*r2             - mu*s2;
  di12 =  (beta/m)*s1*v2                        - (gama+mu)*i12;
  di21 =  (beta/m)*s2*v1                        - (gama+mu)*i21;
  dr   =  gama*(i12+i21)                        - mu*r;
  dsv  = -(teta/n)*sv*(i1+i2+phi*(i12+i21))    + xi*m - nu*sv;
  dv1  =  (teta/n)*sv*(i1+phi*i21)             - nu*v1;
  dv2  =  (teta/n)*sv*(i2+phi*i12)             - nu*v2;
  out  = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2];
endfunction

% ---------------------------------------------------------------------------
% Campo vetorial em LOTE — IN é (P × 13), phi é (P × 1), retorna (P × 13)
% Cada linha é uma simulação independente com seu próprio valor de phi.
% ---------------------------------------------------------------------------
function OUT = field_batch(IN, t, N, phi)
  s=IN(:,1); i1=IN(:,2); i2=IN(:,3); r1=IN(:,4); r2=IN(:,5);
  s1=IN(:,6); s2=IN(:,7); i12=IN(:,8); i21=IN(:,9);
  r=IN(:,10); sv=IN(:,11); v1=IN(:,12); v2=IN(:,13);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; omega=2*pi*6;
  xi=nu*(1+0.4*cos(omega*t)); beta=2*gama;

  ds   = -(beta./m).*s.*(v1+v2)                    + mu.*(N-s);
  di1  =  (beta./m).*s.*v1                          - (gama+mu).*i1;
  di2  =  (beta./m).*s.*v2                          - (gama+mu).*i2;
  dr1  =  gama.*i1                                  - (alfa+mu).*r1;
  dr2  =  gama.*i2                                  - (alfa+mu).*r2;
  ds1  = -(beta./m).*s1.*v2  + alfa.*r1             - mu.*s1;
  ds2  = -(beta./m).*s2.*v1  + alfa.*r2             - mu.*s2;
  di12 =  (beta./m).*s1.*v2                         - (gama+mu).*i12;
  di21 =  (beta./m).*s2.*v1                         - (gama+mu).*i21;
  dr   =  gama.*(i12+i21)                           - mu.*r;
  dsv  = -(teta./N).*sv.*(i1+i2+phi.*(i12+i21))    + xi.*m - nu.*sv;
  dv1  =  (teta./N).*sv.*(i1+phi.*i21)             - nu.*v1;
  dv2  =  (teta./N).*sv.*(i2+phi.*i12)             - nu.*v2;

  OUT = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2];
endfunction

% ---------------------------------------------------------------------------
% RK4 escalar — phi repassado como parâmetro
% ---------------------------------------------------------------------------
function out = rk(inp, t, dt, n, phi)
  k1=field(inp,              t,      n, phi);
  k2=field(inp+(dt/2)*k1,    t+dt/2, n, phi);
  k3=field(inp+(dt/2)*k2,    t+dt/2, n, phi);
  k4=field(inp+dt*k3,        t+dt,   n, phi);
  out = max(0, inp + (dt/6)*(k1+2*k2+2*k3+k4));
endfunction

% ---------------------------------------------------------------------------
% RK4 em LOTE — phi é (P × 1), repassado para field_batch
% ---------------------------------------------------------------------------
function OUT = rk_batch(IN, t, dt, N, phi)
  k1=field_batch(IN,             t,      N, phi);
  k2=field_batch(IN+(dt/2).*k1,  t+dt/2, N, phi);
  k3=field_batch(IN+(dt/2).*k2,  t+dt/2, N, phi);
  k4=field_batch(IN+dt.*k3,      t+dt,   N, phi);
  OUT = max(0, IN + (dt/6).*(k1+2.*k2+2.*k3+k4));
endfunction

% ===========================================================================
%  SCRIPT PRINCIPAL
% ===========================================================================

qq  = 100;
dt  = 1/365;
kk  = 200*365;
quantidade_anos = 100;
passos_registrados = quantidade_anos*365;
num_post = kk - passos_registrados;

cond_ini_base = [700,200,100,0,0,0,0,0,0,0,9000,500,500];

% ---------------------------------------------------------------------------
% Definição de phi: agora fixo em 0.8 para todas as simulações do lote
% (linspace mantido apenas por estrutura, mas ambos os limites são iguais,
%  logo phi_lote é um vetor constante = phi_ref = 0.8)
% ---------------------------------------------------------------------------
phi_ref  = 0.8;
phi_lote = linspace(phi_ref * 1, phi_ref * 1, qq - 1)';  % (P×1), constante = 0.8
fprintf('phi_ref = %.4f | phi_lote: %.4f .. %.4f\n', ...
        phi_ref, phi_lote(1), phi_lote(end));

% ---------------------------------------------------------------------------
% FASE 1 — Simulação 1 (escalar, sequencial, phi = phi_ref)
% ---------------------------------------------------------------------------
fprintf('=== FASE 1: Simulação 1 (referência) ===\n');
t1 = tic;

in1   = cond_ini_base;
n0_1  = sum(in1(1:10));
tempo = 0;
I12_ref = zeros(1, num_post);
I21_ref = zeros(1, num_post);
res1    = zeros(num_post, 13);

for k = 1:kk
  in1   = rk(in1, tempo, dt, n0_1, phi_ref);
  tempo = tempo + dt;

  if k > passos_registrados
    idx = k - passos_registrados;
    I12_ref(idx) = in1(8);
    I21_ref(idx) = in1(9);
    res1(idx,:)  = in1;
  end
end

fprintf('Fase 1 concluída em %.2f s\n\n', toc(t1));

% ---------------------------------------------------------------------------
% FASE 2 — Sims 2..qq em lote vetorizado, todas com phi = 0.8 (mesmo valor)
% ---------------------------------------------------------------------------
if qq > 1
  fprintf('=== FASE 2: Simulações 2..%d em lote vetorizado ===\n', qq);
  t2 = tic;

  cond_ini_IN = res1(passos_registrados, :);

  P   = qq - 1;
  IN  = repmat(cond_ini_IN, P, 1);
  N   = sum(IN(:,1:10), 2);
  tempo = 0;

  res_lote = zeros(P, num_post, 13);

  fatores = 0.9999 + 0.0002*rand(P, 13);
  IN = max(0, IN .* fatores);
  fprintf('  Sims 2..%d: perturbações ±0.0001%% aplicadas', qq);

  for k = 1:passos_registrados
    IN    = rk_batch(IN, tempo, dt, N, phi_lote);
    tempo = tempo + dt;

    idx = k;
    res_lote(:, idx, :) = reshape(IN, P, 1, 13);
  end

  fprintf('Fase 2 concluída em %.2f s\n\n', toc(t2));
else
  res_lote = zeros(0, num_post, 13);
end

fprintf('Todas as simulações concluídas!\n');

% ---------------------------------------------------------------------------
% Monta matrizes de resultado (qq × num_post)
% ---------------------------------------------------------------------------
function M = montar(res1_col, res_lote_dim, qq, num_post)
  M = zeros(qq, num_post);
  M(1,:) = res1_col;
  for p = 1:qq-1
    M(p+1,:) = res_lote_dim(p,:);
  end
endfunction

S_todas   = montar(res1(:,1)',  squeeze(res_lote(:,:,1)),  qq, num_post);
I1_todas  = montar(res1(:,2)',  squeeze(res_lote(:,:,2)),  qq, num_post);
I2_todas  = montar(res1(:,3)',  squeeze(res_lote(:,:,3)),  qq, num_post);
R1_todas  = montar(res1(:,4)',  squeeze(res_lote(:,:,4)),  qq, num_post);
R2_todas  = montar(res1(:,5)',  squeeze(res_lote(:,:,5)),  qq, num_post);
S1_todas  = montar(res1(:,6)',  squeeze(res_lote(:,:,6)),  qq, num_post);
S2_todas  = montar(res1(:,7)',  squeeze(res_lote(:,:,7)),  qq, num_post);
I12_todas = montar(res1(:,8)',  squeeze(res_lote(:,:,8)),  qq, num_post);
I21_todas = montar(res1(:,9)',  squeeze(res_lote(:,:,9)),  qq, num_post);
R_todas   = montar(res1(:,10)', squeeze(res_lote(:,:,10)), qq, num_post);
SV_todas  = montar(res1(:,11)', squeeze(res_lote(:,:,11)), qq, num_post);
V1_todas  = montar(res1(:,12)', squeeze(res_lote(:,:,12)), qq, num_post);
V2_todas  = montar(res1(:,13)', squeeze(res_lote(:,:,13)), qq, num_post);

% ---------------------------------------------------------------------------
% Salva resultados
% ---------------------------------------------------------------------------
save('resultados_multiplas_simulacoes.mat', ...
     'S_todas','I1_todas','I2_todas','R1_todas','R2_todas', ...
     'S1_todas','S2_todas','I12_todas','I21_todas','R_todas', ...
     'SV_todas','V1_todas','V2_todas', ...
     'phi_lote','phi_ref');
fprintf('Resultados salvos em resultados_multiplas_simulacoes.mat\n');

% ---------------------------------------------------------------------------
% Plots comparativos — idênticos ao original
% ---------------------------------------------------------------------------
passos   = 1:num_post;
cores    = {'b','r','g','m','c'};
legendas = {'Sim 1','Sim 2','Sim 3','Sim 4','Sim 5'};

figure('Position',[100,100,1400,1000]);

vars_plot = {S_todas,I1_todas,I2_todas,R1_todas,R2_todas, ...
             S1_todas,S2_todas,I12_todas,I21_todas,R_todas, ...
             SV_todas,V1_todas,V2_todas};
titulos   = {'S (Naive Humans)','I1 (Primary Inf Strain 1)', ...
             'I2 (Primary Inf Strain 2)','R1 (Cross Immune from 1)', ...
             'R2 (Cross Immune from 2)','S1 (Susceptible to 2)', ...
             'S2 (Susceptible to 1)','I12 (Secondary Inf Strain 2)', ...
             'I21 (Secondary Inf Strain 1)','R (Totally recovered)', ...
             'SV (Susceptible vectors)','V1 (Vectors Strain 1)', ...
             'V2 (Vectors Strain 2)'};
ylabels   = {'S','I1','I2','R1','R2','S1','S2','I12','I21','R','SV','V1','V2'};

% Variáveis sem par (mantêm 1 janela por variável, com todas as amostras)
mainSingles     = {S_todas, R_todas, SV_todas};
mainSinglesTit  = {'S (Naive Humans)','R (Totally recovered)','SV (Susceptible vectors)'};
mainSinglesYlab = {'S','R','SV'};

% Pares de variáveis: 1 única janela por par, plotando o máximo a cada
% passo sobre TODAS as amostras (qq simulações) de AMBAS as variáveis
mainPairs = { {I1_todas, I2_todas}, {R1_todas, R2_todas}, {S1_todas, S2_todas}, ...
              {I12_todas, I21_todas}, {V1_todas, V2_todas} };
mainPairsTit = {'max{I1,I2} (Infecções Primárias)', ...
                'max{R1,R2} (Imunidade Cruzada)', ...
                'max{S1,S2} (Suscetíveis ao outro sorotipo)', ...
                'max{I12,I21} (Infecções Secundárias)', ...
                'max{V1,V2} (Vetores)'};
mainPairsYlab = {'max(I1,I2)','max(R1,R2)','max(S1,S2)','max(I12,I21)','max(V1,V2)'};

idx_subplot = 0;

for vi = 1:numel(mainSingles)
  idx_subplot = idx_subplot + 1;
  subplot(5,5,idx_subplot); hold on;
  for sim = 1:qq
    plot(passos, mainSingles{vi}(sim,:));
  end
  title(mainSinglesTit{vi}); xlabel('passos'); ylabel(mainSinglesYlab{vi});
  hold off;
end

for vi = 1:numel(mainPairs)
  idx_subplot = idx_subplot + 1;
  par = mainPairs{vi};
  combinado = [par{1}; par{2}];        % (2*qq × num_post): todas as amostras de ambas as variáveis
  max_curva = max(combinado, [], 1);   % máximo por passo, sobre amostras e variáveis
  subplot(5,5,idx_subplot);
  plot(passos, max_curva, 'k');
  title(mainPairsTit{vi}); xlabel('passos'); ylabel(mainPairsYlab{vi});
end

% Subplots de erro — maior erro absoluto entre todas as sims (2..qq) vs sim 1
erros_vars = {S_todas,I1_todas,I2_todas,R1_todas,R2_todas,S1_todas, ...
              S2_todas,R_todas,SV_todas,V1_todas,V2_todas};
erros_tit  = {'max |err| S','max |err| I1','max |err| I2', ...
              'max |err| R1','max |err| R2','max |err| S1', ...
              'max |err| S2','max |err| R','max |err| SV', ...
              'max |err| V1','max |err| V2'};
erros_ylab = {'max |diff| S','max |diff| I1','max |diff| I2', ...
              'max |diff| R1','max |diff| R2','max |diff| S1', ...
              'max |diff| S2','max |diff| R','max |diff| SV', ...
              'max |diff| V1','max |diff| V2'};

for vi = 1:11
  M = erros_vars{vi};
  erros_abs = abs(M(2:end, :) - M(1, :));
  max_erro = max(erros_abs, [], 1);

  subplot(5,5,13+vi);
  plot(passos, max_erro, 'k');
  title(erros_tit{vi}); xlabel('passos'); ylabel(erros_ylab{vi});
end

print -dpng 'comparacao_simulacoes.png';
fprintf('Gráfico salvo em comparacao_simulacoes.png\n');

% ---------------------------------------------------------------------------
% Estatísticas finais
% ---------------------------------------------------------------------------
fprintf('\n=== ESTATÍSTICAS FINAIS (último passo) ===\n');
nomes = {'S','I1','I2','R1','R2','S1','S2','I12','I21','R','SV','V1','V2'};
for i = 1:13
  fprintf('\n%s:\n', nomes{i});
  for sim = 1:qq
    fprintf('  Sim %d: %.2f\n', sim, vars_plot{i}(sim,end));
  end
end
