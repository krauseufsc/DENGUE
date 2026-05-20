% ===================================================================
% Modelo epidemiológico com dois sorotipos - versão VETORIZADA (batch)
%
% Estratégia de paralelização:
%   Em vez de usar threads/processos (overhead alto no Octave), todas as
%   simulações rodam SIMULTANEAMENTE em uma única matriz (qq × 13).
%   Cada linha é uma simulação; as operações são vetorizadas via BLAS.
%
%   Fase 1 — Sim 1 corre isolada para produzir I12_ref e I21_ref.
%   Fase 2 — Sims 2..qq correm juntas em batch: estado = matriz (qq-1)×13.
%             A cada passo, colunas 8-9 das linhas dependentes são
%             substituídas por I12_ref/I21_ref (leitura, sem concorrência).
% ===================================================================

% ---------------------------------------------------------------------------
% Campo vetorial escalar (usado na Fase 1)
% ---------------------------------------------------------------------------
inp = zeros(1, 13);
function out = field(inp, t, n)
  s=inp(1); i1=inp(2); i2=inp(3); r1=inp(4); r2=inp(5);
  s1=inp(6); s2=inp(7); i12=inp(8); i21=inp(9);
  r=inp(10); sv=inp(11); v1=inp(12); v2=inp(13);

  m=sv+v1+v2; phi=0.8; mu=1/65; alfa=2; gama=52;
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
% Campo vetorial em LOTE — IN é (P × 13), retorna (P × 13)
% Cada linha é uma simulação independente; operações element-wise via BLAS.
% ---------------------------------------------------------------------------
function OUT = field_batch(IN, t, N)
  % N : vetor coluna (P×1) com população total de cada simulação
  s=IN(:,1); i1=IN(:,2); i2=IN(:,3); r1=IN(:,4); r2=IN(:,5);
  s1=IN(:,6); s2=IN(:,7); i12=IN(:,8); i21=IN(:,9);
  r=IN(:,10); sv=IN(:,11); v1=IN(:,12); v2=IN(:,13);

  m=sv+v1+v2; phi=0.8; mu=1/65; alfa=2; gama=52;
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
% RK4 escalar
% ---------------------------------------------------------------------------
function out = rk(inp, t, dt, n)
  k1=field(inp,t,n);
  k2=field(inp+(dt/2)*k1, t+dt/2, n);
  k3=field(inp+(dt/2)*k2, t+dt/2, n);
  k4=field(inp+dt*k3,      t+dt,   n);
  out = max(0, inp + (dt/6)*(k1+2*k2+2*k3+k4));
endfunction

% ---------------------------------------------------------------------------
% RK4 em LOTE — IN (P×13), retorna (P×13)
% ---------------------------------------------------------------------------
function OUT = rk_batch(IN, t, dt, N)
  k1=field_batch(IN,             t,        N);
  k2=field_batch(IN+(dt/2).*k1,  t+dt/2,   N);
  k3=field_batch(IN+(dt/2).*k2,  t+dt/2,   N);
  k4=field_batch(IN+dt.*k3,       t+dt,    N);
  OUT = max(0, IN + (dt/6).*(k1+2.*k2+2.*k3+k4));
endfunction

% ===========================================================================
%  SCRIPT PRINCIPAL
% ===========================================================================

qq  = 5;
dt  = 1/365;
kk  = 100*365;
num_post = kk - 18250;

cond_ini_base = [700,200,100,0,0,0,0,0,0,0,9000,500,500];

% ---------------------------------------------------------------------------
% FASE 1 — Simulação 1 (escalar, sequencial)
%   Obrigatória: gera I12_ref e I21_ref para as demais.
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
  in1   = rk(in1, tempo, dt, n0_1);
  tempo = tempo + dt;

  if k == 18250
    fatores = 0.95 + 0.10*rand(1,13);
    in1 = max(0, in1 .* fatores);
    fprintf('  Sim 1: perturbação ±5%% aplicada no passo 18250\n');
  end

  if k > 18250
    idx = k - 18250;
    I12_ref(idx) = in1(8);
    I21_ref(idx) = in1(9);
    res1(idx,:)  = in1;
  end
end

fprintf('Fase 1 concluída em %.2f s\n\n', toc(t1));

% ---------------------------------------------------------------------------
% FASE 2 — Sims 2..qq em LOTE vetorizado
%   Estado IN: matriz (P × 13), P = qq-1
%   Todas as simulações avançam juntas num único rk_batch por passo.
%   A cada passo após 18250, colunas 8-9 são sobrepostas por I12/I21 da ref.
% ---------------------------------------------------------------------------
if qq > 1
  fprintf('=== FASE 2: Simulações 2..%d em lote vetorizado ===\n', qq);
  t2 = tic;

  P   = qq - 1;                             % número de sims em lote
  IN  = repmat(cond_ini_base, P, 1);        % condições iniciais iguais
  N   = sum(IN(:,1:10), 2);                 % pop. total por simulação (col)
  tempo = 0;

  % Pré-aloca resultados do lote: cell array (uma célula por sim)
  res_lote = zeros(P, num_post, 13);        % (sims × passos × vars)

  for k = 1:kk
    IN    = rk_batch(IN, tempo, dt, N);
    tempo = tempo + dt;

    if k == 18250
      % Perturbação independente para cada simulação do lote
      fatores = 0.95 + 0.10*rand(P, 13);
      IN = max(0, IN .* fatores);
      fprintf('  Sims 2..%d: perturbações ±5%% aplicadas no passo 18250\n', qq);
    end

    if k > 18250
      idx = k - 18250;
      % Força I12/I21 da referência em todas as sims do lote
      IN(:,8) = I12_ref(idx);
      IN(:,9) = I21_ref(idx);
      res_lote(:, idx, :) = reshape(IN, P, 1, 13);
    end
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
     'SV_todas','V1_todas','V2_todas');
fprintf('Resultados salvos em resultados_multiplas_simulacoes.mat\n');

% ---------------------------------------------------------------------------
% Plots comparativos
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

for vi = 1:13
  subplot(5,5,vi); hold on;
  for sim = 1:qq
    plot(passos, vars_plot{vi}(sim,:), cores{sim});
  end
  title(titulos{vi}); xlabel('passos'); ylabel(ylabels{vi});
  if vi == 1
    legend(legendas(1:qq), 'Location','best','FontSize',6);
  end
  hold off;
end

% Subplots de erro — maior erro absoluto entre todas as sims (2..qq) vs sim 1
% Para cada variável e cada passo: max_sim( |var(sim,:) - var(1,:)| )
erros_vars = {S_todas,I1_todas,I2_todas,R1_todas,R2_todas,S1_todas, ...
              S2_todas,R_todas,SV_todas,V1_todas,V2_todas};
erros_tit  = {'max err S','max err I1','max err I2', ...
              'max err R1','max err R2','max err S1', ...
              'max err S2','max err R','max err SV', ...
              'max err V1','max err V2'};
erros_ylab = {'max |diff| S','max |diff| I1','max |diff| I2', ...
              'max |diff| R1','max |diff| R2','max |diff| S1', ...
              'max |diff| S2','max |diff| R','max |diff| SV', ...
              'max |diff| V1','max |diff| V2'};

for vi = 1:11
  M = erros_vars{vi};
  % abs. error de cada sim (2..qq) em relação à sim 1 → matriz (qq-1) × num_post
  erros_abs = abs(M(2:end, :) - M(1, :));
  % maior valor em cada passo de tempo → vetor 1 × num_post
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
