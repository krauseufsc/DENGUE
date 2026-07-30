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

% ---------------------------------------------------------------------------
% Sistema ESTENDIDO (14 variáveis) para o cálculo do maior expoente de
% Lyapunov pelo método de Benettin.
%
% A 14ª variável, z, é a fase da forçante sazonal: z' = w (constante), com
% xi = nu*(1+0.4*cos(z)) no lugar de nu*(1+0.4*cos(w*t)). Isso torna o
% sistema autônomo (não depende mais de t explicitamente), permitindo medir
% corretamente como perturbações na FASE também se propagam e interagem
% com as demais variáveis — algo que o sistema de 13 variáveis (que usa t
% diretamente) não captura.
% ---------------------------------------------------------------------------
function out = field_lyap(inp, n, phi)
  s=inp(1); i1=inp(2); i2=inp(3); r1=inp(4); r2=inp(5);
  s1=inp(6); s2=inp(7); i12=inp(8); i21=inp(9);
  r=inp(10); sv=inp(11); v1=inp(12); v2=inp(13); z=inp(14);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; w=2*pi*6; beta=2*gama;
  xi=nu*(1+0.4*cos(z));

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
  dz   =  w;
  out  = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2,dz];
endfunction

function OUT = field_lyap_batch(IN, N, phi)
  s=IN(:,1); i1=IN(:,2); i2=IN(:,3); r1=IN(:,4); r2=IN(:,5);
  s1=IN(:,6); s2=IN(:,7); i12=IN(:,8); i21=IN(:,9);
  r=IN(:,10); sv=IN(:,11); v1=IN(:,12); v2=IN(:,13); z=IN(:,14);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; w=2*pi*6; beta=2*gama;
  xi=nu*(1+0.4*cos(z));

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
  dz   =  w * ones(size(z));

  OUT = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2,dz];
endfunction

% RK4 para o sistema estendido — apenas as 13 primeiras componentes são
% forçadas a ser não-negativas (z é uma fase, não deve ser limitada em 0)
function out = rk_lyap(inp, dt, n, phi)
  k1=field_lyap(inp,           n, phi);
  k2=field_lyap(inp+(dt/2)*k1, n, phi);
  k3=field_lyap(inp+(dt/2)*k2, n, phi);
  k4=field_lyap(inp+dt*k3,     n, phi);
  out = inp + (dt/6)*(k1+2*k2+2*k3+k4);
  out(1:13) = max(0, out(1:13));
endfunction

function OUT = rk_lyap_batch(IN, dt, N, phi)
  k1=field_lyap_batch(IN,            N, phi);
  k2=field_lyap_batch(IN+(dt/2).*k1, N, phi);
  k3=field_lyap_batch(IN+(dt/2).*k2, N, phi);
  k4=field_lyap_batch(IN+dt.*k3,     N, phi);
  OUT = IN + (dt/6).*(k1+2.*k2+2.*k3+k4);
  OUT(:,1:13) = max(0, OUT(:,1:13));
endfunction

% ---------------------------------------------------------------------------
% Produto interno e norma PONDERADOS (não-euclidianos), usados para medir o
% tamanho da perturbação entre a trajetória de referência e a perturbada.
% As 10 variáveis humanas são adimensionalizadas por dS, as 3 variáveis de
% mosquito por dSv, e a fase z por dZ — sem essa ponderação, a norma seria
% dominada pelas variáveis de maior magnitude (populações, na casa dos
% milhares) e a contribuição da fase (em radianos) seria irrelevante.
% Versão em lote: V1, V2 são matrizes (P × 14), retorno é (P × 1).
% ---------------------------------------------------------------------------
function out = produ_batch(V1, V2, dS, dSv, dZ)
  out =       sum(V1(:,1:10)  .* V2(:,1:10),  2) / dS^2;
  out = out + sum(V1(:,11:13) .* V2(:,11:13), 2) / dSv^2;
  out = out +      (V1(:,14)  .* V2(:,14))       / dZ^2;
endfunction

function out = norma_batch(V, dS, dSv, dZ)
  out = sqrt(produ_batch(V, V, dS, dSv, dZ));
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
% FASE 3 — Maior expoente de Lyapunov (método de Benettin, sistema completo)
%
% Duas trajetórias por amostra: X (referência) e Y (perturbada), ambas no
% sistema estendido de 14 variáveis (13 originais + fase z). A cada
% intervalo de normalização (1 mês = 30 passos), mede-se a separação
% ponderada entre X e Y, acumula-se ln(d/eps), e Y é reposicionada de volta
% à distância eps de X (mesma direção), evitando que a perturbação saia do
% regime linear. No fim: lambda_max ≈ soma / (N_normalizações * DT).
% ---------------------------------------------------------------------------
fprintf('=== FASE 3: Expoente de Lyapunov (Benettin, 14 variáveis) ===\n');
t3 = tic;

DT_passos = 30;                     % intervalo entre normalizações: 1 mês = 30 passos
DT_tempo  = DT_passos * dt;         % o mesmo intervalo, em anos
w_forcante = 2*pi*6;

eps_lyap = 1e-4;                    % mesma magnitude da perturbação da Fase 2 (0.0002*rand, meia-largura 1e-4)

% Ponto de partida: mesmo checkpoint pós-burn-in usado na Fase 2
X0_13 = res1(passos_registrados, :);                       % 1×13
z0    = mod(w_forcante * (passos_registrados*dt), 2*pi);   % fase acumulada até o checkpoint
X0    = [X0_13, z0];                                        % 1×14

dS_lyap  = n0_1;                    % população humana total (constante)
dSv_lyap = sum(X0_13(11:13));       % população de mosquito no checkpoint (referência de escala)
dZ_lyap  = w_forcante * DT_tempo;   % escala de fase entre normalizações

P_lyap = qq - 1;                    % mesmo número de amostras perturbadas da Fase 2

% Direções iniciais de perturbação, normalizadas para ||.|| = eps_lyap
DIR0   = randn(P_lyap, 14);
DIR0   = DIR0 ./ norma_batch(DIR0, dS_lyap, dSv_lyap, dZ_lyap);
DELTA0 = eps_lyap * DIR0;

% X_lyap: UMA única trajetória de referência (1×14) — todas as amostras
% comparam-se com a MESMA referência, então ela não precisa ser repetida.
% Y_lyap: P trajetórias perturbadas (P×14), uma por amostra.
X_lyap = X0;
Y_lyap = X0 + DELTA0;   % soma broadcast: (1×14) + (P×14) = (P×14)

n_normalizacoes = floor(num_post / DT_passos);
soma_log = zeros(P_lyap, 1);

for bloco = 1:n_normalizacoes
  for passo = 1:DT_passos
    X_lyap = rk_lyap(X_lyap, dt, dS_lyap, phi_ref);         % 1 trajetória escalar
    Y_lyap = rk_lyap_batch(Y_lyap, dt, dS_lyap, phi_ref);   % P trajetórias em lote
  end

  DELTA   = Y_lyap - X_lyap;                     % broadcast: (P×14) − (1×14)
  d_norma = norma_batch(DELTA, dS_lyap, dSv_lyap, dZ_lyap);
  soma_log = soma_log + log(d_norma / eps_lyap);

  Y_lyap = X_lyap + eps_lyap * (DELTA ./ d_norma);   % renormalização (broadcast)
end

lambda_max_amostras = soma_log / (n_normalizacoes * DT_tempo);   % P×1
lambda_max_medio    = mean(lambda_max_amostras);

fprintf('Fase 3 concluída em %.2f s\n\n', toc(t3));
fprintf('=== MAIOR EXPOENTE DE LYAPUNOV (Benettin) ===\n');
fprintf('  Intervalo de normalização: %d passos (%.6f anos, ~1 mês)\n', DT_passos, DT_tempo);
fprintf('  eps = %.1e | %d amostras | %d normalizações\n', eps_lyap, P_lyap, n_normalizacoes);
for p = 1:P_lyap
  fprintf('  Amostra %d: lambda_max = %+.6f (1/ano)\n', p, lambda_max_amostras(p));
end
fprintf('  Média sobre as amostras: lambda_max = %+.6f (1/ano)\n\n', lambda_max_medio);

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
% Salva resultados — organizados em structs, para facilitar a leitura
% (em vez de dezenas de variáveis soltas no .mat)
%
%   parametros  — configuração da simulação (qq, dt, phi, condições iniciais)
%   simulacoes  — as 13 variáveis do modelo, cada uma (qq × num_post)
%   lyapunov    — configuração e resultado do método de Benettin (Fase 3)
% ---------------------------------------------------------------------------
parametros.qq                  = qq;
parametros.dt                  = dt;
parametros.kk                  = kk;
parametros.quantidade_anos     = quantidade_anos;
parametros.passos_registrados  = passos_registrados;
parametros.num_post            = num_post;
parametros.phi_ref             = phi_ref;
parametros.phi_lote            = phi_lote;
parametros.cond_ini_base       = cond_ini_base;

simulacoes.S   = S_todas;
simulacoes.I1  = I1_todas;
simulacoes.I2  = I2_todas;
simulacoes.R1  = R1_todas;
simulacoes.R2  = R2_todas;
simulacoes.S1  = S1_todas;
simulacoes.S2  = S2_todas;
simulacoes.I12 = I12_todas;
simulacoes.I21 = I21_todas;
simulacoes.R   = R_todas;
simulacoes.SV  = SV_todas;
simulacoes.V1  = V1_todas;
simulacoes.V2  = V2_todas;

lyapunov.DT_passos            = DT_passos;
lyapunov.DT_tempo             = DT_tempo;
lyapunov.eps                  = eps_lyap;
lyapunov.dS                   = dS_lyap;
lyapunov.dSv                  = dSv_lyap;
lyapunov.dZ                   = dZ_lyap;
lyapunov.n_normalizacoes      = n_normalizacoes;
lyapunov.lambda_max_amostras  = lambda_max_amostras;
lyapunov.lambda_max_medio     = lambda_max_medio;

save('resultados_multiplas_simulacoes.mat', 'parametros', 'simulacoes', 'lyapunov');
fprintf('Resultados salvos (organizados em structs) em resultados_multiplas_simulacoes.mat\n');

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

idx_subplot = 0;

% Subplots de erro — maior erro absoluto entre todas as sims (2..qq) vs sim 1
% Variáveis sem par: erro individual (máximo sobre as amostras)
errSingles     = {S_todas, R_todas, SV_todas};
errSinglesTit  = {'max |err| S','max |err| R','max |err| SV'};
errSinglesYlab = {'max |diff| S','max |diff| R','max |diff| SV'};

% Pares de variáveis de cepas diferentes: 1 única janela por par, com o
% máximo do erro absoluto (cada variável comparada com sua própria sim 1)
% sobre TODAS as amostras (qq-1 sims) de AMBAS as variáveis.
errPairs = { {I1_todas, I2_todas}, {R1_todas, R2_todas}, {S1_todas, S2_todas}, ...
             {I12_todas, I21_todas}, {V1_todas, V2_todas} };
errPairsTit  = {'max |err| I1,I2','max |err| R1,R2','max |err| S1,S2', ...
                'max |err| I12,I21','max |err| V1,V2'};
errPairsYlab = {'max |diff|(I1,I2)','max |diff|(R1,R2)','max |diff|(S1,S2)', ...
                'max |diff|(I12,I21)','max |diff|(V1,V2)'};

for vi = 1:numel(errSingles)
  idx_subplot = idx_subplot + 1;
  M = errSingles{vi};
  erros_abs = abs(M(2:end, :) - M(1, :));
  max_erro  = max(erros_abs, [], 1);

  subplot(5,5,idx_subplot);
  plot(passos, max_erro, 'k');
  title(errSinglesTit{vi}); xlabel('passos'); ylabel(errSinglesYlab{vi});
end

for vi = 1:numel(errPairs)
  idx_subplot = idx_subplot + 1;
  par = errPairs{vi};
  erro1 = abs(par{1}(2:end, :) - par{1}(1, :));   % erro de x1 vs sua sim 1
  erro2 = abs(par{2}(2:end, :) - par{2}(1, :));   % erro de x2 vs sua sim 1
  erros_abs = [erro1; erro2];                     % (2*(qq-1) × num_post)
  max_erro  = max(erros_abs, [], 1);              % máximo sobre amostras e variáveis

  subplot(5,5,idx_subplot);
  plot(passos, max_erro, 'k');
  title(errPairsTit{vi}); xlabel('passos'); ylabel(errPairsYlab{vi});
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
