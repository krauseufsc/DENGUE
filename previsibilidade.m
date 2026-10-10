function previsibilidade()
    % ===================================================================
    % Cálculo do Espectro de Lyapunov (Método de Benettin com Gram-Schmidt)
    % Integração Tangente Acoplada com Diferenças Finitas Centrais
    % ===================================================================

    dt = 1/365;
    anos_burnin = 100;
    passos_burnin = anos_burnin * 365;
    anos_lyap = 150;
    num_post = anos_lyap * 365;

    cond_ini_base = [700,200,100,0,0,0,0,0,0,0,9000,500,500];
    phi_ref = 0.8;

    % ---------------------------------------------------------------------------
    % FASE 1 — Burn-in (Atingir o atrator)
    % ---------------------------------------------------------------------------
    fprintf('=== FASE 1: Burn-in (%d anos) ===\n', anos_burnin);
    t1 = tic;

    in1  = cond_ini_base;
    n0_1 = sum(in1(1:10));
    tempo = 0;

    for k = 1:passos_burnin
      in1 = rk_burnin(in1, tempo, dt, n0_1, phi_ref);
      tempo = tempo + dt;
    end
    X0_13 = in1;

    fprintf('Burn-in concluído em %.2f s\n\n', toc(t1));

    % ---------------------------------------------------------------------------
    % FASE 2 — Espectro de Lyapunov (Benettin-Gram-Schmidt Direcional)
    % ---------------------------------------------------------------------------
    fprintf('=== FASE 2: Espectro de Lyapunov ===\n');
    t2 = tic;

    DT_passos = 30;
    DT_tempo  = DT_passos * dt;
    w_forcante = 2*pi*6;

    num_exp = 14;
    anos_burnin_pert = 60;

    z0 = mod(w_forcante * (passos_burnin * dt), 2*pi);
    X_lyap = [X0_13, z0];

    dS_lyap  = n0_1;
    dSv_lyap = sum(X0_13(11:13));
    dZ_lyap  = w_forcante * DT_tempo;

    % Inicialização da base ortonormal com DUPLA ortogonalização
    W_mat = randn(num_exp, 14);
    for j = 1:num_exp
        U_j = W_mat(j, :);
        for i = 1:(j-1)
            proj = produ_batch(U_j, W_mat(i, :), dS_lyap, dSv_lyap, dZ_lyap);
            U_j  = U_j - proj .* W_mat(i, :);
        end
        for i = 1:(j-1) % Segunda passada garante ortogonalidade perfeita
            proj = produ_batch(U_j, W_mat(i, :), dS_lyap, dSv_lyap, dZ_lyap);
            U_j  = U_j - proj .* W_mat(i, :);
        end
        W_mat(j, :) = U_j ./ norma_batch(U_j, dS_lyap, dSv_lyap, dZ_lyap);
    end

    n_norm_total   = floor(num_post / DT_passos);
    n_norm_burnin  = floor((anos_burnin_pert * 365) / DT_passos);
    n_norm_medicao = n_norm_total - n_norm_burnin;

    soma_log = zeros(1, num_exp);

    for bloco = 1:n_norm_total

        for passo = 1:DT_passos
            [X_lyap, W_mat] = rk_lyap_tangent(X_lyap, W_mat, dt, dS_lyap, phi_ref);
        end

        W_novo = zeros(num_exp, 14);

        for j = 1:num_exp
            U_j = W_mat(j, :);

            % Dupla Ortogonalização (Previne que W_2 herde energia de W_1)
            % Passada 1
            for i = 1:(j-1)
                W_i = W_novo(i, :);
                proj = produ_batch(U_j, W_i, dS_lyap, dSv_lyap, dZ_lyap);
                U_j  = U_j - proj .* W_i;
            end
            % Passada 2 (Refinamento de Kahan)
            for i = 1:(j-1)
                W_i = W_novo(i, :);
                proj = produ_batch(U_j, W_i, dS_lyap, dSv_lyap, dZ_lyap);
                U_j  = U_j - proj .* W_i;
            end

            alpha_j = norma_batch(U_j, dS_lyap, dSv_lyap, dZ_lyap);

            if bloco > n_norm_burnin
                soma_log(j) = soma_log(j) + log(alpha_j);
            end

            W_novo(j, :) = U_j ./ alpha_j;
        end

        W_mat = W_novo;
    end

    lambda_medio = soma_log / (n_norm_medicao * DT_tempo);

    fprintf('Cálculo de Lyapunov concluído em %.2f s\n\n', toc(t2));
    fprintf('=== RESULTADOS DO ESPECTRO ===\n');
    fprintf('  Burn-in descartado: %d anos (%d normalizações)\n', anos_burnin_pert, n_norm_burnin);
    fprintf('  Tempo de medição real: %d anos (%d normalizações)\n', (anos_lyap - anos_burnin_pert), n_norm_medicao);
    fprintf('  Espectro calculado (1/ano):\n');
    for j = 1:num_exp
        fprintf('    lambda_%d = %+.6f\n', j, lambda_medio(j));
    end
    fprintf('\n');

    % ---------------------------------------------------------------------------
    % Salvamento dos Resultados
    % ---------------------------------------------------------------------------
    parametros.dt              = dt;
    parametros.anos_burnin     = anos_burnin;
    parametros.anos_lyap       = anos_lyap;
    parametros.phi_ref         = phi_ref;
    parametros.cond_ini_base   = cond_ini_base;

    lyapunov.DT_passos       = DT_passos;
    lyapunov.DT_tempo        = DT_tempo;
    lyapunov.dS              = dS_lyap;
    lyapunov.dSv             = dSv_lyap;
    lyapunov.dZ              = dZ_lyap;
    lyapunov.n_norm_total = n_norm_total;
    lyapunov.num_exp         = num_exp;
    lyapunov.lambda_medio    = lambda_medio;

    save('resultados_lyapunov.mat', 'parametros', 'lyapunov');
    fprintf('Resultados estruturados salvos em resultados_lyapunov.mat\n');

end


% ===========================================================================
%  FUNÇÕES AUXILIARES
% ===========================================================================

function out = field_burnin(inp, t, n, phi)
  s=inp(1); i1=inp(2); i2=inp(3); r1=inp(4); r2=inp(5);
  s1=inp(6); s2=inp(7); i12=inp(8); i21=inp(9);
  r=inp(10); sv=inp(11); v1=inp(12); v2=inp(13);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; omega=2*pi*6;
  xi=nu*(1+0.4*cos(omega*t)); beta=2*gama;

  ds   = -(beta/m)*s*(v1+v2)                   + mu*(n-s);
  di1  =  (beta/m)*s*v1                          - (gama+mu)*i1;
  di2  =  (beta/m)*s*v2                          - (gama+mu)*i2;
  dr1  =  gama*i1                                - (alfa+mu)*r1;
  dr2  =  gama*i2                                - (alfa+mu)*r2;
  ds1  = -(beta/m)*s1*v2  + alfa*r1              - mu*s1;
  ds2  = -(beta/m)*s2*v1  + alfa*r2              - mu*s2;
  di12 =  (beta/m)*s1*v2                         - (gama+mu)*i12;
  di21 =  (beta/m)*s2*v1                         - (gama+mu)*i21;
  dr   =  gama*(i12+i21)                         - mu*r;
  dsv  = -(teta/n)*sv*(i1+i2+phi*(i12+i21))      + xi*m - nu*sv;
  dv1  =  (teta/n)*sv*(i1+phi*i21)               - nu*v1;
  dv2  =  (teta/n)*sv*(i2+phi*i12)               - nu*v2;
  out  = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2];
end

function out = rk_burnin(inp, t, dt, n, phi)
  k1=field_burnin(inp,              t,      n, phi);
  k2=field_burnin(inp+(dt/2)*k1,    t+dt/2, n, phi);
  k3=field_burnin(inp+(dt/2)*k2,    t+dt/2, n, phi);
  k4=field_burnin(inp+dt*k3,        t+dt,   n, phi);
  out = max(0, inp + (dt/6)*(k1+2*k2+2*k3+k4));
end

function out = field_lyap(inp, n, phi)
  s=inp(1); i1=inp(2); i2=inp(3); r1=inp(4); r2=inp(5);
  s1=inp(6); s2=inp(7); i12=inp(8); i21=inp(9);
  r=inp(10); sv=inp(11); v1=inp(12); v2=inp(13); z=inp(14);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; w=2*pi*6; beta=2*gama;
  xi=nu*(1+0.4*cos(z));

  ds   = -(beta/m)*s*(v1+v2)                   + mu*(n-s);
  di1  =  (beta/m)*s*v1                          - (gama+mu)*i1;
  di2  =  (beta/m)*s*v2                          - (gama+mu)*i2;
  dr1  =  gama*i1                                - (alfa+mu)*r1;
  dr2  =  gama*i2                                - (alfa+mu)*r2;
  ds1  = -(beta/m)*s1*v2  + alfa*r1              - mu*s1;
  ds2  = -(beta/m)*s2*v1  + alfa*r2              - mu*s2;
  di12 =  (beta/m)*s1*v2                         - (gama+mu)*i12;
  di21 =  (beta/m)*s2*v1                         - (gama+mu)*i21;
  dr   =  gama*(i12+i21)                         - mu*r;
  dsv  = -(teta/n)*sv*(i1+i2+phi*(i12+i21))      + xi*m - nu*sv;
  dv1  =  (teta/n)*sv*(i1+phi*i21)               - nu*v1;
  dv2  =  (teta/n)*sv*(i2+phi*i12)               - nu*v2;
  dz   =  w;
  out  = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2,dz];
end

function OUT = field_lyap_batch(IN, N, phi)
  s=IN(:,1); i1=IN(:,2); i2=IN(:,3); r1=IN(:,4); r2=IN(:,5);
  s1=IN(:,6); s2=IN(:,7); i12=IN(:,8); i21=IN(:,9);
  r=IN(:,10); sv=IN(:,11); v1=IN(:,12); v2=IN(:,13); z=IN(:,14);

  m=sv+v1+v2; mu=1/65; alfa=2; gama=52;
  nu=36.5; teta=2*nu; w=2*pi*6; beta=2*gama;
  xi=nu*(1+0.4*cos(z));

  ds   = -(beta./m).*s.*(v1+v2)                   + mu.*(N-s);
  di1  =  (beta./m).*s.*v1                          - (gama+mu).*i1;
  di2  =  (beta./m).*s.*v2                          - (gama+mu).*i2;
  dr1  =  gama.*i1                                  - (alfa+mu).*r1;
  dr2  =  gama.*i2                                  - (alfa+mu).*r2;
  ds1  = -(beta./m).*s1.*v2  + alfa.*r1             - mu.*s1;
  ds2  = -(beta./m).*s2.*v1  + alfa.*r2             - mu.*s2;
  di12 =  (beta./m).*s1.*v2                         - (gama+mu).*i12;
  di21 =  (beta./m).*s2.*v1                         - (gama+mu).*i21;
  dr   =  gama.*(i12+i21)                           - mu.*r;
  dsv  = -(teta./N).*sv.*(i1+i2+phi.*(i12+i21))     + xi.*m - nu.*sv;
  dv1  =  (teta./N).*sv.*(i1+phi.*i21)              - nu.*v1;
  dv2  =  (teta./N).*sv.*(i2+phi.*i12)              - nu.*v2;
  dz   =  w * ones(size(z));

  OUT = [ds,di1,di2,dr1,dr2,ds1,ds2,di12,di21,dr,dsv,dv1,dv2,dz];
end

% ---------------------------------------------------------------------------
% INTEGRADOR TANGENTE COM JACOBIANA ANALÍTICA (Octave/MATLAB)
% ---------------------------------------------------------------------------
function [X_out, W_out] = rk_lyap_tangent(X, W, dt, N, phi)

  % Converte o vetor X em uma lista separada por vírgulas para a Jacobiana
  args_X = num2cell(X);

  % Passo 1
  F_X1 = field_lyap(X, N, phi);
  J1   = jacobiana_lyap(args_X{:}, N, phi);
  k1_X = F_X1;
  k1_W = W * J1';

  % Passo 2
  X2 = X + (dt/2)*k1_X;
  X2(1:13) = max(0, X2(1:13));
  W2 = W + (dt/2)*k1_W;

  args_X2 = num2cell(X2);
  F_X2 = field_lyap(X2, N, phi);
  J2   = jacobiana_lyap(args_X2{:}, N, phi);
  k2_X = F_X2;
  k2_W = W2 * J2';

  % Passo 3
  X3 = X + (dt/2)*k2_X;
  X3(1:13) = max(0, X3(1:13));
  W3 = W + (dt/2)*k2_W;

  args_X3 = num2cell(X3);
  F_X3 = field_lyap(X3, N, phi);
  J3   = jacobiana_lyap(args_X3{:}, N, phi);
  k3_X = F_X3;
  k3_W = W3 * J3';

  % Passo 4
  X4 = X + dt*k3_X;
  X4(1:13) = max(0, X4(1:13));
  W4 = W + dt*k3_W;

  args_X4 = num2cell(X4);
  F_X4 = field_lyap(X4, N, phi);
  J4   = jacobiana_lyap(args_X4{:}, N, phi);
  k4_X = F_X4;
  k4_W = W4 * J4';

  % Combinação Final
  X_out = X + (dt/6)*(k1_X + 2*k2_X + 2*k3_X + k4_X);
  X_out(1:13) = max(0, X_out(1:13));

  W_out = W + (dt/6)*(k1_W + 2*k2_W + 2*k3_W + k4_W);
end

function out = produ_batch(V1, V2, dS, dSv, dZ)
  out =        sum(V1(:,1:10)  .* V2(:,1:10),  2) / dS^2;
  out = out + sum(V1(:,11:13) .* V2(:,11:13), 2) / dSv^2;
  out = out +      (V1(:,14)  .* V2(:,14))       / dZ^2;
end

function out = norma_batch(V, dS, dSv, dZ)
  out = sqrt(produ_batch(V, V, dS, dSv, dZ));
end
