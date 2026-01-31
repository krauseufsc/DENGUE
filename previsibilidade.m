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

function out = rk(inp, t, dt, n)
    k1 = field(inp, t, n);
    k2 = field(inp + (dt / 2) * k1, t + dt / 2, n);
    k3 = field(inp + (dt / 2) * k2, t + dt / 2, n);
    k4 = field(inp + dt * k3, t + dt, n);

    temp = k1 + 2 * k2 + 2 * k3 + k4;
    out = max(0, inp + (dt / 6) * temp);
endfunction

% Main script
in = zeros(1, 13);
in(1) = 700;
in(2) = 200;
in(3) = 100;
in(4) = 0;
in(5) = 0;
in(6) = 0;
in(7) = 0;
in(8) = 0;
in(9) = 0;
in(10) = 0;
in(11) = 9000;
in(12) = 500;
in(13) = 500;

tempo = 0;
dt = 1 / 365;
acum = 0;
for i = 1:10
    acum += in(i);
endfor
n0 = acum;  % numero de pessoas

kk = 100 * 365;

% cria listas para inserir os resultados de cada variável passo a passo
S = zeros(1,kk);
I1 = zeros(1,kk);
I2 = zeros(1,kk);
R1 = zeros(1,kk);
R2 = zeros(1,kk);
S1 = zeros(1,kk);
S2 = zeros(1,kk);
I12 = zeros(1,kk);
I21 = zeros(1,kk);
R = zeros(1,kk);
SV = zeros(1,kk);
V1 = zeros(1,kk);
V2 = zeros(1,kk);

 for k = 1:kk
     in = rk(in, tempo, dt, n0);
     tempo = tempo + dt;
     %insere os resultados de cada variável em uma lista para depois plotar
     S(k) = in(1);
     I1(k) = in(2);
     I2(k) = in(3);
     R1(k) = in(4);
     R2(k) = in(5);
     S1(k) = in(6);
     S2(k) = in(7);
     I12(k) = in(8);
     I21(k) = in(9);
     R(k) = in(10);
     SV(k) = in(11);
     V1(k) = in(12);
     V2(k) = in(13);
     % Para kk passos

endfor

%plots

passos = 1:kk;

figure

subplot(4,4,1)
plot(passos, S);
title ('S (Naive Humans)');
xlabel ('passos');
ylabel ('S');

subplot(4,4,2)
plot(passos, I1);
title('I1 (Primary Inf Strain 1)');
xlabel ('passos');
ylabel ('I1');

subplot(4,4,3)
plot(passos, I2);
title('I2 (Primary Inf Strain 2)');
xlabel ('passos');
ylabel ('I2');

subplot(4,4,4)
plot(passos, R1);
title('R1 (Cross Immune from 1)');
xlabel ('passos');
ylabel ('R1');

subplot(4,4,5)
plot(passos, R2);
title('R2 (Cross Immune from 2)');
xlabel ('passos');
ylabel ('R2');

subplot(4,4,6)
plot(passos, S1);
title('S1 (Susceptible to 2)');
xlabel ('passos');
ylabel ('S1');

subplot(4,4,7)
plot(passos, S2);
title('S2 (Susceptible to 1)');
xlabel ('passos');
ylabel ('S2');

subplot(4,4,8)
plot(passos, I12);
title('Secondary Inf Strain 2');
xlabel ('passos');
ylabel ('I12');

subplot(4,4,9)
plot(passos, I21);
title('I21 (Secondary Inf Strain 1)');
xlabel ('passos');
ylabel ('I21');

subplot(4,4,10)
plot(passos, R);
title('R (Totally recovered)');
xlabel ('passos');
ylabel ('R');

subplot(4,4,11)
plot(passos, SV);
title('SV (Susceptible vectors)');
xlabel ('passos');
ylabel ('SV');

subplot(4,4,12)
plot(passos, V1);
title('V1 (Vectors Strain 1)');
xlabel ('passos');
ylabel ('V1');

subplot(4,4,13)
plot(passos, V2);
title('V2 (Vectors Strain 2)');
xlabel ('passos');
ylabel ('V2');

result = rk(in, tempo, dt, n0);
disp(result);

