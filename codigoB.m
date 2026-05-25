%% PARTE B — INTERPOLACIÓN COMPLETA
% B1: Polinómica (Matricial + Lagrange) | B2: Splines cúbicos

clc; clear; close all;

%% =========================================================
%% DATOS
%% =========================================================
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, ...
     460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, ...
     1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];

Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, ...
     132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, ...
     135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

n  = length(f);
f_eval = 1000;           % punto de evaluación pedido
f_fine = linspace(min(f), max(f), 2000);   % malla fina para graficar

%% =========================================================
%% B1 — INTERPOLACIÓN POLINÓMICA MATRICIAL (grado 29)
%% =========================================================
% Normalizar frecuencias para mejorar condicionamiento
f_norm = (f - mean(f)) / std(f);
f_fine_norm = (f_fine - mean(f)) / std(f);
f_eval_norm = (f_eval - mean(f)) / std(f);

% Construir matriz de Vandermonde manualmente
grado29 = n - 1;
V29 = zeros(n, n);
for j = 1:n
    V29(:, j) = f_norm'.^(j-1);
end

% Resolver sistema: V * c = Z
c29 = V29 \ Z';          % coeficientes

% Evaluar polinomio grado 29
Z_poly29      = polyval_vander(c29, f_fine_norm);
Z_poly29_eval = polyval_vander(c29, f_eval_norm);

fprintf('===========================================\n');
fprintf('  B1 — INTERPOLACIÓN MATRICIAL (grado 29)\n');
fprintf('===========================================\n');
fprintf('|Z|(1000 Hz) por polinomio grado 29 = %.4f Ω\n\n', Z_poly29_eval);

%% =========================================================
%% B1 — POLINOMIOS ESCALONADOS (grado 5, 10, 15)
%%      usando polyfit de MATLAB sobre datos normalizados
%% =========================================================
grados = [5, 10, 15];
colores = {[0.2 0.6 0.3], [0.8 0.4 0.0], [0.6 0.1 0.6]};
Z_escalonado = zeros(length(grados), length(f_fine));
Z_esc_eval   = zeros(1, length(grados));

for k = 1:length(grados)
    p = polyfit(f_norm, Z, grados(k));
    Z_escalonado(k,:) = polyval(p, f_fine_norm);
    Z_esc_eval(k)     = polyval(p, f_eval_norm);
    fprintf('Grado %2d → |Z|(1000 Hz) = %.4f Ω\n', grados(k), Z_esc_eval(k));
end
fprintf('\n');

%% =========================================================
%% B1 — INTERPOLACIÓN DE LAGRANGE (grado 29, todos los puntos)
%% =========================================================
Z_lagrange_eval = lagrange_interp(f, Z, f_eval);
Z_lagrange_fine = zeros(1, length(f_fine));
for i = 1:length(f_fine)
    Z_lagrange_fine(i) = lagrange_interp(f, Z, f_fine(i));
end

fprintf('|Z|(1000 Hz) por Lagrange (grado 29) = %.4f Ω\n\n', Z_lagrange_eval);

%% =========================================================
%% B1 — VALIDACIÓN LOO (Leave-One-Out) con 5 puntos al azar
%% =========================================================
rng(42);   % semilla para reproducibilidad
idx_loo = sort(randperm(n, 5));
errores_loo = zeros(1, 5);

fprintf('--- Validación LOO (Leave-One-Out) con 5 puntos ---\n');
fprintf('Puntos seleccionados: ');
fprintf('%d ', idx_loo); fprintf('\n\n');

for k = 1:5
    % Excluir punto k del conjunto de entrenamiento
    idx_train = setdiff(1:n, idx_loo(k));
    f_train = f(idx_train);
    Z_train = Z(idx_train);
    
    % Polinomio de grado 10 con datos reducidos
    f_tr_norm = (f_train - mean(f_train)) / std(f_train);
    f_te_norm = (f(idx_loo(k)) - mean(f_train)) / std(f_train);
    
    p_loo = polyfit(f_tr_norm, Z_train, 10);
    Z_pred = polyval(p_loo, f_te_norm);
    Z_real = Z(idx_loo(k));
    
    errores_loo(k) = abs(Z_pred - Z_real) / Z_real * 100;
    fprintf('Punto %2d: f=%5.0f Hz | Z_real=%.2f Ω | Z_pred=%.2f Ω | Error=%.4f%%\n', ...
            idx_loo(k), f(idx_loo(k)), Z_real, Z_pred, errores_loo(k));
end
fprintf('\nError relativo LOO promedio: %.4f%%\n\n', mean(errores_loo));

%% =========================================================
%% B2 — SPLINE CÚBICO NATURAL
%% =========================================================
cs = spline(f, Z);                     % spline cúbico de MATLAB
Z_spline_fine = ppval(cs, f_fine);    % evaluar en malla fina
Z_spline_eval = ppval(cs, f_eval);    % evaluar en 1000 Hz

fprintf('===========================================\n');
fprintf('  B2 — SPLINE CÚBICO NATURAL\n');
fprintf('===========================================\n');
fprintf('|Z|(1000 Hz) por spline cúbico = %.4f Ω\n\n', Z_spline_eval);

% Polinomio seleccionado para comparación (grado 10)
p10 = polyfit(f_norm, Z, 10);
Z_poly10_fine = polyval(p10, f_fine_norm);
Z_poly10_eval = polyval(p10, f_eval_norm);
fprintf('|Z|(1000 Hz) por polinomio grado 10 = %.4f Ω\n', Z_poly10_eval);

fprintf('\n--- Comparación en f = 1000 Hz ---\n');
fprintf('Spline cúbico : %.4f Ω\n', Z_spline_eval);
fprintf('Polinomio g10 : %.4f Ω\n', Z_poly10_eval);
fprintf('Diferencia    : %.4f Ω (%.4f%%)\n\n', ...
        abs(Z_spline_eval - Z_poly10_eval), ...
        abs(Z_spline_eval - Z_poly10_eval)/Z_spline_eval*100);

%% =========================================================
%% GRÁFICOS
%% =========================================================

%% Figura 1: Comparación Runge (grado 29 vs escalonados)
figure('Name','B1 — Efecto Runge','Position',[50 50 1100 650]);

subplot(2,1,1);
plot(f, Z, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5, 'DisplayName','Datos'); hold on;
plot(f_fine, Z_poly29, 'r-', 'LineWidth',1.5, 'DisplayName','Grado 29 (Matricial)');
ylim([125 175]);
grid on;
xlabel('f (Hz)'); ylabel('|Z| (Ω)');
title('Polinomio grado 29 — Oscilaciones de Runge visibles en extremos','FontSize',12);
legend('Location','north'); hold off;

subplot(2,1,2);
plot(f, Z, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5, 'DisplayName','Datos'); hold on;
nombres = {'Grado 5','Grado 10','Grado 15'};
for k = 1:3
    plot(f_fine, Z_escalonado(k,:), '-', 'Color', colores{k}, ...
         'LineWidth', 1.8, 'DisplayName', nombres{k});
end
ylim([128 165]);
grid on;
xlabel('f (Hz)'); ylabel('|Z| (Ω)');
title('Polinomios escalonados — Ajuste estable','FontSize',12);
legend('Location','northeast'); hold off;

sgtitle('B1 — Evidencia del Efecto Runge','FontWeight','bold','FontSize',14);

%% Figura 2: Lagrange vs Matricial (ambos grado 29)
figure('Name','B1 — Lagrange vs Matricial','Position',[100 80 900 500]);
plot(f, Z, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6, 'DisplayName','Datos'); hold on;
plot(f_fine, Z_lagrange_fine, 'b--', 'LineWidth',1.8, 'DisplayName','Lagrange (g29)');
plot(f_fine, Z_poly29, 'r-', 'LineWidth',1.5, 'DisplayName','Matricial (g29)');
plot(f_eval, Z_lagrange_eval, 'b^', 'MarkerFaceColor','b', 'MarkerSize',9, ...
     'DisplayName',sprintf('Lagrange: %.2f Ω',Z_lagrange_eval));
plot(f_eval, Z_poly29_eval, 'r^', 'MarkerFaceColor','r', 'MarkerSize',9, ...
     'DisplayName',sprintf('Matricial: %.2f Ω',Z_poly29_eval));
xline(f_eval, 'k--', 'f=1000 Hz', 'LabelVerticalAlignment','bottom','FontSize',9);
ylim([125 175]); grid on;
xlabel('f (Hz)'); ylabel('|Z| (Ω)');
title('B1 — Comparación: Lagrange vs Método Matricial (grado 29)','FontSize',12);
legend('Location','north','FontSize',9); hold off;

%% Figura 3: Spline vs Polinomio grado 10
figure('Name','B2 — Spline vs Polinomio','Position',[150 100 950 550]);
plot(f, Z, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7, 'DisplayName','Datos'); hold on;
plot(f_fine, Z_spline_fine, 'b-', 'LineWidth',2.2, 'DisplayName','Spline cúbico natural');
plot(f_fine, Z_poly10_fine, '--', 'Color',[0.8 0.4 0.0], 'LineWidth',1.8, ...
     'DisplayName','Polinomio grado 10');
plot(f_eval, Z_spline_eval, 'bs', 'MarkerFaceColor','b', 'MarkerSize',10, ...
     'DisplayName',sprintf('Spline: %.4f Ω', Z_spline_eval));
plot(f_eval, Z_poly10_eval, 's', 'Color',[0.8 0.4 0.0], 'MarkerFaceColor',[0.8 0.4 0.0], ...
     'MarkerSize',10, 'DisplayName',sprintf('Polinomio g10: %.4f Ω', Z_poly10_eval));
xline(f_eval, 'k--', 'f = 1000 Hz', 'LabelVerticalAlignment','bottom','FontSize',9);
grid on;
xlabel('f (Hz)'); ylabel('|Z| (Ω)');
title('B2 — Spline cúbico vs Polinomio grado 10','FontSize',12);
legend('Location','northeast','FontSize',9); hold off;

%% =========================================================
%% FUNCIONES AUXILIARES
%% =========================================================

function Z_out = polyval_vander(c, x_norm)
% Evalúa polinomio dado por coeficientes de Vandermonde (potencias crecientes)
    n = length(c);
    Z_out = zeros(size(x_norm));
    for j = 1:n
        Z_out = Z_out + c(j) .* x_norm.^(j-1);
    end
end

function Zq = lagrange_interp(x, y, xq)
% Interpolación de Lagrange clásica
    n  = length(x);
    Zq = 0;
    for i = 1:n
        Li = 1;
        for j = 1:n
            if j ~= i
                Li = Li * (xq - x(j)) / (x(i) - x(j));
            end
        end
        Zq = Zq + y(i) * Li;
    end
end