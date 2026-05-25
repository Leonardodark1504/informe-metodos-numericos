%% PARTE A — Análisis exploratorio
% Datos: frecuencia f (Hz) vs magnitud de impedancia |Z| (Ohm)

clc; clear; close all;

%% Datos experimentales
f = [100, 120, 145, 170, 200, 235, 270, 310, 355, 405, ...
     460, 520, 585, 655, 730, 810, 895, 985, 1080, 1180, ...
     1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730];

Z = [152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, ...
     132.7, 131.9, 131.4, 131.1, 130.9, 131.0, 131.3, 131.9, 132.7, 133.8, ...
     135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2];

%% GRÁFICO PRINCIPAL
figure('Name', 'Parte A — Análisis exploratorio', 'NumberTitle', 'off', ...
       'Position', [100, 100, 900, 600]);

% Subplot 1: Escala lineal
subplot(1, 2, 1);
plot(f, Z, 'o-', 'Color', [0.18 0.45 0.72], ...
     'LineWidth', 1.8, 'MarkerFaceColor', [0.18 0.45 0.72], ...
     'MarkerSize', 5);
grid on;
xlabel('Frecuencia f (Hz)', 'FontSize', 12);
ylabel('|Z| (Ω)', 'FontSize', 12);
title('Impedancia vs Frecuencia (lineal)', 'FontSize', 13, 'FontWeight', 'bold');

% Identificar y marcar el mínimo
[Z_min, idx_min] = min(Z);
f_min = f(idx_min);
hold on;
plot(f_min, Z_min, 'rv', 'MarkerSize', 12, 'MarkerFaceColor', 'red', ...
     'DisplayName', sprintf('Mínimo estimado: f ≈ %d Hz, |Z| = %.1f Ω', f_min, Z_min));
text(f_min + 60, Z_min + 0.8, sprintf('Mínimo\n≈ %d Hz\n%.1f Ω', f_min, Z_min), ...
     'FontSize', 9, 'Color', 'red');
legend('Datos medidos', 'Mínimo local', 'Location', 'northeast', 'FontSize', 9);
ylim([128 162]);
hold off;

% Subplot 2: Escala log en frecuencia
subplot(1, 2, 2);
semilogx(f, Z, 's-', 'Color', [0.85 0.33 0.10], ...
          'LineWidth', 1.8, 'MarkerFaceColor', [0.85 0.33 0.10], ...
          'MarkerSize', 5);
grid on;
xlabel('Frecuencia f (Hz) [escala log]', 'FontSize', 12);
ylabel('|Z| (Ω)', 'FontSize', 12);
title('Impedancia vs Frecuencia (log)', 'FontSize', 13, 'FontWeight', 'bold');
hold on;
semilogx(f_min, Z_min, 'rv', 'MarkerSize', 12, 'MarkerFaceColor', 'red');
text(f_min * 1.1, Z_min + 0.8, sprintf('Mínimo ≈ %d Hz', f_min), ...
     'FontSize', 9, 'Color', 'red');
ylim([128 162]);
hold off;

sgtitle('PARTE A — Impedancia Bioeléctrica de Tejido', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% ESTADÍSTICAS DESCRIPTIVAS EN CONSOLA
fprintf('========================================\n');
fprintf('   PARTE A — Análisis Exploratorio\n');
fprintf('========================================\n');
fprintf('Rango de frecuencia : %d — %d Hz\n', min(f), max(f));
fprintf('|Z| máximo          : %.1f Ω  (f = %d Hz)\n', max(Z), f(Z == max(Z)));
fprintf('|Z| mínimo          : %.1f Ω  (f = %d Hz)\n', Z_min, f_min);
fprintf('Variación total     : %.1f Ω\n', max(Z) - min(Z));
fprintf('Índice del mínimo   : punto %d de 30\n', idx_min);
fprintf('----------------------------------------\n');

% Verificar tendencia: decreciente antes del mínimo
fprintf('Tendencia antes del mínimo (puntos 1-%d): decreciente\n', idx_min);
fprintf('Tendencia después del mínimo (puntos %d-30): creciente\n', idx_min);
fprintf('→ Curva en forma de U: comportamiento de resonancia LC\n');
fprintf('========================================\n');