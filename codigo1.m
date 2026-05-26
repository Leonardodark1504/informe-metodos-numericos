%% =========================================================
%  PARTE 1: Interpolación de Lagrange (grado 2) y Spline Cúbico Natural
%  Sistema de telemetría biomédica — Biosensor resistivo-capacitivo
%  =========================================================
clc; clear; close all;

%% ── Datos del ensayo (50 mediciones) ────────────────────
f = [10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, ...
     35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, ...
     60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, ...
     85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5];

V = [0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551, ...
     1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041, ...
    -0.057,-0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809, ...
     0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004];

Z = [182.4, 178.9, 175.1, 171.0, 166.8, 162.7, 158.9, 155.4, 152.0, 149.0, ...
     146.1, 145.2, 145.8, 147.3, 149.9, 153.5, 158.0, 163.2, 168.9, 174.8, ...
     180.5, 186.2, 191.5, 196.2, 200.1, 203.1, 205.2, 206.3, 206.1, 204.7, ...
     198.0, 194.4, 190.9, 187.8, 185.1, 183.0, 181.6, 180.8, 180.6, 180.9];

%% ── Función de Interpolación de Lagrange grado 2 ─────────
% Recibe: xi (puntos), yi (valores), x (punto a interpolar)
lagrange2 = @(xi, yi, x) ...
    yi(1)*((x-xi(2))*(x-xi(3)))/((xi(1)-xi(2))*(xi(1)-xi(3))) + ...
    yi(2)*((x-xi(1))*(x-xi(3)))/((xi(2)-xi(1))*(xi(2)-xi(3))) + ...
    yi(3)*((x-xi(1))*(x-xi(2)))/((xi(3)-xi(1))*(xi(3)-xi(2)));

%% ── Puntos a estimar ─────────────────────────────────────
f1 = 41.0;   % kHz
f2 = 73.0;   % kHz

%% ══════════════════════════════════════════════════════════
%  MÉTODO 1: Lagrange de grado 2 — 3 puntos más cercanos
%  ══════════════════════════════════════════════════════════

% ── Para f1 = 41.0 kHz → vecinos: 37.5, 40.0, 42.5
idx1 = [12, 13, 14];   % índices en el vector (base 1)
xi1 = f(idx1);   yi1_V = V(idx1);   yi1_Z = Z(idx1);

V_lag_41  = lagrange2(xi1, yi1_V, f1);
Z_lag_41  = lagrange2(xi1, yi1_Z, f1);

% ── Para f2 = 73.0 kHz → vecinos: 70.0, 72.5, 75.0
idx2 = [25, 26, 27];
xi2 = f(idx2);   yi2_V = V(idx2);   yi2_Z = Z(idx2);

V_lag_73  = lagrange2(xi2, yi2_V, f2);
Z_lag_73  = lagrange2(xi2, yi2_Z, f2);

fprintf('\n══════ RESULTADOS LAGRANGE (grado 2) ══════\n');
fprintf('V(41.0 kHz)  = %.4f V\n',   V_lag_41);
fprintf('|Z|(41.0 kHz)= %.4f Ω\n',   Z_lag_41);
fprintf('V(73.0 kHz)  = %.4f V\n',   V_lag_73);
fprintf('|Z|(73.0 kHz)= %.4f Ω\n',   Z_lag_73);

%% ══════════════════════════════════════════════════════════
%  MÉTODO 2: Spline Cúbico Natural (toda la tabla)
%  ══════════════════════════════════════════════════════════
% MATLAB: spline() usa condición not-a-knot por defecto.
% Para spline natural usamos csape con condición 'variational'.

cs_V = csape(f, V, 'variational');   % Spline cúbico natural para V
cs_Z = csape(f, Z, 'variational');   % Spline cúbico natural para |Z|

V_spl_41  = fnval(cs_V, f1);
Z_spl_41  = fnval(cs_Z, f1);
V_spl_73  = fnval(cs_V, f2);
Z_spl_73  = fnval(cs_Z, f2);

fprintf('\n══════ RESULTADOS SPLINE CÚBICO NATURAL ══════\n');
fprintf('V(41.0 kHz)  = %.4f V\n',   V_spl_41);
fprintf('|Z|(41.0 kHz)= %.4f Ω\n',   Z_spl_41);
fprintf('V(73.0 kHz)  = %.4f V\n',   V_spl_73);
fprintf('|Z|(73.0 kHz)= %.4f Ω\n',   Z_spl_73);

%% ── Tabla comparativa ────────────────────────────────────
fprintf('\n══════ COMPARACIÓN DE MÉTODOS ══════════════════════════════\n');
fprintf('%-20s %-15s %-15s\n','Cantidad','Lagrange L2','Spline Cúbico');
fprintf('%-20s %-15.4f %-15.4f\n','V(41.0 kHz) [V]',  V_lag_41, V_spl_41);
fprintf('%-20s %-15.4f %-15.4f\n','|Z|(41.0 kHz) [Ω]',Z_lag_41, Z_spl_41);
fprintf('%-20s %-15.4f %-15.4f\n','V(73.0 kHz) [V]',  V_lag_73, V_spl_73);
fprintf('%-20s %-15.4f %-15.4f\n','|Z|(73.0 kHz) [Ω]',Z_lag_73, Z_spl_73);

%% ══════════════════════════════════════════════════════════
%  GRÁFICAS
%  ══════════════════════════════════════════════════════════
f_fine = linspace(min(f), max(f), 1000);   % eje frecuencia fino

V_spline_fine = fnval(cs_V, f_fine);
Z_spline_fine = fnval(cs_Z, f_fine);

% Colores corporativos
c_data   = [0.15 0.35 0.70];   % azul datos
c_lag    = [0.85 0.33 0.10];   % naranja Lagrange
c_spl    = [0.10 0.65 0.40];   % verde Spline
c_mark41 = [0.90 0.10 0.20];   % rojo punto 41 kHz
c_mark73 = [0.60 0.10 0.80];   % violeta punto 73 kHz

%% ── Figura 1: V(f) ───────────────────────────────────────
figure('Name','Parte 1 — Voltaje V(f)','Color','w','Position',[50 50 900 500]);

% Spline sobre todo el rango
plot(f_fine, V_spline_fine, '-', 'Color', c_spl, 'LineWidth', 2.2); hold on;

% Datos originales
scatter(f, V, 40, c_data, 'filled', 'MarkerEdgeColor','w', 'LineWidth',0.8);

% Puntos Lagrange (marcadores especiales)
scatter(f1, V_lag_41, 120, c_lag, 'd', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
scatter(f2, V_lag_73, 120, c_lag, 'd', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);

% Puntos Spline
scatter(f1, V_spl_41, 120, c_spl, 's', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
scatter(f2, V_spl_73, 120, c_spl, 's', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);

% Líneas verticales de referencia
xline(f1, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'Alpha', 0.7);
xline(f2, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'Alpha', 0.7);
yline(0,  '-',  'Color', [0.3 0.3 0.3], 'LineWidth', 0.8);

% Anotaciones
text(f1+0.5, V_lag_41+0.05, sprintf('Lag: %.4f V', V_lag_41), ...
    'Color', c_lag, 'FontSize', 8, 'FontWeight','bold');
text(f1+0.5, V_spl_41-0.07, sprintf('Spl: %.4f V', V_spl_41), ...
    'Color', c_spl, 'FontSize', 8, 'FontWeight','bold');
text(f2+0.5, V_lag_73+0.05, sprintf('Lag: %.4f V', V_lag_73), ...
    'Color', c_lag, 'FontSize', 8, 'FontWeight','bold');
text(f2+0.5, V_spl_73-0.07, sprintf('Spl: %.4f V', V_spl_73), ...
    'Color', c_spl, 'FontSize', 8, 'FontWeight','bold');

xlabel('Frecuencia f (kHz)', 'FontSize', 12, 'FontWeight','bold');
ylabel('Voltaje V (V)',       'FontSize', 12, 'FontWeight','bold');
title('Parte 1 — Voltaje de salida V(f): Lagrange vs Spline Cúbico Natural', ...
      'FontSize', 13, 'FontWeight','bold');
legend('Spline cúbico natural','Datos medidos', ...
       'Lagrange f=41 kHz','Lagrange f=73 kHz', ...
       'Spline f=41 kHz',  'Spline f=73 kHz', ...
       'Location','southwest','FontSize',9);
grid on; grid minor;
xlim([min(f)-1, max(f)+1]);
set(gca,'FontSize',10,'Box','on');
hold off;

%% ── Figura 2: |Z|(f) ─────────────────────────────────────
figure('Name','Parte 1 — Impedancia |Z|(f)','Color','w','Position',[60 60 900 500]);

plot(f_fine, Z_spline_fine, '-', 'Color', c_spl, 'LineWidth', 2.2); hold on;
scatter(f, Z, 40, c_data, 'filled', 'MarkerEdgeColor','w', 'LineWidth',0.8);
scatter(f1, Z_lag_41, 120, c_lag, 'd', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
scatter(f2, Z_lag_73, 120, c_lag, 'd', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
scatter(f1, Z_spl_41, 120, c_spl, 's', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
scatter(f2, Z_spl_73, 120, c_spl, 's', 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
xline(f1, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'Alpha', 0.7);
xline(f2, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'Alpha', 0.7);

text(f1+0.5, Z_lag_41+1.5, sprintf('Lag: %.2f Ω', Z_lag_41), ...
    'Color', c_lag, 'FontSize', 8, 'FontWeight','bold');
text(f1+0.5, Z_spl_41-2.5, sprintf('Spl: %.2f Ω', Z_spl_41), ...
    'Color', c_spl, 'FontSize', 8, 'FontWeight','bold');
text(f2+0.5, Z_lag_73+1.5, sprintf('Lag: %.2f Ω', Z_lag_73), ...
    'Color', c_lag, 'FontSize', 8, 'FontWeight','bold');
text(f2+0.5, Z_spl_73-2.5, sprintf('Spl: %.2f Ω', Z_spl_73), ...
    'Color', c_spl, 'FontSize', 8, 'FontWeight','bold');

xlabel('Frecuencia f (kHz)',       'FontSize', 12, 'FontWeight','bold');
ylabel('Impedancia |Z| (Ω)',       'FontSize', 12, 'FontWeight','bold');
title('Parte 1 — Impedancia |Z(f)|: Lagrange vs Spline Cúbico Natural', ...
      'FontSize', 13, 'FontWeight','bold');
legend('Spline cúbico natural','Datos medidos', ...
       'Lagrange f=41 kHz','Lagrange f=73 kHz', ...
       'Spline f=41 kHz',  'Spline f=73 kHz', ...
       'Location','southeast','FontSize',9);
grid on; grid minor;
xlim([min(f)-1, max(f)+1]);
set(gca,'FontSize',10,'Box','on');
hold off;

%% ── Figura 3: Ampliación zona 41 kHz ────────────────────
figure('Name','Zoom 41 kHz','Color','w','Position',[70 70 700 400]);

idx_zoom1 = (f >= 35) & (f <= 50);
f_zoom1   = linspace(35, 50, 300);
V_zoom1   = fnval(cs_V, f_zoom1);

plot(f_zoom1, V_zoom1, '-', 'Color', c_spl, 'LineWidth', 2); hold on;
scatter(f(idx_zoom1), V(idx_zoom1), 60, c_data, 'filled','MarkerEdgeColor','w');
scatter(f1, V_lag_41, 150, c_lag, 'd', 'filled','MarkerEdgeColor','k','LineWidth',1.5);
scatter(f1, V_spl_41, 150, c_spl, 's', 'filled','MarkerEdgeColor','k','LineWidth',1.5);
xline(f1,'--','Color',[0.5 0.5 0.5],'LineWidth',1.2);
xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('V (V)',  'FontSize',11,'FontWeight','bold');
title(sprintf('Zoom zona 41 kHz  |  Lag=%.4f V  Spl=%.4f V', V_lag_41, V_spl_41), ...
      'FontSize',11,'FontWeight','bold');
legend('Spline','Datos','Lagrange','Spline','Location','northeast','FontSize',9);
grid on; grid minor; hold off;

%% ── Figura 4: Ampliación zona 73 kHz ────────────────────
figure('Name','Zoom 73 kHz','Color','w','Position',[80 80 700 400]);

idx_zoom2 = (f >= 65) & (f <= 80);
f_zoom2   = linspace(65, 80, 300);
V_zoom2   = fnval(cs_V, f_zoom2);

plot(f_zoom2, V_zoom2, '-', 'Color', c_spl, 'LineWidth', 2); hold on;
scatter(f(idx_zoom2), V(idx_zoom2), 60, c_data, 'filled','MarkerEdgeColor','w');
scatter(f2, V_lag_73, 150, c_lag, 'd', 'filled','MarkerEdgeColor','k','LineWidth',1.5);
scatter(f2, V_spl_73, 150, c_spl, 's', 'filled','MarkerEdgeColor','k','LineWidth',1.5);
xline(f2,'--','Color',[0.5 0.5 0.5],'LineWidth',1.2);
xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('V (V)',  'FontSize',11,'FontWeight','bold');
title(sprintf('Zoom zona 73 kHz  |  Lag=%.4f V  Spl=%.4f V', V_lag_73, V_spl_73), ...
      'FontSize',11,'FontWeight','bold');
legend('Spline','Datos','Lagrange','Spline','Location','northwest','FontSize',9);
grid on; grid minor; hold off;

fprintf('\n✔ Gráficas generadas correctamente.\n');

%% ══════════════════════════════════════════════════════════
%  DISCUSIÓN: Preguntas 3 y 4
%  ══════════════════════════════════════════════════════════
fprintf('\n══════ DISCUSIÓN TÉCNICA ══════════════════════════════════\n');
fprintf(['Pregunta 3 — Comparación de métodos:\n' ...
'  Ambos métodos producen resultados muy cercanos para puntos interiores\n' ...
'  bien respaldados por datos vecinos. El spline cúbico es más confiable\n' ...
'  porque minimiza la curvatura global (condición natural S''''=0 en extremos)\n' ...
'  y garantiza continuidad C2 en cada nodo, evitando oscilaciones locales.\n' ...
'  Lagrange de grado 2 puede generar errores si la función no es localmente\n' ...
'  cuadrática, especialmente cerca de inflexiones como la región 55-65 kHz.\n\n']);

fprintf(['Pregunta 4 — Spline vs polinomio global de alto grado:\n' ...
'  El polinomio global de grado n-1 (con n=40 puntos → grado 39) sufre\n' ...
'  el fenómeno de Runge: oscilaciones incontroladas en los extremos del\n' ...
'  intervalo. El spline cúbico usa polinomios de grado 3 por tramo,\n' ...
'  controlando la suavidad localmente sin propagación de error global.\n' ...
'  Esto es crítico en señales biomédicas donde la estabilidad numérica\n' ...
'  y la interpretación física son igualmente importantes.\n']);