%% =========================================================
%  PARTE 2: Derivación Numérica de V(f)
%  Diferencias centradas O(h²) y O(h⁴), progresiva O(h²) en extremo
%  Sistema de telemetría biomédica
%  =========================================================
clc; clear; close all;

%% ── Datos del ensayo ─────────────────────────────────────
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

h = 2.5;   % paso uniforme [kHz]
n = length(f);

%% ══════════════════════════════════════════════════════════
%  FUNCIONES DE DIFERENCIAS FINITAS
%  ══════════════════════════════════════════════════════════

% Diferencia centrada orden 2: O(h²)
% dV/df ≈ (V_{i+1} - V_{i-1}) / (2h)
cent_O2 = @(V, i) (V(i+1) - V(i-1)) / (2*h);

% Diferencia centrada orden 4: O(h⁴)
% dV/df ≈ (-V_{i+2} + 8V_{i+1} - 8V_{i-1} + V_{i-2}) / (12h)
cent_O4 = @(V, i) (-V(i+2) + 8*V(i+1) - 8*V(i-1) + V(i-2)) / (12*h);

% Diferencia progresiva orden 2: O(h²)  [para extremo inferior]
% dV/df ≈ (-3V_0 + 4V_1 - V_2) / (2h)
prog_O2 = @(V, i) (-3*V(i) + 4*V(i+1) - V(i+2)) / (2*h);

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 1: dV/df en f = 40.0, 70.0 y 100.0 kHz
%  ══════════════════════════════════════════════════════════

% Índices (base 1 en MATLAB)
i_40  = find(f == 40.0);    % índice 13
i_70  = find(f == 70.0);    % índice 25
i_100 = find(f == 100.0);   % índice 37

%% ── f = 40.0 kHz ─────────────────────────────────────────
dV_40_O2 = cent_O2(V, i_40);
dV_40_O4 = cent_O4(V, i_40);

fprintf('\n══════ PREGUNTA 1: dV/df en puntos interiores ══════════\n\n');
fprintf('► f = 40.0 kHz\n');
fprintf('  Puntos usados O2: f=%.1f(%.3f), f=%.1f(%.3f)\n', ...
    f(i_40-1), V(i_40-1), f(i_40+1), V(i_40+1));
fprintf('  Puntos usados O4: f=%.1f(%.3f), f=%.1f(%.3f), f=%.1f(%.3f), f=%.1f(%.3f)\n', ...
    f(i_40-2),V(i_40-2), f(i_40-1),V(i_40-1), f(i_40+1),V(i_40+1), f(i_40+2),V(i_40+2));
fprintf('  dV/df [O2] = (%.3f - %.3f) / (2×%.1f) = %.6f V/kHz\n', ...
    V(i_40+1), V(i_40-1), h, dV_40_O2);
fprintf('  dV/df [O4] = (-%.3f + 8×%.3f - 8×%.3f + %.3f) / (12×%.1f) = %.6f V/kHz\n', ...
    V(i_40+2), V(i_40+1), V(i_40-1), V(i_40-2), h, dV_40_O4);

%% ── f = 70.0 kHz ─────────────────────────────────────────
dV_70_O2 = cent_O2(V, i_70);
dV_70_O4 = cent_O4(V, i_70);

fprintf('\n► f = 70.0 kHz\n');
fprintf('  dV/df [O2] = (%.3f - %.3f) / (2×%.1f) = %.6f V/kHz\n', ...
    V(i_70+1), V(i_70-1), h, dV_70_O2);
fprintf('  dV/df [O4] = (-%.3f + 8×%.3f - 8×%.3f + %.3f) / (12×%.1f) = %.6f V/kHz\n', ...
    V(i_70+2), V(i_70+1), V(i_70-1), V(i_70-2), h, dV_70_O4);

%% ── f = 100.0 kHz ────────────────────────────────────────
dV_100_O2 = cent_O2(V, i_100);
dV_100_O4 = cent_O4(V, i_100);

fprintf('\n► f = 100.0 kHz\n');
fprintf('  dV/df [O2] = (%.3f - %.3f) / (2×%.1f) = %.6f V/kHz\n', ...
    V(i_100+1), V(i_100-1), h, dV_100_O2);
fprintf('  dV/df [O4] = (-%.3f + 8×%.3f - 8×%.3f + %.3f) / (12×%.1f) = %.6f V/kHz\n', ...
    V(i_100+2), V(i_100+1), V(i_100-1), V(i_100-2), h, dV_100_O4);

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 2: dV/df en f = 10.0 kHz — Progresiva O(h²)
%  ══════════════════════════════════════════════════════════
i_10 = find(f == 10.0);   % índice 1

dV_10_prog = prog_O2(V, i_10);

fprintf('\n══════ PREGUNTA 2: dV/df en extremo f = 10.0 kHz ══════\n');
fprintf('  Fórmula: (-3V₀ + 4V₁ - V₂) / (2h)\n');
fprintf('  = (-3×%.3f + 4×%.3f - %.3f) / (2×%.1f)\n', ...
    V(i_10), V(i_10+1), V(i_10+2), h);
fprintf('  = (%.4f) / %.1f = %.6f V/kHz\n', ...
    -3*V(i_10)+4*V(i_10+1)-V(i_10+2), 2*h, dV_10_prog);

%% ── Tabla resumen ────────────────────────────────────────
fprintf('\n══════ TABLA RESUMEN dV/df ══════════════════════════════\n');
fprintf('%-12s %-18s %-18s %-15s\n','f (kHz)','O2 (V/kHz)','O4 (V/kHz)','Método');
fprintf('%-12.1f %-18.6f %-18s %-15s\n', 10.0, dV_10_prog,'    —','Progresiva');
fprintf('%-12.1f %-18.6f %-18.6f %-15s\n', 40.0, dV_40_O2, dV_40_O4,'Centrada');
fprintf('%-12.1f %-18.6f %-18.6f %-15s\n', 70.0, dV_70_O2, dV_70_O4,'Centrada');
fprintf('%-12.1f %-18.6f %-18.6f %-15s\n',100.0,dV_100_O2,dV_100_O4,'Centrada');

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 4: Derivada del spline cúbico (fnder)
%  ══════════════════════════════════════════════════════════
cs_V = csape(f, V, 'variational');   % spline cúbico natural
ds_V = fnder(cs_V, 1);               % derivada analítica del spline

dV_spl_10  = fnval(ds_V, 10.0);
dV_spl_40  = fnval(ds_V, 40.0);
dV_spl_70  = fnval(ds_V, 70.0);
dV_spl_100 = fnval(ds_V,100.0);

fprintf('\n══════ PREGUNTA 4: Derivada del Spline Cúbico ══════════\n');
fprintf('%-12s %-18s %-18s %-18s\n','f (kHz)','Spline (V/kHz)','O2 (V/kHz)','Dif. Spl-O2');
fprintf('%-12.1f %-18.6f %-18.6f %-18.6f\n', 10.0,dV_spl_10, dV_10_prog,  dV_spl_10 -dV_10_prog);
fprintf('%-12.1f %-18.6f %-18.6f %-18.6f\n', 40.0,dV_spl_40, dV_40_O2,   dV_spl_40 -dV_40_O2);
fprintf('%-12.1f %-18.6f %-18.6f %-18.6f\n', 70.0,dV_spl_70, dV_70_O2,   dV_spl_70 -dV_70_O2);
fprintf('%-12.1f %-18.6f %-18.6f %-18.6f\n',100.0,dV_spl_100,dV_100_O2,  dV_spl_100-dV_100_O2);

%% ══════════════════════════════════════════════════════════
%  DERIVADA NUMÉRICA PARA TODO EL DOMINIO (visualización)
%  ══════════════════════════════════════════════════════════

% Diferencias centradas O2 para todos los puntos interiores
dV_num = zeros(1, n);
dV_num(1)   = prog_O2(V, 1);                                    % extremo inf
dV_num(n)   = (3*V(n) - 4*V(n-1) + V(n-2)) / (2*h);           % extremo sup (regresiva O2)
for i = 2:n-1
    dV_num(i) = cent_O2(V, i);
end

% Derivada del spline en todo el dominio
f_fine     = linspace(min(f), max(f), 1000);
dV_spl_all = fnval(ds_V, f_fine);

%% ══════════════════════════════════════════════════════════
%  GRÁFICAS
%  ══════════════════════════════════════════════════════════

% Paleta de colores
c_data = [0.15 0.35 0.70];
c_O2   = [0.85 0.33 0.10];
c_O4   = [0.50 0.10 0.70];
c_spl  = [0.10 0.65 0.40];
c_zero = [0.40 0.40 0.40];

%% ── Figura 1: V(f) con derivadas marcadas ────────────────
figure('Name','V(f) y puntos de derivación','Color','w','Position',[50 50 950 420]);

cs_V_full = csape(f, V, 'variational');
V_fine    = fnval(cs_V_full, f_fine);

plot(f_fine, V_fine, '-', 'Color', c_spl, 'LineWidth', 2); hold on;
scatter(f, V, 35, c_data, 'filled', 'MarkerEdgeColor','w');
yline(0, '-', 'Color', c_zero, 'LineWidth', 1);

% Marcar puntos donde se calcula la derivada
pts_f = [10, 40, 70, 100];
pts_V = arrayfun(@(fi) fnval(cs_V_full, fi), pts_f);
scatter(pts_f, pts_V, 120, [0.9 0.1 0.2], 'p', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.2);

for k = 1:length(pts_f)
    text(pts_f(k)+0.5, pts_V(k)+0.06, sprintf('f=%g', pts_f(k)), ...
        'FontSize', 8, 'FontWeight','bold', 'Color',[0.7 0.1 0.1]);
end

xlabel('Frecuencia f (kHz)', 'FontSize',12, 'FontWeight','bold');
ylabel('V (V)',               'FontSize',12, 'FontWeight','bold');
title('V(f) — Puntos donde se estima dV/df', 'FontSize',13, 'FontWeight','bold');
legend('Spline cúbico','Datos medidos','Puntos de derivación', ...
       'Location','southwest', 'FontSize',9);
grid on; grid minor;
xlim([8 110]); set(gca,'FontSize',10,'Box','on'); hold off;

%% ── Figura 2: dV/df — Diferencias numéricas vs Spline ───
figure('Name','dV/df numérica vs spline','Color','w','Position',[60 60 950 480]);

% Spline derivada (continua)
plot(f_fine, dV_spl_all, '-', 'Color', c_spl, 'LineWidth', 2.2); hold on;

% Diferencias centradas O2 (puntos interiores)
scatter(f, dV_num, 45, c_O2, 'o', 'filled', 'MarkerEdgeColor','w');

% Diferencias centradas O4 (solo puntos donde se tiene i±2)
f_O4   = [40.0, 70.0, 100.0];
dV_O4  = [dV_40_O4, dV_70_O4, dV_100_O4];
scatter(f_O4, dV_O4, 80, c_O4, 'd', 'filled', 'MarkerEdgeColor','k','LineWidth',1);

% Extremo inferior progresiva
scatter(10, dV_10_prog, 100, [0.9 0.6 0.1], 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.2);

% Línea cero
yline(0, '--', 'Color', c_zero, 'LineWidth', 1, 'Alpha', 0.8);

% Anotaciones en los 4 puntos clave
puntos_f  = [10.0,    40.0,    70.0,    100.0];
puntos_O2 = [dV_10_prog, dV_40_O2, dV_70_O2, dV_100_O2];
puntos_O4 = [NaN,        dV_40_O4, dV_70_O4, dV_100_O4];
offset_y  = [0.003, -0.006,  0.003,   0.002];

for k = 1:4
    texto = sprintf('f=%g\nO2=%.4f', puntos_f(k), puntos_O2(k));
    text(puntos_f(k)+1, puntos_O2(k)+offset_y(k), texto, ...
        'FontSize',7.5, 'Color', c_O2, 'FontWeight','bold');
end

xlabel('Frecuencia f (kHz)', 'FontSize',12, 'FontWeight','bold');
ylabel('dV/df (V/kHz)',       'FontSize',12, 'FontWeight','bold');
title('Derivada dV/df: Diferencias Numéricas vs Spline Cúbico', ...
      'FontSize',13, 'FontWeight','bold');
legend('Spline d/df (analítica del spline)', ...
       'Centrada O(h²) — todos los puntos', ...
       'Centrada O(h⁴) — f=40, 70, 100 kHz', ...
       'Progresiva O(h²) — f=10 kHz (extremo)', ...
       'Location','northeast', 'FontSize',9);
grid on; grid minor;
xlim([8 110]); set(gca,'FontSize',10,'Box','on'); hold off;

%% ── Figura 3: Zoom en los 4 puntos clave ────────────────
figure('Name','Zoom derivadas','Color','w','Position',[70 70 950 440]);

f_pts = [10, 40, 70, 100];
col   = {[0.9 0.6 0.1], c_O2, c_O2, c_O2};   % colores por método
lab   = {'Prog O2','Cent O2','Cent O2','Cent O2'};
dv_O2 = [dV_10_prog, dV_40_O2, dV_70_O2, dV_100_O2];
dv_O4 = [NaN,        dV_40_O4, dV_70_O4, dV_100_O4];
dv_sp = [dV_spl_10,  dV_spl_40, dV_spl_70, dV_spl_100];
labels = {'f = 10 kHz','f = 40 kHz','f = 70 kHz','f = 100 kHz'};

x_cat = categorical(labels, labels);
bar_data = [dv_O2; dv_O4; dv_sp]';   % filas=puntos, cols=métodos

b = bar(x_cat, bar_data, 'grouped');
b(1).FaceColor = c_O2;
b(2).FaceColor = c_O4;
b(3).FaceColor = c_spl;

% Etiquetas de valor sobre cada barra
for j = 1:3
    for k = 1:4
        val = bar_data(k, j);
        if ~isnan(val)
            text(b(j).XEndPoints(k), val + 0.0005*sign(val), ...
                sprintf('%.4f', val), ...
                'HorizontalAlignment','center', 'FontSize', 7.5, ...
                'FontWeight','bold', 'Color', [0.2 0.2 0.2]);
        end
    end
end

yline(0, '-k', 'LineWidth', 1);
ylabel('dV/df (V/kHz)', 'FontSize',12, 'FontWeight','bold');
title('Comparación de métodos de derivación en los 4 puntos', ...
      'FontSize',13, 'FontWeight','bold');
legend('Centrada / Prog. O(h²)','Centrada O(h⁴)','Spline cúbico', ...
       'Location','northeast','FontSize',9);
grid on; grid minor;
set(gca,'FontSize',10,'Box','on');

%% ── Figura 4: Error relativo entre métodos ───────────────
figure('Name','Error relativo O2 vs Spline','Color','w','Position',[80 80 750 380]);

err_rel = abs((dV_num - fnval(ds_V, f)) ./ (fnval(ds_V, f) + 1e-12)) * 100;

bar(f, err_rel, 'FaceColor', c_O2, 'EdgeColor','none', 'FaceAlpha',0.75);
xlabel('Frecuencia f (kHz)', 'FontSize',12, 'FontWeight','bold');
ylabel('Error relativo (%)',  'FontSize',12, 'FontWeight','bold');
title('Error relativo: Diferencia centrada O(h²) vs Spline cúbico (referencia)', ...
      'FontSize',12, 'FontWeight','bold');
grid on; grid minor;
set(gca,'FontSize',10,'Box','on');
xlim([8 110]);

fprintf('\n✔ Gráficas Parte 2 generadas correctamente.\n');

%% ══════════════════════════════════════════════════════════
%  DISCUSIÓN FINAL — Preguntas 3 y 4
%  ══════════════════════════════════════════════════════════
fprintf('\n══════ PREGUNTA 3: Interpretación física ════════════════\n');
fprintf([
'  f=10 kHz (+0.0264 V/kHz): V crece con f → sensor en rampa ascendente.\n'...
'  El front-end gana amplitud al aumentar frecuencia. Zona de alta\n'...
'  sensibilidad positiva, no óptima para operar el sistema.\n\n'...
'  f=40 kHz (-0.0718 V/kHz): V cae rápidamente → zona post-pico (pico ~32.5 kHz).\n'...
'  La señal es muy sensible a cambios de frecuencia (mayor |dV/df|).\n'...
'  Operar aquí implica INESTABILIDAD ante variaciones de excitación.\n\n'...
'  f=70 kHz (+0.0444 V/kHz): V crece de nuevo tras el cruce por cero.\n'...
'  El sensor se recupera en esta zona. Derivada positiva moderada.\n\n'...
'  f=100 kHz (+0.0120 V/kHz): V crece muy lentamente → zona casi plana.\n'...
'  Menor sensibilidad a cambios de f = MÁXIMA ESTABILIDAD del sistema.\n'...
'  Zona recomendada para operación del biosensor.\n'
]);

fprintf('\n══════ PREGUNTA 4: Derivada spline vs diferencias finitas ══\n');
fprintf([
'  Los valores del spline y las diferencias centradas O2/O4 son muy\n'...
'  similares en zonas suaves (f=70, f=100), donde la función varía\n'...
'  de manera regular y los h=2.5 kHz son adecuados.\n\n'...
'  Diferencias notables ocurren cerca de:\n'...
'    - f=32.5 kHz (pico local): la función tiene curvatura alta y las\n'...
'      diferencias finitas no capturan bien la tangente exacta.\n'...
'    - f=57.5–62.5 kHz (cruce por cero): cambio de signo brusco.\n\n'...
'  El spline deriva analíticamente cada polinomio cúbico del tramo,\n'...
'  lo que da una estimación más suave y continua de la derivada.\n'...
'  Las diferencias finitas son simples y robustas, pero discretas.\n'...
'  RECOMENDACIÓN: usar el spline como referencia y las diferencias\n'...
'  finitas como verificación rápida.\n'
]);