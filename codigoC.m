%% PARTE C — DERIVACIÓN NUMÉRICA CON SPLINE CÚBICO
% C1: Primera derivada d|Z|/df → localizar mínimo con precisión
% C2: Segunda derivada d²|Z|/df² en el mínimo → estabilidad

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

%% =========================================================
%% CONSTRUCCIÓN DEL SPLINE CÚBICO
%% =========================================================
cs = spline(f, Z);          % spline cúbico natural (not-a-knot por defecto en MATLAB)
f_fine = linspace(min(f), max(f), 10000);   % malla fina para derivación

%% =========================================================
%% C1 — PRIMERA DERIVADA ANALÍTICA DEL SPLINE
%% Derivar el spline tramo a tramo: si S(f) = a+b(f-fi)+c(f-fi)²+d(f-fi)³
%% entonces S'(f) = b + 2c(f-fi) + 3d(f-fi)²
%% =========================================================
cs_d1 = fnder(cs, 1);       % derivada primera del spline (analítica)
cs_d2 = fnder(cs, 2);       % derivada segunda del spline (analítica)

% Evaluar derivadas en malla fina
dZ_df_fine  = ppval(cs_d1, f_fine);   % d|Z|/df en malla fina
d2Z_df2_fine = ppval(cs_d2, f_fine);  % d²|Z|/df² en malla fina

% Evaluar derivadas en los 30 puntos de datos originales
dZ_df_nodos  = ppval(cs_d1, f);
d2Z_df2_nodos = ppval(cs_d2, f);

%% =========================================================
%% LOCALIZACIÓN PRECISA DEL MÍNIMO
%% El mínimo ocurre donde d|Z|/df = 0 y cambia de negativo a positivo
%% =========================================================
% Buscar cambio de signo en la malla fina
signo = sign(dZ_df_fine);
idx_cambio = find(diff(signo) > 0, 1);  % primer cambio - a +

% Refinar con bisección sobre el spline derivado
f_a = f_fine(idx_cambio);
f_b = f_fine(idx_cambio + 1);
tol = 1e-6;
while (f_b - f_a) > tol
    f_mid = (f_a + f_b) / 2;
    if ppval(cs_d1, f_mid) < 0
        f_a = f_mid;
    else
        f_b = f_mid;
    end
end
f_minimo = (f_a + f_b) / 2;
Z_minimo = ppval(cs, f_minimo);

% Segunda derivada en el mínimo
d2Z_minimo = ppval(cs_d2, f_minimo);

fprintf('===========================================\n');
fprintf('  PARTE C — DERIVACIÓN CON SPLINE CÚBICO\n');
fprintf('===========================================\n');
fprintf('\n--- C1: Primera derivada ---\n');
fprintf('Frecuencia del mínimo (d|Z|/df = 0): f* = %.4f Hz\n', f_minimo);
fprintf('Impedancia en el mínimo : |Z|(f*) = %.4f Ω\n', Z_minimo);
fprintf('\n--- C2: Segunda derivada en el mínimo ---\n');
fprintf('d²|Z|/df² en f* = %.6f Ω/Hz²\n', d2Z_minimo);
if d2Z_minimo > 0
    fprintf('Signo POSITIVO → mínimo ESTABLE (concavidad hacia arriba)\n');
else
    fprintf('Signo NEGATIVO → sería un máximo (revisar)\n');
end

% Mostrar derivadas en todos los nodos
fprintf('\n--- Primera derivada en los 30 nodos ---\n');
fprintf('%-6s %-10s %-12s %-12s\n','Punto','f (Hz)','|Z| (Ω)','d|Z|/df');
fprintf('%s\n', repmat('-',1,44));
for i = 1:length(f)
    fprintf('%-6d %-10.0f %-12.2f %-12.6f\n', i, f(i), Z(i), dZ_df_nodos(i));
end

%% =========================================================
%% DEPENDENCIA DEL ERROR CON EL ESPACIAMIENTO
%% Para diferencias finitas (referencia teórica):
%% Error ~ O(h²) para diferencias centradas
%% Para el spline: el error de derivación es O(h^4) en el interior
%% =========================================================
fprintf('\n--- Análisis del error de derivación ---\n');
h_promedio = mean(diff(f));
h_max = max(diff(f));
h_min = min(diff(f));
fprintf('Espaciamiento promedio : %.1f Hz\n', h_promedio);
fprintf('Espaciamiento máximo   : %.1f Hz\n', h_max);
fprintf('Espaciamiento mínimo   : %.1f Hz\n', h_min);
fprintf('Zona del mínimo (±1 nodo): Δf = %.0f Hz\n', ...
        f(16) - f(14));

%% =========================================================
%% GRÁFICOS
%% =========================================================

%% Figura 1: |Z|(f) con mínimo marcado precisamente
figure('Name','C — Impedancia y mínimo preciso','Position',[50 50 950 420]);
Z_fine = ppval(cs, f_fine);
plot(f, Z, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5, 'DisplayName','Datos'); hold on;
plot(f_fine, Z_fine, 'b-', 'LineWidth',1.8, 'DisplayName','Spline cúbico');
plot(f_minimo, Z_minimo, 'rv', 'MarkerFaceColor','r', 'MarkerSize',12, ...
     'DisplayName', sprintf('Mínimo: f* = %.1f Hz, |Z| = %.2f Ω', f_minimo, Z_minimo));
xline(f_minimo, 'r--', sprintf('f* = %.1f Hz', f_minimo), ...
      'LabelVerticalAlignment','bottom', 'FontSize',9);
grid on; xlabel('f (Hz)'); ylabel('|Z| (Ω)');
title('Spline cúbico — Localización precisa del mínimo','FontSize',12);
legend('Location','northeast','FontSize',9); hold off;

%% Figura 2: Primera derivada d|Z|/df vs f
figure('Name','C1 — Primera derivada','Position',[80 80 950 450]);
subplot(1,1,1);
plot(f_fine, dZ_df_fine, 'b-', 'LineWidth',2, 'DisplayName','d|Z|/df (spline)'); hold on;
plot(f, dZ_df_nodos, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5, ...
     'DisplayName','Derivada en nodos');
yline(0, 'k--', 'LineWidth',1.2);
plot(f_minimo, 0, 'rv', 'MarkerFaceColor','r', 'MarkerSize',12, ...
     'DisplayName', sprintf('Cruce por cero: f* = %.1f Hz', f_minimo));
xline(f_minimo, 'r--', sprintf('f* = %.1f Hz', f_minimo), ...
      'LabelVerticalAlignment','bottom','FontSize',9);

% Sombrear zona negativa (decreciente)
idx_neg = dZ_df_fine < 0;
fill([f_fine(idx_neg), fliplr(f_fine(idx_neg))], ...
     [dZ_df_fine(idx_neg), zeros(1,sum(idx_neg))], ...
     'b', 'FaceAlpha',0.08, 'EdgeColor','none', 'DisplayName','Zona decreciente');
% Sombrear zona positiva (creciente)
idx_pos = dZ_df_fine > 0;
fill([f_fine(idx_pos), fliplr(f_fine(idx_pos))], ...
     [dZ_df_fine(idx_pos), zeros(1,sum(idx_pos))], ...
     'r', 'FaceAlpha',0.08, 'EdgeColor','none', 'DisplayName','Zona creciente');

grid on; xlabel('f (Hz)'); ylabel('d|Z|/df  (Ω/Hz)');
title('C1 — Primera derivada analítica del spline cúbico','FontSize',12);
legend('Location','northwest','FontSize',9); hold off;

%% Figura 3: Segunda derivada d²|Z|/df²
figure('Name','C2 — Segunda derivada','Position',[100 100 950 420]);
plot(f_fine, d2Z_df2_fine, 'm-', 'LineWidth',2, 'DisplayName','d²|Z|/df²'); hold on;
plot(f, d2Z_df2_nodos, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5, ...
     'DisplayName','Segunda deriv. en nodos');
plot(f_minimo, d2Z_minimo, 'r^', 'MarkerFaceColor','r', 'MarkerSize',12, ...
     'DisplayName', sprintf('En mínimo: %.6f Ω/Hz²', d2Z_minimo));
yline(0,'k--','LineWidth',1);
xline(f_minimo,'r--',sprintf('f* = %.1f Hz',f_minimo), ...
      'LabelVerticalAlignment','bottom','FontSize',9);
grid on; xlabel('f (Hz)'); ylabel('d²|Z|/df²  (Ω/Hz²)');
title('C2 — Segunda derivada: estabilidad del mínimo','FontSize',12);
legend('Location','northeast','FontSize',9); hold off;

%% Figura 4: Panel resumen (4 subplots)
figure('Name','C — Panel resumen','Position',[50 30 1200 700]);

subplot(2,2,1);
plot(f_fine, Z_fine,'b-','LineWidth',1.8); hold on;
plot(f, Z,'ko','MarkerFaceColor','k','MarkerSize',4);
plot(f_minimo, Z_minimo,'rv','MarkerFaceColor','r','MarkerSize',10);
xline(f_minimo,'r--'); grid on;
xlabel('f (Hz)'); ylabel('|Z| (Ω)');
title('Spline cúbico |Z|(f)');

subplot(2,2,2);
plot(f_fine, dZ_df_fine,'b-','LineWidth',1.8); hold on;
plot(f, dZ_df_nodos,'ko','MarkerFaceColor','k','MarkerSize',4);
yline(0,'k--'); plot(f_minimo,0,'rv','MarkerFaceColor','r','MarkerSize',10);
grid on; xlabel('f (Hz)'); ylabel('d|Z|/df (Ω/Hz)');
title('Primera derivada');

subplot(2,2,3);
plot(f_fine, d2Z_df2_fine,'m-','LineWidth',1.8); hold on;
plot(f, d2Z_df2_nodos,'ko','MarkerFaceColor','k','MarkerSize',4);
yline(0,'k--');
plot(f_minimo, d2Z_minimo,'r^','MarkerFaceColor','r','MarkerSize',10);
grid on; xlabel('f (Hz)'); ylabel('d²|Z|/df² (Ω/Hz²)');
title('Segunda derivada');

subplot(2,2,4);
% Zoom zona del mínimo ±300 Hz
idx_zoom = f_fine >= 400 & f_fine <= 1100;
yyaxis left;
plot(f_fine(idx_zoom), Z_fine(idx_zoom),'b-','LineWidth',2);
ylabel('|Z| (Ω)'); 
yyaxis right;
plot(f_fine(idx_zoom), dZ_df_fine(idx_zoom),'r-','LineWidth',2);
ylabel('d|Z|/df (Ω/Hz)');
yline(0,'k--'); xline(f_minimo,'k:');
grid on; xlabel('f (Hz)');
title(sprintf('Zoom mínimo (f* ≈ %.0f Hz)', f_minimo));

sgtitle('PARTE C — Derivación analítica del spline cúbico','FontWeight','bold','FontSize',13);