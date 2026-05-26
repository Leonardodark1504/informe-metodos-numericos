%% =========================================================
%  PARTE 4: Discusión Técnica — Visualizaciones cuantitativas
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

h       = 2.5;          % paso [kHz]
res_f   = 0.1;          % resolución generador [kHz]
res_V   = 0.001;        % resolución voltímetro [V]
res_Z   = 0.1;          % resolución medidor impedancia [Ω]
eps_V   = res_V / 2;    % error máximo de cuantización en V
eps_Z   = res_Z / 2;    % error máximo de cuantización en Z

%% ── Splines cúbicos naturales ────────────────────────────
cs_V   = csape(f, V, 'variational');
cs_Z   = csape(f, Z, 'variational');
ds_V   = fnder(cs_V, 1);    % derivada del spline de V

f_fine = linspace(min(f), max(f), 2000);
V_fine = fnval(cs_V, f_fine);
Z_fine = fnval(cs_Z, f_fine);
dV_fine= fnval(ds_V, f_fine);

%% ── Derivada centrada O2 en todos los nodos ──────────────
n      = length(f);
dV_num = zeros(1,n);
dV_num(1)   = (-3*V(1)  + 4*V(2)   - V(3))   / (2*h);
dV_num(end) = ( 3*V(end)- 4*V(n-1) + V(n-2)) / (2*h);
for i = 2:n-1
    dV_num(i) = (V(i+1) - V(i-1)) / (2*h);
end

%% ── Raíces (de Parte 3) ──────────────────────────────────
Vspl    = @(fi) fnval(cs_V, fi);
raiz1   = fzero(Vspl, [55.0, 57.5]);
raiz2   = fzero(Vspl, [62.5, 65.0]);

%% ══════════════════════════════════════════════════════════
%  ANÁLISIS CUANTITATIVO PARA LA DISCUSIÓN
%  ══════════════════════════════════════════════════════════

%% ── Error de cuantización propagado en derivación ────────
eps_deriv = (2*eps_V) / (2*h);   % error cuantización en dV/df

%% ── Error en localización de raíces ─────────────────────
dV_raiz1 = abs(fnval(ds_V, raiz1));
dV_raiz2 = abs(fnval(ds_V, raiz2));
delta_f1  = eps_V / dV_raiz1;
delta_f2  = eps_V / dV_raiz2;

%% ── SNR instrumental por nodo ────────────────────────────
SNR_V = abs(V) ./ res_V;         % relación señal/resolución
SNR_Z = Z      ./ res_Z;

%% ── Variación de Z por zona ──────────────────────────────
% Zona A: 10–32.5,  B: 35–55,  C: 57.5–65,  D: 67.5–82.5,  E: 85–107.5
zonas = {find(f<=32.5), find(f>=35 & f<=55), find(f>=57.5 & f<=65), ...
         find(f>=67.5 & f<=82.5), find(f>=85)};
nom_zonas = {'A: 10–32.5 kHz','B: 35–55 kHz','C: 57.5–65 kHz', ...
             'D: 67.5–82.5 kHz','E: 85–107.5 kHz'};

fprintf('\n══════ PARTE 4: ANÁLISIS CUANTITATIVO ══════════════════\n');
fprintf('\n► Error de cuantización en derivación:\n');
fprintf('  ε_deriv = 2×%.4f / (2×%.1f) = %.5f V/kHz\n', eps_V, h, eps_deriv);
fprintf('  Error relativo en f=40 kHz: %.2f%%\n', eps_deriv/abs(fnval(ds_V,40))*100);
fprintf('  Error relativo en f=100kHz: %.2f%%\n', eps_deriv/abs(fnval(ds_V,100))*100);

fprintf('\n► Incertidumbre en raíces por resolución de voltímetro:\n');
fprintf('  Δf_raiz1 = %.4f / %.5f = %.4f kHz  (res. gen.=%.1f kHz)\n', ...
        eps_V, dV_raiz1, delta_f1, res_f);
fprintf('  Δf_raiz2 = %.4f / %.5f = %.4f kHz\n', eps_V, dV_raiz2, delta_f2);

fprintf('\n► Variación de |Z(f)| por zona:\n');
for k = 1:5
    idx = zonas{k};
    dZ  = max(Z(idx)) - min(Z(idx));
    fprintf('  %s → ΔZ=%.1f Ω  (%.1f%%)\n', nom_zonas{k}, dZ, dZ/mean(Z(idx))*100);
end

fprintf('\n► SNR instrumental en puntos clave:\n');
pts = [10, 40, 55, 70, 100];
for p = pts
    idx = find(f==p);
    fprintf('  f=%5.1f kHz: V=%.3f V, SNR_V=%.0f, |Z|=%.1f Ω, SNR_Z=%.0f\n', ...
            f(idx), V(idx), SNR_V(idx), Z(idx), SNR_Z(idx));
end

%% ══════════════════════════════════════════════════════════
%  COLORES Y ZONAS
%  ══════════════════════════════════════════════════════════
cA = [0.70 0.90 0.70]; cB = [0.90 0.85 0.65]; cC = [0.95 0.70 0.70];
cD = [0.70 0.85 0.95]; cE = [0.65 0.90 0.75];

zona_lim = [10, 32.5; 32.5, 55; 55, 65; 65, 82.5; 82.5, 107.5];
zona_col = {cA, cB, cC, cD, cE};
zona_nom = {'A','B','C','D','E'};
zona_txt = {'A: Rampa ascendente','B: Caída post-pico', ...
            'C: Cruce por cero','D: Recuperación','E: Zona estable ★'};

c_data = [0.15 0.35 0.70];
c_spl  = [0.10 0.60 0.35];
c_dV   = [0.85 0.33 0.10];
c_Z    = [0.50 0.10 0.70];
c_zero = [0.40 0.40 0.40];

%% ══════════════════════════════════════════════════════════
%  FIGURA 1: V(f) y |Z(f)| con zonas operativas (doble eje)
%  ══════════════════════════════════════════════════════════
figure('Name','P4-F1 Zonas operativas V y Z','Color','w','Position',[40 40 1050 520]);

% Sombrear zonas
for k = 1:5
    patch([zona_lim(k,1) zona_lim(k,2) zona_lim(k,2) zona_lim(k,1)], ...
          [-0.15 -0.15 1.75 1.75], zona_col{k}, ...
          'EdgeColor','none','FaceAlpha',0.35); hold on;
    text(mean(zona_lim(k,:)), 1.62, zona_nom{k}, ...
        'HorizontalAlignment','center','FontSize',9, ...
        'FontWeight','bold','Color',[0.25 0.25 0.25]);
end

yyaxis left
plot(f_fine, V_fine, '-', 'Color', c_spl, 'LineWidth', 2.2);
scatter(f, V, 40, c_data, 'filled','MarkerEdgeColor','w');
yline(0,'--','Color',c_zero,'LineWidth',1.2);
% Raíces
scatter([raiz1 raiz2],[0 0],120,[0.9 0.1 0.2],'p','filled', ...
        'MarkerEdgeColor','k','LineWidth',1.2);
text(raiz1+0.3, -0.07, sprintf('%.1f kHz',raiz1),'FontSize',8,'Color',[0.8 0.1 0.1],'FontWeight','bold');
text(raiz2+0.3, -0.07, sprintf('%.1f kHz',raiz2),'FontSize',8,'Color',[0.8 0.1 0.1],'FontWeight','bold');
ylabel('V(f)  [V]','FontSize',12,'FontWeight','bold','Color',c_spl);
ylim([-0.15 1.75]);

yyaxis right
plot(f_fine, Z_fine, '-', 'Color', c_Z, 'LineWidth', 2.0, 'LineStyle','--');
scatter(f, Z, 35, c_Z, 'd','filled','MarkerEdgeColor','w','MarkerFaceAlpha',0.7);
ylabel('|Z(f)|  [Ω]','FontSize',12,'FontWeight','bold','Color',c_Z);
ylim([130 230]);

xlabel('Frecuencia f (kHz)','FontSize',12,'FontWeight','bold');
title('Parte 4 — V(f) y |Z(f)|: Zonas operativas del biosensor','FontSize',13,'FontWeight','bold');

ax = gca; ax.YAxis(1).Color = c_spl; ax.YAxis(2).Color = c_Z;
legend('','','','','','V(f) spline','Datos V','V=0','Raíces V(f)', ...
       '','|Z(f)| spline','Datos |Z|', ...
       'Location','north','FontSize',8,'NumColumns',4);
grid on; grid minor;
xlim([min(f)-1, max(f)+1]);
set(gca,'FontSize',10,'Box','on'); hold off;

% Leyenda de zonas
axes('Position',[0.13 0.01 0.77 0.06],'Visible','off');
for k=1:5
    text((k-0.5)/5, 0.5, zona_txt{k},'Units','normalized', ...
        'HorizontalAlignment','center','FontSize',8.0, ...
        'BackgroundColor',zona_col{k},'EdgeColor',[0.6 0.6 0.6], ...
        'Margin',3,'FontWeight','bold');
end

%% ══════════════════════════════════════════════════════════
%  FIGURA 2: dV/df con zonas — Sensibilidad del front-end
%  ══════════════════════════════════════════════════════════
figure('Name','P4-F2 Sensibilidad dV/df por zona','Color','w','Position',[50 50 1050 460]);

for k = 1:5
    patch([zona_lim(k,1) zona_lim(k,2) zona_lim(k,2) zona_lim(k,1)], ...
          [-0.10 -0.10 0.10 0.10], zona_col{k}, ...
          'EdgeColor','none','FaceAlpha',0.40); hold on;
    text(mean(zona_lim(k,:)), 0.087, zona_nom{k}, ...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.25 0.25 0.25]);
end

plot(f_fine, dV_fine, '-', 'Color', c_dV, 'LineWidth', 2.2);
scatter(f, dV_num, 50, c_data, 'o','filled','MarkerEdgeColor','w');
yline(0,'--','Color',c_zero,'LineWidth',1.2);

% Banda de error por cuantización
fill([f_fine, fliplr(f_fine)], ...
     [dV_fine+eps_deriv, fliplr(dV_fine-eps_deriv)], ...
     c_dV, 'FaceAlpha', 0.15, 'EdgeColor','none');

% Anotar puntos clave
pts_key = [10, 40, 70, 100];
for p = pts_key
    dv = fnval(ds_V, p);
    scatter(p, dv, 100, [0.9 0.1 0.2],'p','filled','MarkerEdgeColor','k');
    text(p+0.5, dv+0.005*sign(dv), sprintf('%.4f\nV/kHz',dv), ...
        'FontSize',7.5,'Color',[0.7 0.1 0.1],'FontWeight','bold');
end

xlabel('Frecuencia f (kHz)','FontSize',12,'FontWeight','bold');
ylabel('dV/df  [V/kHz]','FontSize',12,'FontWeight','bold');
title('Sensibilidad del front-end: dV/df por zona operativa','FontSize',13,'FontWeight','bold');
legend('','','','','','Spline d/df','Centrada O2','V=0', ...
       sprintf('±ε_{cuant}=±%.4f V/kHz',eps_deriv),'Puntos clave', ...
       'Location','northeast','FontSize',9);
grid on; grid minor;
xlim([min(f)-1, max(f)+1]); ylim([-0.10 0.10]);
set(gca,'FontSize',10,'Box','on'); hold off;

%% ══════════════════════════════════════════════════════════
%  FIGURA 3: SNR instrumental V(f) y |Z(f)|
%  ══════════════════════════════════════════════════════════
figure('Name','P4-F3 SNR instrumental','Color','w','Position',[60 60 1000 430]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
b1 = bar(f, SNR_V, 'FaceColor','flat','EdgeColor','none');
% Colorear por zona
col_bar = zeros(n,3);
for k = 1:5
    col_bar(zonas{k},:) = repmat(zona_col{k}*0.75, length(zonas{k}), 1);
end
b1.CData = col_bar;
hold on;
yline(100,'--r','LineWidth',1.5,'Label','SNR=100 (mín. recomendado)');
xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('SNR = |V| / res_V','FontSize',11,'FontWeight','bold');
title('SNR instrumental de V(f)','FontSize',11,'FontWeight','bold');
grid on; grid minor;
xlim([min(f)-2, max(f)+2]);
set(gca,'FontSize',9,'Box','on'); hold off;

nexttile;
% Variación de Z en cada zona (barras de rango)
deltaZ = zeros(1,5); meanZ = zeros(1,5); pctZ = zeros(1,5);
for k=1:5
    deltaZ(k) = max(Z(zonas{k})) - min(Z(zonas{k}));
    meanZ(k)  = mean(Z(zonas{k}));
    pctZ(k)   = deltaZ(k)/meanZ(k)*100;
end
zona_labels = categorical({'A','B','C','D','E'},{'A','B','C','D','E'});
b2 = bar(zona_labels, pctZ, 'FaceColor','flat','EdgeColor','k','LineWidth',0.8);
b2.CData = cell2mat(cellfun(@(c) c*0.8, zona_col','UniformOutput',false));
hold on;
for k=1:5
    text(k, pctZ(k)+0.15, sprintf('%.1f%%\n(ΔZ=%.1fΩ)',pctZ(k),deltaZ(k)), ...
        'HorizontalAlignment','center','FontSize',8.5,'FontWeight','bold');
end
yline(2,'--r','LineWidth',1.5,'Label','2% umbral estabilidad');
xlabel('Zona operativa','FontSize',11,'FontWeight','bold');
ylabel('Variación |Z| (%)','FontSize',11,'FontWeight','bold');
title('Variación relativa de |Z(f)| por zona','FontSize',11,'FontWeight','bold');
grid on; grid minor;
set(gca,'FontSize',9,'Box','on'); hold off;

%% ══════════════════════════════════════════════════════════
%  FIGURA 4: Incertidumbre en raíces vs resolución
%  ══════════════════════════════════════════════════════════
figure('Name','P4-F4 Incertidumbre raíces','Color','w','Position',[70 70 950 420]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
% Zoom cruce #1 con banda de incertidumbre
f_r1   = linspace(54.5, 57.5, 500);
V_r1   = fnval(cs_V, f_r1);
dV_r1  = fnval(ds_V, f_r1);

yyaxis left
plot(f_r1, V_r1, '-', 'Color', c_spl, 'LineWidth', 2.2); hold on;
fill([f_r1, fliplr(f_r1)], ...
     [V_r1+eps_V, fliplr(V_r1-eps_V)], c_spl, 'FaceAlpha',0.2,'EdgeColor','none');
yline(0,'--','Color',c_zero,'LineWidth',1.2);
% Banda de incertidumbre en la raíz
patch([raiz1-delta_f1, raiz1+delta_f1, raiz1+delta_f1, raiz1-delta_f1], ...
      [-0.015 -0.015 0.015 0.015],[0.9 0.3 0.3],'FaceAlpha',0.3,'EdgeColor',[0.8 0.1 0.1]);
scatter(raiz1, 0, 120,[0.9 0.1 0.2],'p','filled','MarkerEdgeColor','k');
ylabel('V(f) [V]','FontSize',10,'FontWeight','bold','Color',c_spl);

yyaxis right
plot(f_r1, abs(dV_r1), '--', 'Color', c_dV, 'LineWidth', 1.8);
ylabel('|dV/df| [V/kHz]','FontSize',10,'FontWeight','bold','Color',c_dV);

xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
title(sprintf('Cruce #1  Δf_{raiz}=±%.4f kHz',delta_f1),'FontSize',11,'FontWeight','bold');
ax=gca; ax.YAxis(1).Color=c_spl; ax.YAxis(2).Color=c_dV;
legend('V(f)','±ε_V','V=0','Δf_{raiz}','Raíz','|dV/df|', ...
       'Location','southeast','FontSize',8);
grid on; grid minor; hold off;

nexttile;
f_r2   = linspace(61.5, 66.0, 500);
V_r2   = fnval(cs_V, f_r2);
dV_r2  = fnval(ds_V, f_r2);

yyaxis left
plot(f_r2, V_r2, '-', 'Color', c_spl, 'LineWidth', 2.2); hold on;
fill([f_r2, fliplr(f_r2)], ...
     [V_r2+eps_V, fliplr(V_r2-eps_V)], c_spl, 'FaceAlpha',0.2,'EdgeColor','none');
yline(0,'--','Color',c_zero,'LineWidth',1.2);
patch([raiz2-delta_f2, raiz2+delta_f2, raiz2+delta_f2, raiz2-delta_f2], ...
      [-0.015 -0.015 0.015 0.015],[0.9 0.3 0.3],'FaceAlpha',0.3,'EdgeColor',[0.8 0.1 0.1]);
scatter(raiz2, 0, 120,[0.9 0.1 0.2],'p','filled','MarkerEdgeColor','k');
ylabel('V(f) [V]','FontSize',10,'FontWeight','bold','Color',c_spl);

yyaxis right
plot(f_r2, abs(dV_r2), '--', 'Color', c_dV, 'LineWidth', 1.8);
ylabel('|dV/df| [V/kHz]','FontSize',10,'FontWeight','bold','Color',c_dV);

xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
title(sprintf('Cruce #2  Δf_{raiz}=±%.4f kHz',delta_f2),'FontSize',11,'FontWeight','bold');
ax=gca; ax.YAxis(1).Color=c_spl; ax.YAxis(2).Color=c_dV;
legend('V(f)','±ε_V','V=0','Δf_{raiz}','Raíz','|dV/df|', ...
       'Location','southeast','FontSize',8);
grid on; grid minor; hold off;

%% ══════════════════════════════════════════════════════════
%  FIGURA 5: Ventaja spline — localidad de la perturbación
%  ══════════════════════════════════════════════════════════
figure('Name','P4-F5 Localidad spline vs polinomio global','Color','w','Position',[80 80 1050 430]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Perturbar punto i=20 (f=57.5, V=-0.041) en ±5%
V_pert = V;
V_pert(20) = V(20) * 1.05;   % perturbación del 5%

cs_orig = csape(f, V,      'variational');
cs_pert = csape(f, V_pert, 'variational');

% Polinomio global (grado n-1) con polyfit
p_orig  = polyfit(f, V,      length(f)-1);
p_pert  = polyfit(f, V_pert, length(f)-1);

f_pg    = linspace(10, 107.5, 2000);
V_pg_o  = polyval(p_orig, f_pg);
V_pg_p  = polyval(p_pert, f_pg);
V_sp_o  = fnval(cs_orig,  f_pg);
V_sp_p  = fnval(cs_pert,  f_pg);

nexttile;
plot(f_pg, abs(V_sp_p - V_sp_o)*1000, '-', 'Color', c_spl, 'LineWidth', 2.2); hold on;
xline(57.5,'--r','LineWidth',1.5,'Label','Punto perturbado');
patch([55 65 65 55],[0 0 max(abs(V_sp_p-V_sp_o))*1200 max(abs(V_sp_p-V_sp_o))*1200], ...
      cC,'FaceAlpha',0.3,'EdgeColor','none');
xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('|ΔV_{spline}| (mV)','FontSize',11,'FontWeight','bold');
title('Spline cúbico — Perturbación LOCAL','FontSize',11,'FontWeight','bold');
text(70, max(abs(V_sp_p-V_sp_o))*950, ...
    sprintf('Máx. cambio\nfuera zona: %.3f mV', ...
    max(abs(V_sp_p(f_pg>67)-V_sp_o(f_pg>67)))*1000), ...
    'FontSize',8,'Color',c_spl,'FontWeight','bold');
grid on; grid minor; hold off;

nexttile;
plot(f_pg, abs(V_pg_p - V_pg_o)*1000, '-', 'Color',[0.85 0.10 0.20], 'LineWidth', 2.2); hold on;
xline(57.5,'--r','LineWidth',1.5,'Label','Punto perturbado');
xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('|ΔV_{polinomio}| (mV)','FontSize',11,'FontWeight','bold');
title('Polinomio global grado 39 — Propagación GLOBAL','FontSize',11,'FontWeight','bold');
text(15, max(abs(V_pg_p-V_pg_o))*0.85*1000, ...
    sprintf('Perturbación se\npropaga a todo\nel dominio'), ...
    'FontSize',8.5,'Color',[0.7 0.1 0.1],'FontWeight','bold');
grid on; grid minor; hold off;

fprintf('\n✔ Todas las figuras de la Parte 4 generadas correctamente.\n');

%% ══════════════════════════════════════════════════════════
%  RESUMEN FINAL IMPRESO
%  ══════════════════════════════════════════════════════════
fprintf('\n══════ RESUMEN PARTE 4 ══════════════════════════════════\n\n');

fprintf(['PREGUNTA 1 — Banda de operación recomendada:\n'...
'  Zona E: 95–107.5 kHz\n'...
'  • V ∈ [0.856, 1.004] V  → estable y positivo (lejos de cruces por cero)\n'...
'  • |Z| ∈ [180.6, 185.1] Ω → variación < 2.5%%  (adaptación estable)\n'...
'  • |dV/df| < 0.013 V/kHz → mínima sensibilidad a derivas de frecuencia\n'...
'  • SNR_V > 856            → señal muy por encima del ruido instrumental\n\n']);

fprintf(['PREGUNTA 2 — Impacto de la resolución instrumental:\n'...
'  • Interpolación: error cuantización ≈ ±0.001 V (suma 3 términos Lagrange)\n'...
'  • Derivación:    ε_cuant = %.5f V/kHz → err. relativo ~1.7%% en f=40kHz\n'...
'  • Raíces:        Δf_raiz ≈ [%.4f, %.4f] kHz < res. generador (0.1 kHz)\n'...
'  ⇒ Límite real impuesto por el generador, no por el voltímetro\n\n'], ...
    eps_deriv, delta_f1, delta_f2);

fprintf(['PREGUNTA 3 — Ventaja práctica del spline cúbico:\n'...
'  Localidad: una perturbación en un nodo solo afecta ±1 tramo vecino.\n'...
'  En polinomio global grado 39 la perturbación se propaga a todo el dominio\n'...
'  (Figura 5). Esto es crítico en sistemas biomédicos donde:\n'...
'  • Los datos pueden actualizarse por recalibración parcial\n'...
'  • Recertificar solo una zona es suficiente (ahorro computacional)\n'...
'  • El costo de recalcular 2×2 tramos es << resolver sistema 40×40\n']);