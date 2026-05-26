%% =========================================================
%  PARTE 3: Raíces por Cambio de Signo y Bisección
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

%% ── Spline cúbico natural (función de referencia) ────────
cs_V   = csape(f, V, 'variational');
Vspl   = @(fi) fnval(cs_V, fi);     % función anónima para evaluar

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 1: Identificar intervalos con cambio de signo
%  ══════════════════════════════════════════════════════════
fprintf('\n══════ PREGUNTA 1: Cambios de signo en V(f) ════════════\n\n');
fprintf('%-10s %-10s %-10s %-15s\n','i','f(kHz)','V(V)','Cambio de signo');

cruces = [];   % guarda [a, b] de cada intervalo con cruce
for i = 1:length(f)-1
    flag = '';
    if V(i)*V(i+1) < 0
        flag = sprintf('  ◄ CRUCE #%d: [%.1f, %.1f]', length(cruces)+1, f(i), f(i+1));
        cruces = [cruces; f(i), f(i+1)];
    end
    fprintf('%-10d %-10.1f %-10.3f %s\n', i, f(i), V(i), flag);
end
fprintf('%-10d %-10.1f %-10.3f\n', length(f), f(end), V(end));

fprintf('\n► Se encontraron %d intervalos con cambio de signo:\n', size(cruces,1));
for k = 1:size(cruces,1)
    i_a = find(f == cruces(k,1));
    i_b = find(f == cruces(k,2));
    fprintf('  Cruce #%d: [%.1f, %.1f] kHz   V(%.1f)=%.3f  V(%.1f)=%.3f\n', ...
        k, cruces(k,1), cruces(k,2), ...
        cruces(k,1), V(i_a), cruces(k,2), V(i_b));
end

%% ══════════════════════════════════════════════════════════
%  FUNCIÓN DE BISECCIÓN (genérica)
%  ══════════════════════════════════════════════════════════
function [raiz, tabla] = biseccion(func, a, b, tol, max_iter, nombre)
    % Retorna la raíz y una tabla con el historial de iteraciones
    tabla = [];
    fprintf('\n══════ BISECCIÓN — %s ════════════════════════════════\n', nombre);
    fprintf('  Intervalo inicial: [%.4f, %.4f]   |b-a|=%.4f\n', a, b, b-a);
    fprintf('  Tolerancia: %.4f kHz\n\n', tol);
    fprintf('  %-5s %-10s %-10s %-12s %-12s %-10s %-8s\n', ...
            'It','a','b','c=(a+b)/2','V(c)','|b-a|','Decisión');
    fprintf('  %s\n', repmat('-',70,1));

    for k = 1:max_iter
        c  = (a + b) / 2;
        Vc = func(c);
        Va = func(a);
        ba = b - a;

        if Va * Vc < 0
            decision = 'b ← c';
            b = c;
        elseif Va * Vc > 0
            decision = 'a ← c';
            a = c;
        else
            decision = 'EXACTO';
        end

        tabla = [tabla; k, a, b, c, Vc, ba];
        fprintf('  %-5d %-10.5f %-10.5f %-12.5f %-12.6f %-10.5f %-8s\n', ...
                k, a, b, c, Vc, ba, decision);

        if ba < tol
            fprintf('\n  ✔ Convergencia en iteración %d\n', k);
            break;
        end
    end
    raiz = (a + b) / 2;
    fprintf('  ► Raíz aproximada: %.5f kHz   |V(raíz)|=%.2e\n', raiz, abs(func(raiz)));
end

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 2: Bisección — Primer cruce [55.0, 57.5]
%  ══════════════════════════════════════════════════════════
tol = 0.1;   % resolución del instrumento [kHz]

[raiz1_bis, tabla1] = biseccion(Vspl, cruces(1,1), cruces(1,2), tol, 50, 'Cruce #1');

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 3: Bisección — Segundo cruce [62.5, 65.0]
%  ══════════════════════════════════════════════════════════
[raiz2_bis, tabla2] = biseccion(Vspl, cruces(2,1), cruces(2,2), tol, 50, 'Cruce #2');

%% ══════════════════════════════════════════════════════════
%  PREGUNTA 4: Refinamiento con Spline (fzero)
%  ══════════════════════════════════════════════════════════
raiz1_spl = fzero(Vspl, [cruces(1,1), cruces(1,2)]);
raiz2_spl = fzero(Vspl, [cruces(2,1), cruces(2,2)]);

fprintf('\n══════ PREGUNTA 4: Refinamiento con fzero (spline) ════\n');
fprintf('  Raíz 1 — Bisección: %.5f kHz   Spline fzero: %.6f kHz\n', raiz1_bis, raiz1_spl);
fprintf('  Raíz 2 — Bisección: %.5f kHz   Spline fzero: %.6f kHz\n', raiz2_bis, raiz2_spl);
fprintf('  Diferencia raíz 1: %.5f kHz  (%.2f%%)\n', ...
    abs(raiz1_bis-raiz1_spl), abs(raiz1_bis-raiz1_spl)/raiz1_spl*100);
fprintf('  Diferencia raíz 2: %.5f kHz  (%.2f%%)\n', ...
    abs(raiz2_bis-raiz2_spl), abs(raiz2_bis-raiz2_spl)/raiz2_spl*100);

fprintf('\n══════ TABLA COMPARATIVA FINAL ═════════════════════════\n');
fprintf('  %-25s %-15s %-15s %-15s\n','Raíz','Bisección (kHz)','Spline (kHz)','Dif. (kHz)');
fprintf('  %-25s %-15.5f %-15.6f %-15.5f\n','1er cruce (V: + → −)',raiz1_bis,raiz1_spl,abs(raiz1_bis-raiz1_spl));
fprintf('  %-25s %-15.5f %-15.6f %-15.5f\n','2do cruce (V: − → +)',raiz2_bis,raiz2_spl,abs(raiz2_bis-raiz2_spl));

%% ══════════════════════════════════════════════════════════
%  GRÁFICAS
%  ══════════════════════════════════════════════════════════
c_data  = [0.15 0.35 0.70];
c_spl   = [0.10 0.65 0.40];
c_bis   = [0.85 0.33 0.10];
c_raiz  = [0.90 0.10 0.20];
c_zero  = [0.40 0.40 0.40];
c_zona  = [1.00 0.90 0.80];   % fondo zona de cruce

f_fine  = linspace(min(f), max(f), 2000);
V_fine  = fnval(cs_V, f_fine);

%% ── Figura 1: V(f) completo con raíces marcadas ──────────
figure('Name','Raíces de V(f)','Color','w','Position',[50 50 1000 480]);

% Sombrear zonas de cruce
patch([54 58 58 54], [-0.12 -0.12 0.12 0.12], c_zona, ...
      'EdgeColor','none','FaceAlpha',0.5); hold on;
patch([61 66 66 61], [-0.12 -0.12 0.12 0.12], c_zona, ...
      'EdgeColor','none','FaceAlpha',0.5);

plot(f_fine, V_fine, '-', 'Color', c_spl, 'LineWidth', 2.2);
scatter(f, V, 40, c_data, 'filled', 'MarkerEdgeColor','w');
yline(0, '-', 'Color', c_zero, 'LineWidth', 1.5);

% Raíces bisección
xline(raiz1_bis,'--','Color',c_bis,'LineWidth',1.5,'Alpha',0.8);
xline(raiz2_bis,'--','Color',c_bis,'LineWidth',1.5,'Alpha',0.8);

% Raíces spline
scatter(raiz1_spl, 0, 150, c_raiz, 'p', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5);
scatter(raiz2_spl, 0, 150, c_raiz, 'p', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5);

% Anotaciones
text(raiz1_bis+0.3, 0.02, sprintf('Bis: %.3f kHz', raiz1_bis), ...
    'Color',c_bis,'FontSize',8.5,'FontWeight','bold');
text(raiz1_spl+0.3,-0.025,sprintf('Spl: %.4f kHz', raiz1_spl), ...
    'Color',c_raiz,'FontSize',8.5,'FontWeight','bold');
text(raiz2_bis+0.3, 0.02, sprintf('Bis: %.3f kHz', raiz2_bis), ...
    'Color',c_bis,'FontSize',8.5,'FontWeight','bold');
text(raiz2_spl+0.3,-0.025,sprintf('Spl: %.4f kHz', raiz2_spl), ...
    'Color',c_raiz,'FontSize',8.5,'FontWeight','bold');

text(55.2, 0.10, 'Zona cruce #1', 'FontSize',8,'Color',[0.6 0.4 0.0],'FontWeight','bold');
text(62.0, 0.10, 'Zona cruce #2', 'FontSize',8,'Color',[0.6 0.4 0.0],'FontWeight','bold');

xlabel('Frecuencia f (kHz)', 'FontSize',12,'FontWeight','bold');
ylabel('V (V)',               'FontSize',12,'FontWeight','bold');
title('V(f) — Cruces por cero | Bisección vs Spline refinado', ...
      'FontSize',13,'FontWeight','bold');
legend('Zona cruce #1','Zona cruce #2','Spline cúbico','Datos medidos', ...
       'V = 0','Bisección','','Raíces spline (fzero)', ...
       'Location','northwest','FontSize',8.5);
grid on; grid minor;
xlim([min(f)-1, max(f)+1]); ylim([-0.12, 1.65]);
set(gca,'FontSize',10,'Box','on'); hold off;

%% ── Figura 2: Zoom cruce #1 + convergencia bisección ─────
figure('Name','Zoom Cruce #1 + Bisección','Color','w','Position',[60 60 1050 440]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel izquierdo: zoom cruce #1
nexttile;
f_z1  = linspace(54.5, 57.5, 500);
V_z1  = fnval(cs_V, f_z1);

patch([cruces(1,1) cruces(1,2) cruces(1,2) cruces(1,1)], ...
      [-0.07 -0.07 0.03 0.03], c_zona, 'EdgeColor','none','FaceAlpha',0.6); hold on;
plot(f_z1, V_z1, '-', 'Color', c_spl, 'LineWidth', 2.2);
scatter([cruces(1,1), cruces(1,2)], [V(f==cruces(1,1)), V(f==cruces(1,2))], ...
        70, c_data, 'filled','MarkerEdgeColor','w');
yline(0,'-','Color',c_zero,'LineWidth',1.5);
xline(raiz1_bis,'--','Color',c_bis,'LineWidth',1.8);
scatter(raiz1_spl, 0, 130, c_raiz, 'p','filled','MarkerEdgeColor','k');

% Iteraciones de bisección (puntos medios)
c_iter = tabla1(:,4);
scatter(c_iter, zeros(size(c_iter)), 40, ...
    linspace(0.1,0.9,length(c_iter))', 'o', 'filled', ...
    'MarkerEdgeColor','k','LineWidth',0.5, 'CData', ...
    repmat([0.85 0.33 0.10], length(c_iter), 1));

text(raiz1_bis+0.05, -0.055, sprintf('Bis: %.4f', raiz1_bis), ...
    'Color',c_bis,'FontSize',8,'FontWeight','bold');
text(raiz1_spl+0.05,  0.010, sprintf('Spl: %.5f', raiz1_spl), ...
    'Color',c_raiz,'FontSize',8,'FontWeight','bold');

xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('V (V)',  'FontSize',11,'FontWeight','bold');
title('Zoom cruce #1  [55.0 – 57.5 kHz]','FontSize',11,'FontWeight','bold');
legend('Intervalo inicial','Spline','Datos extremos','V=0', ...
       'Bisección','Raíz spline','Puntos medios','FontSize',8,'Location','northeast');
grid on; grid minor; hold off;

% Panel derecho: zoom cruce #2
nexttile;
f_z2  = linspace(61.5, 65.5, 500);
V_z2  = fnval(cs_V, f_z2);

patch([cruces(2,1) cruces(2,2) cruces(2,2) cruces(2,1)], ...
      [-0.05 -0.05 0.03 0.03], c_zona, 'EdgeColor','none','FaceAlpha',0.6); hold on;
plot(f_z2, V_z2, '-', 'Color', c_spl, 'LineWidth', 2.2);
scatter([cruces(2,1), cruces(2,2)], [V(f==cruces(2,1)), V(f==cruces(2,2))], ...
        70, c_data, 'filled','MarkerEdgeColor','w');
yline(0,'-','Color',c_zero,'LineWidth',1.5);
xline(raiz2_bis,'--','Color',c_bis,'LineWidth',1.8);
scatter(raiz2_spl, 0, 130, c_raiz, 'p','filled','MarkerEdgeColor','k');

c_iter2 = tabla2(:,4);
scatter(c_iter2, zeros(size(c_iter2)), 40, ...
    repmat([0.85 0.33 0.10], length(c_iter2), 1), 'o','filled','MarkerEdgeColor','k');

text(raiz2_bis+0.05, -0.040, sprintf('Bis: %.4f', raiz2_bis), ...
    'Color',c_bis,'FontSize',8,'FontWeight','bold');
text(raiz2_spl+0.05,  0.008, sprintf('Spl: %.5f', raiz2_spl), ...
    'Color',c_raiz,'FontSize',8,'FontWeight','bold');

xlabel('f (kHz)','FontSize',11,'FontWeight','bold');
ylabel('V (V)',  'FontSize',11,'FontWeight','bold');
title('Zoom cruce #2  [62.5 – 65.0 kHz]','FontSize',11,'FontWeight','bold');
legend('Intervalo inicial','Spline','Datos extremos','V=0', ...
       'Bisección','Raíz spline','Puntos medios','FontSize',8,'Location','northwest');
grid on; grid minor; hold off;

%% ── Figura 3: Convergencia de la bisección ───────────────
figure('Name','Convergencia Bisección','Color','w','Position',[70 70 900 400]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(tabla1(:,1), abs(tabla1(:,5)), 'o-', ...
    'Color',c_bis,'LineWidth',2,'MarkerFaceColor',c_bis,'MarkerSize',7);
hold on;
yline(1e-3,'--','Color',[0.5 0.5 0.5],'LineWidth',1,'Label','|V|=0.001');
xlabel('Iteración','FontSize',11,'FontWeight','bold');
ylabel('|V(c)| (V)','FontSize',11,'FontWeight','bold');
title('Convergencia Bisección — Cruce #1','FontSize',11,'FontWeight','bold');
grid on; grid minor; hold off;

nexttile;
semilogy(tabla2(:,1), abs(tabla2(:,5)), 's-', ...
    'Color',[0.10 0.65 0.40],'LineWidth',2, ...
    'MarkerFaceColor',[0.10 0.65 0.40],'MarkerSize',7);
hold on;
yline(1e-3,'--','Color',[0.5 0.5 0.5],'LineWidth',1,'Label','|V|=0.001');
xlabel('Iteración','FontSize',11,'FontWeight','bold');
ylabel('|V(c)| (V)','FontSize',11,'FontWeight','bold');
title('Convergencia Bisección — Cruce #2','FontSize',11,'FontWeight','bold');
grid on; grid minor; hold off;

%% ── Figura 4: Ancho del intervalo por iteración ──────────
figure('Name','Reducción del intervalo','Color','w','Position',[80 80 900 400]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(tabla1(:,1), tabla1(:,6), 'o-', 'Color',c_bis,'LineWidth',2, ...
    'MarkerFaceColor',c_bis,'MarkerSize',7); hold on;
yline(tol,'--r','LineWidth',1.5,'Label',sprintf('Tol = %.1f kHz',tol));
xlabel('Iteración','FontSize',11,'FontWeight','bold');
ylabel('|b - a| (kHz)','FontSize',11,'FontWeight','bold');
title('Reducción del intervalo — Cruce #1','FontSize',11,'FontWeight','bold');
grid on; grid minor; hold off;

nexttile;
plot(tabla2(:,1), tabla2(:,6), 's-', 'Color',[0.10 0.65 0.40],'LineWidth',2, ...
    'MarkerFaceColor',[0.10 0.65 0.40],'MarkerSize',7); hold on;
yline(tol,'--r','LineWidth',1.5,'Label',sprintf('Tol = %.1f kHz',tol));
xlabel('Iteración','FontSize',11,'FontWeight','bold');
ylabel('|b - a| (kHz)','FontSize',11,'FontWeight','bold');
title('Reducción del intervalo — Cruce #2','FontSize',11,'FontWeight','bold');
grid on; grid minor; hold off;

fprintf('\n✔ Gráficas Parte 3 generadas correctamente.\n');