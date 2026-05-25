clc;
clear;
close all;

%% Datos del experimento
f = [100 120 145 170 200 235 270 310 355 405 ...
     460 520 585 655 730 810 895 985 1080 1180 ...
     1290 1410 1540 1680 1830 1990 2160 2340 2530 2730];

Z = [152.3 149.1 146.8 144.9 142.0 139.5 137.9 136.1 134.8 133.6 ...
     132.7 131.9 131.4 131.1 130.9 131.0 131.3 131.9 132.7 133.8 ...
     135.2 136.9 138.9 141.1 143.5 146.1 149.0 152.2 155.6 159.2];

Zth = 150;

%% Construcción del spline cúbico natural
M = splineNaturalMomentos(f, Z);

%% Función a resolver: g(f) = Z(f) - Zth
g  = @(x) evalSplineNatural(x, f, Z, M) - Zth;
gp = @(x) evalDerivSplineNatural(x, f, Z, M);

%% Búsqueda de intervalos con cambio de signo
intervalos = [];

for i = 1:length(f)-1
    if g(f(i))*g(f(i+1)) < 0
        intervalos = [intervalos; f(i), f(i+1)];
    end
end

disp('Intervalos donde existen raíces:')
disp(intervalos)

%% Método de bisección
tol = 1e-8;
maxIter = 100;

raices_biseccion = [];
iter_biseccion = [];

for k = 1:size(intervalos,1)
    a = intervalos(k,1);
    b = intervalos(k,2);

    [raiz, iter] = biseccion(g, a, b, tol, maxIter);

    raices_biseccion = [raices_biseccion; raiz];
    iter_biseccion = [iter_biseccion; iter];
end

%% Método de Newton-Raphson
% Aproximaciones iniciales cercanas a las raíces
x0 = [110; 2220];

raices_newton = [];
iter_newton = [];

for k = 1:length(x0)
    [raiz, iter] = newtonRaphson(g, gp, x0(k), tol, maxIter);

    raices_newton = [raices_newton; raiz];
    iter_newton = [iter_newton; iter];
end

%% Mostrar resultados
fprintf('\nRESULTADOS DE LA SECCIÓN D\n');
fprintf('-------------------------------------------\n');

for k = 1:length(raices_biseccion)
    fprintf('Raíz %d por bisección:       %.6f Hz | Iteraciones: %d\n', ...
        k, raices_biseccion(k), iter_biseccion(k));
end

fprintf('\n');

for k = 1:length(raices_newton)
    fprintf('Raíz %d por Newton-Raphson:  %.6f Hz | Iteraciones: %d\n', ...
        k, raices_newton(k), iter_newton(k));
end

%% Sensibilidad df/dZ en la raíz alta
raiz_alta = raices_newton(2);

dZdf = gp(raiz_alta);
dfdZ = 1/dZdf;

fprintf('\nSensibilidad en la raíz alta:\n');
fprintf('f = %.6f Hz\n', raiz_alta);
fprintf('dZ/df = %.8f ohm/Hz\n', dZdf);
fprintf('df/dZ = %.4f Hz/ohm\n', dfdZ);

fprintf('\nInterpretación de sensibilidad:\n');
fprintf('Un error de 0.1 ohm desplaza la raíz aproximadamente %.2f Hz.\n', abs(dfdZ)*0.1);
fprintf('Un error de 0.5 ohm desplaza la raíz aproximadamente %.2f Hz.\n', abs(dfdZ)*0.5);
fprintf('Un error de 1.0 ohm desplaza la raíz aproximadamente %.2f Hz.\n', abs(dfdZ)*1.0);

%% Gráfica del spline y el umbral
ff = linspace(min(f), max(f), 1000);
ZZ = evalSplineNatural(ff, f, Z, M);

figure;
plot(f, Z, 'ko', 'MarkerFaceColor', 'k');
hold on;
plot(ff, ZZ, 'b', 'LineWidth', 1.8);
yline(Zth, 'r--', 'LineWidth', 1.5);

plot(raices_newton, Zth*ones(size(raices_newton)), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);

grid on;
xlabel('Frecuencia f (Hz)');
ylabel('|Z| (ohm)');
title('Búsqueda de raíces usando spline cúbico natural');
legend('Datos experimentales', 'Spline cúbico natural', ...
       'Umbral Z_{th} = 150 \Omega', 'Raíces', ...
       'Location', 'best');

%% Gráfica de g(f) = |Z|(f) - 150
gg = ZZ - Zth;

figure;
plot(ff, gg, 'm', 'LineWidth', 1.8);
hold on;
yline(0, 'k--', 'LineWidth', 1.5);
plot(raices_newton, zeros(size(raices_newton)), 'ro', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);

grid on;
xlabel('Frecuencia f (Hz)');
ylabel('g(f) = |Z|(f) - 150');
title('Cruces de la función g(f) con el eje cero');
legend('g(f)', 'Eje cero', 'Raíces', 'Location', 'best');

%% ---------------- FUNCIONES LOCALES ----------------

function M = splineNaturalMomentos(x, y)

    n = length(x);
    h = diff(x);

    A = zeros(n,n);
    b = zeros(n,1);

    % Condiciones naturales
    A(1,1) = 1;
    A(n,n) = 1;

    for i = 2:n-1
        A(i,i-1) = h(i-1);
        A(i,i)   = 2*(h(i-1) + h(i));
        A(i,i+1) = h(i);

        b(i) = 6*((y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1));
    end

    M = A\b;
end

function S = evalSplineNatural(xx, x, y, M)

    S = zeros(size(xx));

    for j = 1:length(xx)

        if xx(j) <= x(1)
            i = 1;
        elseif xx(j) >= x(end)
            i = length(x)-1;
        else
            i = find(x <= xx(j), 1, 'last');
            if i == length(x)
                i = length(x)-1;
            end
        end

        h = x(i+1) - x(i);

        A = (x(i+1) - xx(j))/h;
        B = (xx(j) - x(i))/h;

        S(j) = M(i)*(x(i+1)-xx(j))^3/(6*h) + ...
               M(i+1)*(xx(j)-x(i))^3/(6*h) + ...
               (y(i)-M(i)*h^2/6)*A + ...
               (y(i+1)-M(i+1)*h^2/6)*B;
    end
end

function dS = evalDerivSplineNatural(xx, x, y, M)

    dS = zeros(size(xx));

    for j = 1:length(xx)

        if xx(j) <= x(1)
            i = 1;
        elseif xx(j) >= x(end)
            i = length(x)-1;
        else
            i = find(x <= xx(j), 1, 'last');
            if i == length(x)
                i = length(x)-1;
            end
        end

        h = x(i+1) - x(i);

        dS(j) = -M(i)*(x(i+1)-xx(j))^2/(2*h) + ...
                 M(i+1)*(xx(j)-x(i))^2/(2*h) + ...
                 (y(i+1)-y(i))/h - ...
                 (M(i+1)-M(i))*h/6;
    end
end

function [raiz, iter] = biseccion(g, a, b, tol, maxIter)

    fa = g(a);
    fb = g(b);

    if fa*fb > 0
        error('No hay cambio de signo en el intervalo.');
    end

    for iter = 1:maxIter
        c = (a+b)/2;
        fc = g(c);

        if abs(fc) < tol || abs(b-a)/2 < tol
            raiz = c;
            return;
        end

        if fa*fc < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end
    end

    raiz = (a+b)/2;
end

function [raiz, iter] = newtonRaphson(g, gp, x0, tol, maxIter)

    x = x0;

    for iter = 1:maxIter
        gx = g(x);
        gpx = gp(x);

        if abs(gpx) < 1e-12
            error('La derivada es muy pequeña. Newton puede fallar.');
        end

        xnew = x - gx/gpx;

        if abs(xnew - x) < tol
            raiz = xnew;
            return;
        end

        x = xnew;
    end

    raiz = x;
end