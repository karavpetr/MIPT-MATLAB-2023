%% Задание 3.1
clear

%----Parameters-----------
N = 100;

S_n = zeros([N 1]);
Phi_n = zeros([N 1]);
S_n(1) = 2;
Phi_n(1) = 3;
for k = 2:N
    S_n(k) = S_n(k- 1) + 1/factorial(k);
    Phi_n(k) = 1/k * 1/factorial(k);
end

S = exp(1);
%----End of parameters----

n_values = 2:10;

figure;
hold on;

differences = S - S_n(n_values);

drawing_1 = plot(n_values, differences, 'b-', 'linewidth', 0.5);
drawing_2 = plot(n_values, Phi_n(n_values), 'r-', 'linewidth', 0.5);
legend([drawing_1, drawing_2], {'S - S_n', '\phi_n'});
xlabel('n');
ylabel('value');
title('Difference between partial sums and actual value');

hold off;


%% Задание 3.2
clear

%----Parameters-----------
f_left = @(x) sqrt(x);
f_right = @(x) tan(x);
%----End of parameters----

f = @(x) f_left(x) - f_right(x);

% Plotting functions
figure;
hold on;
x = linspace(0, 12.0, 1000);
ylim([-2 8]);

drawing_1 = plot(x, f_left(x), 'g-');
drawing_2 = plot(x, f_right(x), 'b-');
drawing_3 = plot(x, f(x), 'r-');
drawing_4 = plot(x, zeros(size(x)), 'k-');
legend([drawing_1, drawing_2, drawing_3], {'sqrt(x)', 'tg(x)', 'sqrt(x) - tg(x)'});
xlabel('x');
ylabel('y = f(x)');
title('Отметьте (курсором) границы отрезка для поиска корня');
hold off;

% Finding roots using ginput
[x_values, y_values] = ginput(2);
first_root = fzero(f, x_values(1));
second_root = fzero(f, x_values(2));

hold on;
if abs(f(first_root)) <= 1e-10
    point_1 = plot(first_root, 0, 'k.', 'MarkerSize', 30);
    extra_point_1 = plot(first_root, f_left(first_root), 'k.', 'MarkerSize', 30);
    fprintf('#1: %f is a root\n', first_root);
    point_1_legend = 'some root  (#1)';
else
    point_1 = plot(first_root, 0, 'k.', 'MarkerSize', 30);
    fprintf('#1: %f is NOT a root\n', first_root);
    point_1_legend = 'NOT a root (#1)';
end

if abs(f(second_root)) <= 1e-10
    point_2 = plot(second_root, 0, 'ko', 'MarkerSize', 8);
    extra_point_2 = plot(second_root, f_left(second_root), 'ko', 'MarkerSize', 8);
    fprintf('#2: %f is a root\n', second_root);
    point_2_legend = 'some root  (#2)';
else
    point_2 = plot(second_root, 0, 'ko', 'MarkerSize', 8);
    fprintf('#2: %f is NOT a root\n', second_root);
    point_2_legend = 'NOT a root (#2)';
end
legend([drawing_1, drawing_2, drawing_3, point_1, point_2], {'sqrt(x)', 'tg(x)', 'sqrt(x) - tg(x)', point_1_legend, point_2_legend});
title('');
hold off;


%% Задание 3.3
clear

%----Parameters-----------
f = @(x) fillmissing(sqrt(abs(x)) .* sin(1./x.^2) .* (x ~= 0), 'constant', 0);
%----End of parameters----

figure;
hold on;
N = 1000;
xVec = linspace(-1, 1, N);
yVec = zeros([1 N]);
for i = 1:N
    yVec(i) = fzero(f, xVec(i));
end
drawing_1 = plot(xVec, f(xVec), 'b-');
drawing_2 = plot(xVec, yVec, 'k-');
drawing_3 = plot(xVec, zeros([1 N]), 'g-');

legend([drawing_1, drawing_2, drawing_3], {'function', 'roots', 'y = 0'});
xlabel('x');
ylabel('y = f(x)');
hold off;

disp(yVec);


%% Задание 3.4
clear

%----Parameters-----------
AMat = rand(3, 3);
%----End of parameters----


% Method 1: partial sum
N = 50;
expA1Mat = eye(size(AMat));
for k = 1:N
    expA1Mat = expA1Mat + (AMat^k / factorial(k));
end

% Method 2: Chaucy ODE
ode_function = @(t, XMat) reshape(AMat * reshape(XMat, 3, []), [], 1);
[t, XMat_values] = ode45(ode_function, [0 1], eye(3));
expA2Mat = reshape(XMat_values(end, :), 3, 3);

% Method 3: in-built expm
expA3Mat = expm(AMat);

% Comparing results
disp("Method 1 (partial sum):")
disp(expA1Mat);

disp("Method 2 (Chauchy ODE):")
disp(expA2Mat);

disp("Mathod 3: (In-built function)")
disp(expA3Mat);

fprintf("Difference between Partial-sum-method and Cauchy-ode-method: %f\n", norm(expA1Mat - expA2Mat, 'fro'));
fprintf("Difference between Partial-sum-method and In-built-method: %f\n", norm(expA1Mat - expA3Mat, 'fro'));
fprintf("Difference between Cauchy-ode-method and In-built-method: %f\n", norm(expA2Mat - expA3Mat, 'fro'));


%% Задание 3.5
clear

%----Parameters-----------
T = 20;
dt = 0.01;
population = 1000000;
vaccinated_speed = 100;

% Первая ситуация
vaccinated = 0.9;
infected = 0.02;
isolated = 0.01;

X_start_1 = [(1 - vaccinated - infected - isolated) * population; vaccinated * population; infected * population; isolated * population; vaccinated_speed];
params_for_first_situation = [1e-5; 3 * 1e-6; 1.5 * 1e-6; 1; 0.1; 0.2; 1; 1];

% Вторая ситуация
vaccinated = 0.1;
infected = 0.5;
isolated = 0.01;

X_start_2 = [(1 - vaccinated - infected - isolated) * population; vaccinated * population; infected * population; isolated * population; vaccinated_speed];
params_for_second_situation = [0.1; 2 * 1e-7; 1e-7; 1e-2; 0.8; 0.2; 1; 1];
%----End of parameters----

f1 = @(t, x) SIR_model(t, x, params_for_first_situation);
f2 = @(t, x) SIR_model(t, x, params_for_second_situation);

figure;

subplot(1, 2, 1);
hold on;
[t, X] = ode45(f1, 0:dt:T, X_start_1);

drawing_S_1 = plot(t, X(:, 1), 'b-');
drawing_S_2 = plot(t, X(:, 2), 'g-');
drawing_I = plot(t, X(:, 3), 'r-');
drawing_R = plot(t, X(:, 4), 'k-');
drawing_population = plot(t, X(:, 1) + X(:, 2) + X(:, 3) + X(:, 4), 'y-');

legend([drawing_S_1, drawing_S_2, drawing_I, drawing_R, drawing_population], {'непривитых (S_1)', 'привитых (S_2)', 'инфицированных (I)', 'изолированных (R)', 'население'});
xlim([0 T]);
ylim([0 population]);
xlabel('time');
ylabel('number of people');
title('Население вымирает');
hold off;

subplot(1, 2, 2);
hold on;
[t, X] = ode45(f2, 0:dt:T, X_start_2);

drawing_S_1 = plot(t, X(:, 1), 'b-');
drawing_S_2 = plot(t, X(:, 2), 'g-');
drawing_I = plot(t, X(:, 3), 'r-');
drawing_R = plot(t, X(:, 4), 'k-');
drawing_population = plot(t, X(:, 1) + X(:, 2) + X(:, 3) + X(:, 4), 'y-');

legend([drawing_S_1, drawing_S_2, drawing_I, drawing_R, drawing_population], {'непривитых (S_1)', 'привитых (S_2)', 'инфицированных (I)', 'изолированных (R)', 'население'});
xlim([0 T]);
ylim([0 population]);
xlabel('time');
ylabel('number of people');
title('Ситуация стабилизировалась');
hold off;


%% Задание 3.6
clear

%----Parameters-----------
alpha = 1.1;
beta = 2.5;

dt = 0.01;
t_1 = 5;
t_2 = 20;

v_0 = [-5; -3];
x_0 = [-2; -1];
r = 10;
%----End of parameters----

figure;
axis([-5 5 -5 5]);
hold on;
plot([-4 -4], [-4  4], 'k-', 'LineWidth', 2);
plot([ 4  4], [-4  4], 'k-', 'LineWidth', 2);
plot([-4  4], [-4 -4], 'k-', 'LineWidth', 2);
plot([-4  4], [ 4  4], 'k-', 'LineWidth', 2);

x = x_0;
v = v_0;
ball = plot(x(1), x(2), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', r);

for t = 0:dt:t_1
    if x(1) < -4
        v(1) = -v(1);
        v = v / alpha;
        x(1) = -4;
    elseif x(1) > 4
        v(1) = -v(1);
        v = v / alpha;
        x(1) = 4;
    end
    if x(2) < -4
        v(2) = -v(2);
        v = v / alpha;
        x(2) = -4;
    elseif x(2) > 4
        v(2) = -v(2);
        v = v / alpha;
        x(2) = 4;
    end
    
    x = x + v*dt;
    set(ball, 'XData', x(1), 'YData', x(2));
    drawnow;
end

plot([0  0], [-4 4], 'k-', 'LineWidth', 2);
is_ball_in_left_part = (x(1) <= 0);

for t = t_1:dt:t_2
    if x(1) < -4
        v(1) = -v(1);
        v = v / alpha;
        x(1) = -4;
    elseif x(1) > 4
        v(1) = -v(1);
        v = v / alpha;
        x(1) = 4;
    end
    if x(2) < -4
        v(2) = -v(2);
        v = v / alpha;
        x(2) = -4;
    elseif x(2) > 4
        v(2) = -v(2);
        v = v / alpha;
        x(2) = 4;
    end
    
    if is_ball_in_left_part && x(1) > 0
        v(1) = -v(1);
        v = v / beta;
        x(1) = 0;
    elseif ~is_ball_in_left_part && x(1) < 0
        v(1) = -v(1);
        v = v / beta;
        x(1) = 0;
    end
    
    x = x + v*dt;
    set(ball, 'XData', x(1), 'YData', x(2));
    drawnow;
end

hold off;


%% Задание 3.7
clear

%----Parameters-----------
G = 6.7;
m_1 = 10; %10 %1.5
m_2 = 10; %10 %3
K = 5; %5 %4
r = 10;

x_start = [1, 0, -2, 0, 0, 2, 0, -1]; %[1, 0, 0, 1, 0, 4, 2, 0]; %[1, 0, -1, 0, -3, 3, 3, -3]; %[1, 0, -2, 0, 0, 2, 0, -1];
T = 100;
dt = 0.01;
%----End of parameters----

f = @(t, x) TwoBodyProblem(t, x, m_1, m_2, G);
time_period = 0:dt:T;
[t, x] = ode45(f, time_period, x_start);

figure;
axis([-K K -K K]);
title('Two-body problem');
xlabel('x');
ylabel('y');

hold on;
first_body = plot(x(1, 1), x(1,2), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', r);
second_body = plot(x(1, 3), x(1, 4), 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'MarkerSize', r);
mass_center = plot((x(1, 1) * m_1 + x(1, 3) * m_2) / (m_1 + m_2), (x(1, 2) * m_1 + x(1, 4) * m_2) / (m_1 + m_2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', r / 2);
legend('body #1', 'body #2', 'center of mass');

for n = 1:length(time_period)
    set(first_body, 'XData', x(n, 1), 'YData', x(n, 2));
    set(second_body, 'XData', x(n, 3), 'YData', x(n, 4));
    set(mass_center, 'XData', (x(n, 1) * m_1 + x(n, 3) * m_2) / (m_1 + m_2), 'YData', (x(n, 2) * m_1 + x(n, 4) * m_2) / (m_1 + m_2));
    drawnow;
end
hold off;


%% 3.8
clear

figure;

% Седло

f = @(t, y) [y(2); -3*y(2) + 4*y(1)];

n = 15;
a = -100;
b = 100;
time_period = 0:0.1:100;

x = linspace(a, b, n);
y = linspace(a, b, n);
[X, Y] = meshgrid(x, y);
U = zeros(size(X));
V = zeros(size(Y));
for i = 1:n^2
    temp = f(time_period, [X(i); Y(i)]);
    U(i) = temp(1);
    V(i) = temp(2);
end

subplot(2, 3, 1);
hold on;
quiver(X, Y, U, V);
hold on;
y_0 = [a/5 b];
[t, y] = ode45(f, time_period, y_0);
plot(y(:,1), y(:,2), 'b', 'LineWidth', 2);

plot(x, x, 'k-');
plot(x, -4*x, 'k-');

xlabel('x');
ylabel('y');
xlim([a b]);
ylim([a b]);
legend('Field', 'One solution');
title('Седло');
hold off;

% Фокус

f = @(t, y) [y(1) - 2*y(2); 2*y(2) + 2*y(1)];

n = 15;
a = -100;
b = 100;
time_period = 0:0.1:100;

x = linspace(a, b, n);
y = linspace(a, b, n);
[X, Y] = meshgrid(x, y);
U = zeros(size(X));
V = zeros(size(Y));
for i = 1:n^2
    temp = f(time_period, [X(i); Y(i)]);
    U(i) = temp(1);
    V(i) = temp(2);
end

subplot(2, 3, 2);
hold on;
quiver(X, Y, U, V);
hold on;
time_period = [0 100];
y_0 = [0 1];
[t, y] = ode45(f, time_period, y_0);
plot(y(:,1), y(:,2), 'b', 'LineWidth', 2);

xlabel('x');
ylabel('y');
xlim([a b]);
ylim([a b]);
legend('Field', 'One solution');
title('Фокус');
hold off;

% Узел

f = @(t, y) [y(2); -3*y(1) - 4*y(2)];

n = 20;
a = -20;
b = 20;
time_period = 0:0.1:100;

x = linspace(a, b, n);
y = linspace(a, b, n);
[X, Y] = meshgrid(x, y);
U = zeros(size(X));
V = zeros(size(Y));
for i = 1:n^2
    temp = f(time_period, [X(i); Y(i)]);
    U(i) = temp(1);
    V(i) = temp(2);
end

subplot(2, 3, 3);
hold on;
quiver(X, Y, U, V);
hold on;
y_0 = [a/2 a];
[t, y] = ode45(f, time_period, y_0);
plot(y(:,1), y(:,2), 'b', 'LineWidth', 2);

plot(x, -x, 'k-');
plot(x, -3*x, 'k-');

xlabel('x');
ylabel('y');
xlim([a b]);
ylim([a b]);
legend('Field', 'One solution');
title('Узел');
hold off;

% Центр

f = @(t, y) [y(2); -y(1)];

n = 15;
a = -100;
b = 100;
time_period = 0:0.1:100;

x = linspace(a, b, n);
y = linspace(a, b, n);
[X, Y] = meshgrid(x, y);
U = zeros(size(X));
V = zeros(size(Y));
for i = 1:n^2
    temp = f(time_period, [X(i); Y(i)]);
    U(i) = temp(1);
    V(i) = temp(2);
end

subplot(2, 3, 4);
hold on;
quiver(X, Y, U, V);
hold on;
y_0 = [a/2 b/2];
[t, y] = ode45(f, time_period, y_0);
plot(y(:,1), y(:,2), 'b', 'LineWidth', 2);

xlabel('x');
ylabel('y');
xlim([a b]);
ylim([a b]);
legend('Field', 'One solution');
title('Центр');
hold off;


% Диктритический узел

f = @(t, y) [y(1); y(2)];

n = 15;
a = -100;
b = 100;
time_period = 0:0.1:100;

x = linspace(a, b, n);
y = linspace(a, b, n);
[X, Y] = meshgrid(x, y);
U = zeros(size(X));
V = zeros(size(Y));
for i = 1:n^2
    temp = f(time_period, [X(i); Y(i)]);
    U(i) = temp(1);
    V(i) = temp(2);
end

subplot(2, 3, 5);
hold on;
quiver(X, Y, U, V);
hold on;
y_0 = [1 1];
[t, y] = ode45(f, time_period, y_0);
plot(y(:,1), y(:,2), 'b', 'LineWidth', 2);

xlabel('x');
ylabel('y');
xlim([a b]);
ylim([a b]);
legend('Field', 'One solution');
title('Диктритический узел');
hold off;


%% Задание 3.11
clear

%----Parameters-----------
m = 1;
L = 100;
g = 9.8;

dt = 0.01;
T = 100;
theta_0 = 2;
derivative_theta_0 = 0;
%----End of parameters----

f = @(t, y) [y(2); - g / L * sin(y(1))];

time_period = 0:dt:T;
y_0 = [theta_0 derivative_theta_0];
[t, thetas] = ode45(f, time_period, y_0);

figure;
K = L * 5 / 4;
axis([-K K -K K]);
hold on;

plot([-K -K], [-K  K], 'k-', 'LineWidth', 2);
plot([ K  K], [-K  K], 'k-', 'LineWidth', 2);
plot([-K  K], [-K -K], 'k-', 'LineWidth', 2);
plot([-K  K], [ K  K], 'k-', 'LineWidth', 2);

x = [L * sin(theta_0), - L * cos(theta_0)];
ball = plot(x(1), x(2), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 10);
line = plot([0 x(1)], [0 x(2)], 'k-', 'LineWidth', 2);
plot(0, 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5);

N = size(thetas);
N = N(1);
for n = 1:N
    theta = thetas(n, 1);
    x = [L * sin(theta), - L * cos(theta)];
    
    delete(line);
    line = plot([0 x(1)], [0 x(2)], 'k-', 'LineWidth', 2);
    set(ball, 'XData', x(1), 'YData', x(2));
    drawnow;
    
    kinetic_energy = m * (L * thetas(n, 2))^2 / 2;
    potential_energy = m * g * (- L * cos(theta) + L);
    title(sprintf('Движение маятника\n\nКинетическая энергия: %d Дж\nПотенциальная энергия: %d Дж\nСуммарная энергия: %d Дж', round(kinetic_energy), round(potential_energy), round(kinetic_energy + potential_energy)));
end

hold off;




