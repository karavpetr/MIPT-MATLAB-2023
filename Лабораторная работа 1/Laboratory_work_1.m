%% Лабораторная работа №1

%% Задание 1
clear

%----Parameters-----------
a = -2;
b = 1;
n = 1000000;
%----End of parameters----

xVec = linspace(a, b, n);
f = @(x) cos(10 * x) + sqrt(abs(x));
fVec = f(xVec);
eps = 1e-6;
index_of_min = find(fVec <= min(fVec) + eps);
index_of_max = find(fVec >= max(fVec) - eps);
    
figure;
hold on;
plot(xVec, fVec, 'b-');
plot(xVec(index_of_min), f(xVec(index_of_min)), '*', 'color', 'red');
plot(xVec(index_of_max), f(xVec(index_of_max)), '*', 'color', 'red');
hold off;

%% Задание 2
clear

%----Parameters-----------
n = input('Input n: ');
%----End of parameters----

% 2.0
if round(n) - n == 0
    fprintf('Input is integer\n');
else
    error('Input in not integer');
end

% 2.1
xVec = 9:18:n;
disp(xVec);

% 2.2
aMat = n * ones(n, n);
bMat = reshape(repelem(1:n, n), n, n).';
resMat = aMat;
resMat(1:2:end, :) = bMat(1:2:end, :);
disp(resMat);

% 2.3
bMat = 1:(n*(n+1));
bMat = reshape(bMat, [n+1, n]).';
disp(bMat);

cVec = bMat(:);   
disp(cVec);

DMat = bMat(:, 2:2:end);
disp(DMat);

%% Задание 3
clear

%----Parameters-----------
n = 6;
%----End of parameters----

aMat = round(rand(n)*100 + 0.5);
disp(aMat);

fprintf('Max diag element: %d\n', max(diag(aMat)));

[dMat, cVec] = sort(aMat(:, 3));
eMat = aMat(cVec, :);
disp(cVec);
disp(eMat);

%% Задание 4
clear

%----Parameters-----------
n = 2;
m = 4;

xVec = randi(5, 1, n);
yVec = randi(5, 1, m);
%----End of parameters----

aMat = xVec.' * yVec(end:-1:1);
disp(xVec);
disp(yVec);
disp(aMat);

%% Задание 5
clear

%----Parameters-----------
n = 3;

AMat = rand(n, n);
bVec = rand(n, 1);
%----End of parameters----

if isprime(n)
    fprintf('Parameter n is prime\n');
else
    error('Parameter n is not prime');
end

if abs(det(AMat)) < 1e-6
    fprintf('Матрица A вырожденная\n');
else
    first_solution = inv(AMat) * bVec;
    fprintf('First solution:\n');
    disp(first_solution);
    fprintf('Is correct:\n');
    disp(pdist2((AMat * first_solution)', bVec') < 1e-9);
    
    helpingMat = rref([AMat bVec]);
    second_solution = helpingMat(:, end);
    fprintf('Second solution:\n');
    disp(second_solution);
    fprintf('Is correct:\n');
    disp(pdist2((AMat * second_solution)', bVec') < 1e-9);
end

%% Задание 6
clear

%----Parameters-----------
n = 10;
k = 10;

aMat = rand(n, k);
%----End of parameters----

disp(aMat);
bMat = dot(aMat.', aMat.');
cMat = reshape(repelem(reshape(bMat, n, 1), n), n, n);

resMat = cMat + cMat' - aMat * aMat.' * 2;
disp(resMat);

%% Задание 7
clear

%----Parameters-----------
n = 4;
%----End of parameters----

AMat = num2cell(dec2bin(0:2^n-1, n));
AMat = reshape(bin2dec(AMat), [2^n, n]);
disp(AMat);

AMat = AMat'; % Симметричных стобцов никогда нет, ищем симметричные строки

A2Mat = fliplr(AMat')';
disp(find(sum(abs(AMat - A2Mat)) == 0));


%% Задание 8
clear

%----Parameters-----------
n = 5;
m = 3;

aMat = rand(m, n);
%----End of parameters----

min_size = min(m, n);
max_size = max(m, n);
disp(aMat);

bMat = diag(ones(1, max_size));
cMat = diag(cat(1, diag(aMat), ones(max_size - min_size, 1)));

resMat = 5*ones(size(aMat)) - 5 * bMat(1:m, 1:n) + cMat(1:m, 1:n);
disp(resMat);

%% Задание 9
clear

%----Parameters-----------
N_values = 20:10:100;
%----End of parameters----

my_time_values = zeros(size(N_values));
ml_time_values = zeros(size(N_values));

for i = 1:size(N_values, 2)
    n = N_values(i);
    
    AMat = rand(n, n);
    BMat = rand(n, n);
 
    f = @() MLMatMul(AMat, BMat);
    ml_time_values(i) = timeit(f);

    f = @() MyMatMul(AMat, BMat);
    my_time_values(i) = timeit(f);
end
    
figure;
hold on;
plot(N_values, my_time_values, 'b-');
plot(N_values, ml_time_values, 'g-');
hold off;

%% Задание 10
clear

XMat = [NaN, 1, 2; NaN, 0, 6; 1, 5, NaN];
disp(nanMean(XMat));

%% Задание 11
clear

%----Parameters-----------
a = 10;
sigma = 2;
n = 1000000;
xVec = randn(1, n) * sigma + a;
%----End of parameters----

array_for_more = xVec >= a - sigma * 3;
array_for_less = xVec <= a + sigma * 3;
array_for_segment = array_for_more & array_for_less;

fprintf('Доля попавших в интервал значений: %d\n', sum(array_for_segment) / n);


%% Задание 12
clear

%----Parameters-----------
a = -50;
b = 50;
n = 1000;

f = @(x) sin(x) ./ x;
%----End of parameters----

x = linspace(a, b, n);
trapz_integral = trapz(x, f(x));
simpson_integral = simpson(f, a, b, n);
rectangles_integral = rectangles(f, a, b, n);

%% (Задание 12) Первообразная
figure;
hold on;
title('Antiderivative of sin(x)/x');

N = 100;
x = linspace(a, b, N);
trapz_antiderivative = zeros([N 1]);
simpson_antiderivative = zeros([N 1]);
rectangles_antiderivative = zeros([N 1]);

for i = 1:N
    x_values = linspace(a, x(i), n);
    trapz_antiderivative(i) = trapz(x_values, f(x_values));
    simpson_antiderivative(i) = simpson(f, a, x(i), n);
    rectangles_antiderivative(i) = rectangles(f, a, x(i), n);
end

plot(x, trapz_antiderivative, 'r-');
plot(x, simpson_antiderivative, 'b-');
plot(x, rectangles_antiderivative, 'k-');
legend(["using trapz function", "using simpson method", "using rectangles method"]);
hold off;

%% (Задание12) Внутренняя скорость сходимости
n_values = 100:10:10000;
h_values = (b - a) ./ n_values;
trapz_errors = zeros([length(n_values) 1]);
simpson_errors = zeros([length(n_values) 1]);
rectangles_errors = zeros([length(n_values) 1]);

for i = 1:length(h_values)
    x_values_first = linspace(a, b, n_values(i));
    x_values_second = linspace(a, b, 2*n_values(i));
    trapz_errors(i) = abs(trapz(x_values_first, f(x_values_first)) - trapz(x_values_second, f(x_values_second)));
    simpson_errors(i) = abs(simpson(f, a, b, n_values(i)) - simpson(f, a, b, 2*n_values(i)));
    rectangles_errors(i) = abs(rectangles(f, a, b, n_values(i)) - rectangles(f, a, b, 2*n_values(i)));
end
    
figure;
hold on;
title("Внутренняя скорость сходимости");
plot(h_values, trapz_errors, 'r-', h_values, simpson_errors, 'b-', h_values, rectangles_errors, 'k-');
legend(["trapz function", "simpson method", "rectangles method"]);
xlabel('h');
ylabel('absolute difference');
hold off;

%% (Задание 12) Время работы

n_values = 1000:10:10000;
trapz_times = zeros([length(n_values) 1]);
simpson_times = zeros([length(n_values) 1]);
rectangles_times = zeros([length(n_values) 1]);

for i = 1:length(n_values)
    x_values = linspace(a, b, n_values(i));

    tic
    trapz(x_values, f(x_values));
    trapz_times(i) = toc;

    tic
    simpson(f, a, b, n_values(i));
    simpson_times(i) = toc;

    tic
    rectangles(f, a, b, n_values(i));
    rectangles_times(i) = toc;
end

figure;
hold on;
title("Time of work");
plot(n_values, trapz_times, 'r-', n_values, simpson_times, 'b-', n_values, rectangles_times, 'k-');
legend(["trapz function", "simpson method", "rectangles method"]);
xlabel('number of points (n)');
ylabel('time');
hold off;


%% Задание 13
clear

%----Parameters-----------
F = @(x) sin(x) + x.^2;
f = @(x) cos(x) + 2 .* x;

x_0 = 0.1;
%----End of parameters----

eps = logspace(-16, -1, 200);
f_right = (F(x_0 + eps) - F(x_0)) ./ eps;
differences_right = abs(f_right - f(x_0));
f_center = (F(x_0 + eps) - F(x_0 - eps)) ./ eps ./ 2;
differences_center = abs(f_center - f(x_0));

loglog(eps, differences_right, 'r-', eps, differences_center, 'k-');
xlabel('epsilon (log scale)');
ylabel('difference (log scale)');
legend('Right derivative Difference', 'Central derivative Difference');
title('Difference between approx and actual derivatives');

