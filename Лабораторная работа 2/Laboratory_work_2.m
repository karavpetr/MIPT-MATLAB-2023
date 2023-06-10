% cd /home/karavpetr/Desktop/MATLAB

%% Задание 1
clear

%----Parameters-----------
f = @(x) x .* sin(1 ./ x.^2);
a = 0.2;
b = 1.2;

x = linspace(a, b, 11);
xx = linspace(a, b, 101);
%----End of parameters----

compareInterp(x,xx,f);

%% Задание 3
clear
format compact

%----Parameters-----------
f = @(x) x .^ 2;
fn = @(n,x) n * atan(x .^ 2 ./ n);

a = 1;
b = 2;
n = 10;  % Number of frames
%----End of parameters----

% possible parameters: pointwise, uniform, rms
convergenceFunc(fn, f, a, b, n, 'rms');


%% Задание 4
clear

%----Parameters-----------
f = @(x) (x<-0.2).*(-1) + (x>=-0.2).*(x<=0.2) .*x ./ 0.2  + (x>0.2).*1;
n = 20;
%----End of parameters----

chebApprox(f, n);


%% Задание 2.5
clear

%----Parameters-----------
x = linspace(4, 10, 1000);
f = @(x) sin(2*x) + 0.5*sin(6*x);
%f = @(x) fillmissing(sqrt(abs(x)) .* sin(1./x.^2) .* (x ~= 0), 'constant', 0);
%----End of parameters----

y = f(x);

[~, indexes_local_minimums] = findpeaks(-y);
[~, indexes_local_maximums] = findpeaks(y);
[~, indexes_global_maximums] = max(y(indexes_local_maximums));

index_of_global_maximum = indexes_local_maximums(indexes_global_maximums(1));
[~, indexes_of_nearest_local_minimums] = min(abs(x(index_of_global_maximum) - x(indexes_local_minimums)));
index_of_nearest_local_minimum = indexes_local_minimums(indexes_of_nearest_local_minimums(1));

figure;
plot(x, y, 'LineWidth', 1);
hold on;
plot(x(indexes_local_minimums), y(indexes_local_minimums), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
plot(x(index_of_global_maximum), y(index_of_global_maximum), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

title('Function with local minimums and global maximum');
xlabel('x');
ylabel('f(x)');

comet_x = linspace(x(index_of_global_maximum), x(index_of_nearest_local_minimum), 100);
comet_y = f(comet_x);

comet(comet_x, comet_y, 0.1);


%% Задание 6
clear

%----Parameters-----------
r = 5.11;
f = @(x) r*x.*(1-sqrt(x));
x0 = 0.01;
number_of_frames = 10;
%----End of parameters----

x = linspace(0,1,100);

Video1(1:number_of_frames) = struct("cdata", [], "colormap", []);
fig = figure;

xlabel('x');
ylabel('y = f(x)');

x_prev = x0;
for i=1:number_of_frames
    hold on;
    x_next = f(x_prev);
    plot([x_prev, x_next], [x_next, x_next], 'r-', 'linewidth', 1);
    plot([x_prev, x_prev], [x_prev, x_next], 'r-', 'linewidth', 1);
    x_prev = x_next;
    
    drawing_1 = plot(x, f(x), 'b-', 'linewidth', 2);
    drawing_2 = plot(x, x, 'k--', 'linewidth', 1);
    legend([drawing_1, drawing_2], {'f(x)', 'y = x'});
    title(sprintf('Ход приближений: (frame #%d/%d)', i, number_of_frames));
    
    hold off;
    pause(0.3);
    Video1(i) = getframe(gcf);
end
pause(1);
close(fig);

vidObj = VideoWriter('./IterativAprox.avi');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj, Video1);
close(vidObj);


%% Задание 7
clear

%----Parameters-----------
x_min = -1;
x_max = 1;
y_min = -2;
y_max = 2;

n = 10; %number of molecues

delta_time = 0.01;
time_duration = 10;
%----End of parameters----

% start features of molecules
x = x_min + (x_max - x_min) * rand(n,1);
y = y_min + (y_max - y_min) * rand(n,1);
vx = randn(n,1);
vy = randn(n,1);

fig = figure;
for t = 0:delta_time:time_duration
    hold on;
    drawing = plot(x, y, 'ro');
    
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    
    % Updating positions
    x = x + vx * delta_time;
    y = y + vy * delta_time;
    
    % Check for collisions with the walls
    x_bad_points = find((x < x_min) | (x > x_max));
    vx(x_bad_points) = -vx(x_bad_points);
    y_bad_points = find((y < y_min) | (y > y_max));
    vy(y_bad_points) = -vy(y_bad_points);

    drawnow;
    hold off;
    
    pause(delta_time);
    delete(drawing);
end
close(fig);


%% Задание 8
clear

% Задание 8.1

data = readtable('incomes.csv', 'Delimiter', ';');

income_data = table2array(data(:, 2:end));
income_data = strrep(income_data, ',', '.');
income_data = str2double(income_data);

names = data.Properties.VariableNames;
dates = datetime(data.("Date"), 'InputFormat', 'dd.MM.yyyy');
total_income_data = sum(income_data, 2);

plot(dates, total_income_data);
title('Total income data in all regions');
xlabel('Date');
ylabel('Total income (in billions of rubles)');

table_for_file = table(dates, total_income_data);
writetable(table_for_file, 'inc.txt');

%% Задание 8.2

years = year(dates);
regions = ["Belgorodskaya_oblast", "Bryanskaya_oblast", "Vladimirskaya_oblast", "Lipetskaya_oblast", "Respublica_Dagestan", "Kemerovskaya_oblast"];
regions_names = strrep(regions, '_', ' ');

unique_years = unique(years);
total_incomes = zeros([length(unique_years) length(regions)]);

for n = 1:length(unique_years)
    for i = 1:length(regions)
        regionData = data.(regions(i));
        regionData = strrep(regionData, ',', '.');
        regionData = str2double(regionData);
        year_indeces = find(years == years(n));
        total_incomes(n, i) = sum(regionData(year_indeces));
    end
end

figure;
hold on;
b = bar(unique_years, total_incomes);
legend(regions_names);

xlabel('year');
ylabel('Total income');
title('Total incomes of different regions');
hold off;

%% Задание 8.3

year2020_indices = find(years == 2020);

shares = zeros(1, length(regions));

for i = 1:length(regions)
    regionData = data.(regions(i));
    regionData = strrep(regionData, ',', '.');
    regionData = str2double(regionData);
    shares(i) = sum(regionData(year2020_indices));
end

figure;
pie(shares, regions_names);
title('Total incomes of regions in 2020');


%% Задание 9
clear

%----Parameters-----------
a = -10;
b = 10;
h = 0.5;

f = @(X, Y, alpha) Y.*sin(X) * 9/10 - alpha.*X.*cos(alpha.*Y);
number_of_frames = 25;
alphas = 1:0.05:10;
%----End of parameters----


Video1(1:number_of_frames) = struct("cdata", [], "colormap", []);
fig = figure;
[X, Y] = meshgrid(a:h:b);

for i=1:number_of_frames
    Z = f(X, Y, alphas(i));
    surface = surf(X,Y,Z,'FaceAlpha',0.5);
    
    local_maximums = find(imregionalmax(Z));
    local_minimums = find(imregionalmin(Z));
    hold on;
    drawing_1 = plot3(X(local_maximums),Y(local_maximums),Z(local_maximums),'r.','MarkerSize',24);
    drawing_2 = plot3(X(local_minimums),Y(local_minimums),Z(local_minimums),'b.','MarkerSize',24);
    
    hold off;
    pause(0.5);
    Video1(i) = getframe(gcf);
    delete(surface);
    delete(drawing_1);
    delete(drawing_2);
end
pause(1);
close(fig);

vidObj = VideoWriter('./SurfaceEvolution.avi');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj, Video1);
close(vidObj);


%% Задание 10

%----Parameters-----------
f = @(x, y) (x.^2 - 11).^2 + y.^3;
%----End of parameters----

x = linspace(-5, 5, 100);
y = linspace(-5, 5, 100);
contour_level = 10;

[X, Y] = ndgrid(x, y);
Z = f(X, Y);

figure;
hold on;
contour(x, y, Z, contour_level, 'LineWidth', 2);

xlabel('x');
ylabel('y');
title('Contour of f(x, y)');
hold off;


%% Задание 11

% polar coordinates
theta = linspace(0, 2*pi, 100);
r = linspace(0, 1, 100);
[theta, r] = meshgrid(theta, r);

x = r.*cos(theta);
y = r.*sin(theta);
z = 5*r;

figure;
surf(x, y, z);
hold on;
surf(x, y, -z);
hold off;


%% Задание 12
clear

alpha = 3.0;
color = 'yellow';
value = 0.8;
delta = 0.1;

drawSphere(alpha, color, value, delta);


%% Задание 13
clear;

alphas = [10, 3, 1.2];
colors = ["red", "yellow", "green"];
edges = ["none", "none", "none"];
delta = 0.2;

drawManySpheres(alphas, colors, edges, delta);

