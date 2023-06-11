%% Задание 4.1
clear

%----Parameters-----------
A = [1, 2; 3*i, 1];
B = [3, 4; 2*i, i];
C = [5, 6; 2, 9+i];
%----End of parameters----

[x1, x2, x3, x4] = biquadrate_solution(A, B, C);

disp(x4);

% Проверка корней на корректность (норма должна быть маленькой)
disp(norm((x1 .^ 4) .* A + (x1 .^ 2) .* B + C, 'fro'));
disp(norm((x2 .^ 4) .* A + (x2 .^ 2) .* B + C, 'fro'));
disp(norm((x3 .^ 4) .* A + (x3 .^ 2) .* B + C, 'fro'));
disp(norm((x4 .^ 4) .* A + (x4 .^ 2) .* B + C, 'fro'));


%% Задание 4.2
clear

A = rand(5, 3);
disp('A = ');
disp(A);

disp('With C++')
[Q1, R1] = qr_c(A);
disp('Q1 = ');
disp(Q1);
disp('R1 = ');
disp(R1);
disp('(C++) ||Q*R - A|| = ');
disp(norm(Q1*R1 - A, 'fro'));

disp('With MatLab')
[Q2, R2] = qr_m(A);
disp("Q2 =");
disp(Q2);
disp("R2 =");
disp(R2);
disp('(MatLab) ||Q*R - A|| = ');
disp(norm(Q2*R2 - A, 'fro'));

disp('qr');
[Q3, R3] = qr(A, 'matrix');
disp("Q3 =");
disp(Q3);
disp("R3 =");
disp(R3);
disp('(qr) ||Q*R - A|| = ');
disp(norm(Q3*R3 - A, 'fro'));

disp('(C++) ||Q1*R1 - A|| = ');
disp(norm(Q1*R1 - A, 'fro'));
disp('(MatLab) ||Q1*R1 - A|| = ');
disp(norm(Q2*R2 - A, 'fro'));
disp('(qr) ||Q1*R1 - A|| = ');
disp(norm(Q3*R3 - A, 'fro'));

%% Задание 4.3

dims = [50 150 250 400 550 700];

frob_qr = zeros(length(dims),1);
frob_qr_c = zeros(length(dims),1);
frob_qr_m = zeros(length(dims),1);

for i = 1:length(dims)
    A = rand(dims(i));

    % qr
    [Q,R] = qr(A);
    frob_qr(i) = norm(Q*R-A, 'fro');

    % qr_c
    [Q_c,R_c] = qr_c(A);
    frob_qr_c(i) = norm(Q_c*R_c-A,'fro');

    % qr_m
    [Q_m,R_m] = qr_m(A);
    frob_qr_m(i) = norm(Q_m*R_m-A, 'fro');
end

% Plot results
plot(dims, frob_qr, '-o', dims, frob_qr_c, '-*', dims, frob_qr_m, '-x');
legend('qr', 'qr\_c', 'qr\_m');
xlabel('Matrix dimension (n*n size)');
ylabel('Value of norms');
title('Frob norm');


%% Задание 4.4

dims = [50 150 250 400 550 700];

t_qr = zeros(length(dims),1);
t_qr_c = zeros(length(dims),1);
t_qr_m = zeros(length(dims),1);

for i = 1:length(dims)
    A = rand(dims(i));

    % qr time
    tic
    [Q,R] = qr(A);
    t_qr(i) = toc;

    % qr_c time
    tic
    [Q_c,R_c] = qr_c(A);
    t_qr_c(i) = toc;

    % qr_m time
    tic
    [Q_m,R_m] = qr_m(A);
    t_qr_m(i) = toc;
end

% Plot results
plot(dims, t_qr, '-o', dims, t_qr_c, '-*', dims, t_qr_m, '-x');
legend('qr', 'qr\_c', 'qr\_m');
xlabel('Matrix dimension (n*n size)');
ylabel('Time (s)');
title('Working time of functions');

