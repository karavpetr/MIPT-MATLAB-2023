function [Q, R] = qr_m(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    
    for j = 1:n
        v = A(:,j);
        for i = 1:j-1
            v = v - dot(Q(:,i), A(:,j)) * Q(:,i);
        end
        Q(:,j) = v / norm(v);
    end
    
    R = Q' * A;
end
