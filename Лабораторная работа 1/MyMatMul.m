function CMat = MyMatMul(AMat, BMat)
    [n, ~] = size(AMat);
    [~, k] = size(BMat);
    CMat = zeros(n, k);

    for i = 1:n
        for j = 1:k
            CMat(i, j) = AMat(i, :) * BMat(:, j);
        end
    end

end
