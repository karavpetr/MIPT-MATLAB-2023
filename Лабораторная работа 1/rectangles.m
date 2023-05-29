function I = rectangles(f, a, b, n)
    h = (b-a)/n;
    x = linspace(a, b, n);
    I = h * sum(f(x(1:end-1)));
end

