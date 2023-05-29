function I = simpson(f, a, b, n)
    if mod(n, 2) == 1
        error('n must be even for the Simpson method');
    end
    
    h = (b - a) / n;
    x = linspace(a, b, n);
    I = (h/3) * (f(x(1)) + 2*sum(f(x(3:2:end-2))) + 4*sum(f(x(2:2:end-1))) + f(x(end)));
end
