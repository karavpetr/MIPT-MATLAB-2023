function drawManySpheres(alphas, colors, edges, delta)
    a = -5;
    b = 5;
    h = delta;
    value = 1;
    [xx3Mat, yy3Mat, zz3Mat] = meshgrid(a:h:b, a:h:b, a:h:b);
    
    figure;
    K = 4/3;
    xlim([-K K]);
    ylim([-K K]);
    zlim([-K K]);
    hold on;
    grid on;
    
    for n = 1:length(alphas)
        alpha_value = alphas(n);
        if alpha_value == Inf
            f = @(x, y, z) max(abs(x), max(abs(y),abs(z)));
        else
            f = @(x, y, z) (abs(x) .^ alpha_value + abs(y) .^ alpha_value + abs(z) .^ alpha_value) .^ (1 / alpha_value);
        end
        
        rr3Mat = f(xx3Mat, yy3Mat, zz3Mat);
        surface = patch(isosurface(xx3Mat, yy3Mat, zz3Mat, rr3Mat, value));
        set(surface,'FaceColor',colors(n),'EdgeColor',edges(n));
    end
    
    alpha(0.5);
    camlight
    view(3)
    hold off;
end