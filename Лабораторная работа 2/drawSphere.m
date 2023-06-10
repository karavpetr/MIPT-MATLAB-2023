function drawSphere(alpha, color, value, delta)
    if alpha == Inf
        f = @(x, y, z) max(abs(x), max(abs(y),abs(z)));
    else
        f = @(x, y, z) (abs(x) .^ alpha + abs(y) .^ alpha + abs(z) .^ alpha) .^ (1 / alpha);
    end
    
    a = -5;
    b = 5;
    h = delta;
    [xx3Mat, yy3Mat, zz3Mat] = meshgrid(a:h:b, a:h:b, a:h:b);
    rr3Mat = f(xx3Mat, yy3Mat, zz3Mat);
    
    figure;
    K = 4/3;
    xlim([-K K]);
    ylim([-K K]);
    zlim([-K K]);
    hold on;
    grid on;
    %patch('EdgeColor', color);
    surface = patch(isosurface(xx3Mat, yy3Mat, zz3Mat, rr3Mat, value, color));
    set(surface,'FaceColor',color,'EdgeColor','none');
    camlight
    view(3)
    hold off;
end