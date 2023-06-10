function compareInterp(x,xx,f)
    yy = f(xx);
    y = f(x);
    
    % Interpolation using different methods
    yy_nearest = interp1(x, y, xx, 'nearest');
    yy_linear = interp1(x, y, xx, 'linear');
    yy_spline = interp1(x, y, xx, 'spline');
    yy_cubic = interp1(x, y, xx, 'PCHIP');
    
    % Plot the original function and interpolated functions
    figure;
    hold on;
    plot(xx, yy, 'b');
    plot(xx, yy_nearest, 'r--');
    plot(xx, yy_linear, 'g--');
    plot(xx, yy_spline, 'm--');
    plot(xx, yy_cubic, 'k--');
    legend('f', 'nearest', 'linear', 'spline', 'cubic');
    xlabel('x');
    ylabel('y = f(x)');
    title('Interpolation Methods');
    hold off;

end