function convergenceFunc(fn, f, a, b, n, convType)
    x = linspace(a, b, 1000);
    y = f(x);
    
    Video1(1:n) = struct("cdata", [], "colormap", []);
    
    fig = figure;
    for i = 1:n
        hold on;
        f_n = fn(i, x);
        drawing = plot(x, f_n, 'b-', x, y, 'r-');
        legend('f_n(x)', 'f (x)');
        
        if convType == "pointwise"
            title(sprintf('Pointwise convergence: (frame #%d/%d)', i, n))
        elseif convType == "uniform"
            abs_diff = max(abs(f_n - y));
            title(sprintf('Uniform convergence: (frame #%d/%d)\nMax absolute difference = %.4f', i, n, abs_diff));
        elseif convType == "rms"
            rms_diff = sqrt(mean((f_n - y).^2));
            title(sprintf('Uniform convergence: (frame #%d/%d)\nRMS difference = %.4f', i, n, rms_diff));
        else
            error('Invalid convergence type');
        end
        
        pause(0.5);
        hold off;

        Video1(i) = getframe(gcf);
        delete(drawing);
    end
    close(fig);
    
    vidObj = VideoWriter('./Convergence.avi');
    vidObj.FrameRate = 5;
    vidObj.Quality = 100;
    open(vidObj);
    writeVideo(vidObj, Video1);
    close(vidObj);
    
end