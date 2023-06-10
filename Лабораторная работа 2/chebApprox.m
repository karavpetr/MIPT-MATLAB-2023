function chebApprox(f, n)
    T = getFunc(n);
    x = linspace(-1, 1, 1000);
    y = f(x);
    
    Video1(1:n) = struct("cdata", [], "colormap", []);
    
    fig = figure;
    
    for k = 1:n
        hold on;
        
        % Calculate the k-th partial sum
        S = zeros(size(x));
        for j = 1:k
            sz = length(x);
            coeff = 2/sz * sum(f(cos(pi*(0:sz-1 + 0.5)/sz)) .* (cos(pi*j*(0:sz-1 + 0.5)/sz)));
            S = S + coeff * T{j + 1}(x);
        end
        
        drawing = plot(x, S, 'g-', x, y, 'r-');
        legend('Partial sum', 'f');
        abs_diff = max(abs(S - f(x)));
        title(sprintf('Convergence: (frame #%d/%d)\nMax absolute difference = %.4f', k, n, abs_diff));
        
        pause(0.5);
        hold off;
        
        Video1(k) = getframe(gcf);
        delete(drawing);
    end
    
    vidObj = VideoWriter('./chebApprox.avi');
    vidObj.FrameRate = 5;
    vidObj.Quality = 100;
    open(vidObj);
    writeVideo(vidObj, Video1);
    close(vidObj);
end

function T = getFunc(n)
    T{1} = @(x) ones(size(x));
    T{2} = @(x) x;

    for k = 3:n+1
        T{k} = @(x) 2*x.*T{k-1}(x) - T{k-2}(x);
    end
end
