function dx = TwoBodyProblem(t, x, m_1, m_2, G)
    dx = zeros([8 1]); % [x1 y1 x2 y2 x1' y1' x2' y2']
    norma = norm([x(1) - x(3), x(2) - x(4)]) ^ 3;
    
    dx(1) = x(5);
    dx(2) = x(6);
    dx(3) = x(7);
    dx(4) = x(8);
    
    dx(5) = G * m_2 * (x(3) - x(1)) / norma;
    dx(6) = G * m_2 * (x(4) - x(2)) / norma;
    dx(7) = G * m_1 * (x(1) - x(3)) / norma;
    dx(8) = G * m_1 * (x(2) - x(4)) / norma;
end