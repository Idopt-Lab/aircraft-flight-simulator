function coeff = QuadAeroLookup(x, ~, ~)
    alt = max(-x(3), 0);
    [~, ~, ~, rho] = atmosisa(alt);
    
    u = x(4); v = x(5); w = x(6);
    V = sqrt(u^2 + v^2 + w^2);
    
    CD0 = 0.8;
    
    CL = 0;
    CY = 0;
    Cl = 0;
    Cm = 0;
    Cn = 0;
    CD = CD0;
    
    coeff = struct('CL', CL, 'CD', CD, 'CY', CY, 'Cl', Cl, 'Cm', Cm, 'Cn', Cn);
end