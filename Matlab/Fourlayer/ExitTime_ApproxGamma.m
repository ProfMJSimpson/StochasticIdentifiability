function f = ExitTime_ApproxGamma(l,x0,p)

    % Inputs:
    %   l       = [l(1),l(2),l(3),l(4)]
    %   x0      = initial location
    %   p       = [p(1),p(2),p(3),p(4)] = [D1,D2,D3,D4]
    %
    % Outputs:
    %   f  = @(t) ....
    
    %% Get raw moments
    [M1,M2] = MomentsRaw4(l(1),l(2),l(3),l(4));
    
    %% Get central moments
    [M,V] = MomentsCentral4(M1,M2);

    %m = M(x0,p(1),p(2),p(3),p(4));
    %v = V(x0,p(1),p(2),p(3),p(4));
    
    %% Get Gamma parameters
    %a = m^2 / v;
    %b = v / m;
    [a,b] = GammaParams4(M,V);
    
    %% Return interpolating function
    f  = @(t) gampdf(t,a(x0,p(1),p(2),p(3),p(4)),b(x0,p(1),p(2),p(3),p(4)));
    
end
