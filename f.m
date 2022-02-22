function scores = f(T, B)

    switch (B)
        case 1
            scores = ackleyfcn(T);
        case 2
            scores = ackleyn2fcn(T);
        case 3
            scores = ackleyn3fcn(T);
        case 4
            scores = ackleyn4fcn(T);
        case 5
            scores = adjimanfcn(T);
        case 6
            scores = alpinen1fcn(T);
        case 7
            scores = alpinen2fcn(T);
        case 8
            scores = bartelsconnfcn(T);
        case 9
            scores = bealefcn(T);
        case 10
            scores = birdfcn(T);
        case 11
            scores = bohachevskyn1fcn(T);
        case 12
            scores = bohachevskyn2fcn(T);
        case 13
            scores = boothfcn(T);
        case 14
            scores = brentfcn(T);
        case 15
            scores = brownfcn(T);
        case 16
            scores = bukinn6fcn(T);
        case 17
            scores = crossintrayfcn(T);
        case 18
            scores = deckkersaartsfcn(T);
        case 19
            scores = dropwavefcn(T);
        case 20
            scores = easomfcn(T);
        case 21
            scores = eggcratefcn(T);
        case 22
            scores = exponentialfcn(T);
        case 23
            scores = goldsteinpricefcn(T);
        case 24
            scores = gramacyleefcn(T);
        case 25
            scores = griewankfcn(T);
        case 26
            scores = happycatfcn(T);
        case 27
            scores = himmelblaufcn(T);
        case 28
            scores = holdertablefcn(T);
        case 29
            scores = keanefcn(T);
        case 30
            scores = leonfcn(T);
        case 31
            scores = levin13fcn(T);
        case 32
            scores = matyasfcn(T);
        case 33
            scores = mccormickfcn(T);
        case 34
            scores = periodicfcn(T);
        case 35
            scores = powellsumfcn(T);
        case 36
            scores = qingfcn(T);
        case 37
            scores = quarticfcn(T);
        case 38
            scores = rastriginfcn(T);
        case 39
            scores = ridgefcn(T);
        case 40
            scores = rosenbrockfcn(T);
        case 41
            scores = salomonfcn(T);
        case 42
            scores = schaffern1fcn(T);
        case 43
            scores = schaffern2fcn(T);
        case 44
            scores = schaffern3fcn(T);
        case 45
            scores = schaffern4fcn(T);
        case 46
            scores = schwefel220fcn(T);
        case 47
            scores = schwefel221fcn(T);
        case 48
            scores = schwefel222fcn(T);
        case 49
            scores = schwefel223fcn(T);
        case 50
            scores = schwefelfcn(T);
        case 51
            scores = shubert3fcn(T);
        case 52
            scores = shubert4fcn(T);
        case 53
            scores = shubertfcn(T);
        case 54
            scores = spherefcn(T);
        case 55
            scores = styblinskitankfcn(T);
        case 56
            scores = sumsquaresfcn(T);
        case 57
            scores = threehumpcamelfcn(T);
        case 58
            scores = wolfefcn(T);
        case 59
            scores = xinsheyangn1fcn(T);
        case 60
            scores = xinsheyangn2fcn(T);
        case 61
            scores = xinsheyangn3fcn(T);
        case 62
            scores = xinsheyangn4fcn(T);
        case 63
            scores = zakharovfcn(T);
    end

end

function scores = ackleyfcn(x)
    n = size(x, 2);
    ninverse = 1 / n;
    sum1 = sum(x.^2, 2);
    sum2 = sum(cos(2 * pi * x), 2);

    scores = 20 + exp(1) - (20 * exp(-0.2 * sqrt(ninverse * sum1))) - exp(ninverse * sum2);
end

function scores = ackleyn2fcn(x)
%     n = size(x, 2);
    %     assert(n == 2, 'Ackley N. 2 function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = -200 * exp(-0.02 * sqrt((X.^2) + (Y.^2)));
end

function scores = ackleyn3fcn(x)

    n = size(x, 2);
    assert(n == 2, 'Ackley N. 3 function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = -200 * exp(-0.2 * sqrt((X.^2) + (Y.^2))) + ...
        5 * exp(cos(3 * X) + sin(3 * Y));
end

function scores = ackleyn4fcn(x)
    [m, n] = size(x);

    scores = zeros(m, 1);

    for i = 1:m

        for j = 1:(n - 1)
            scores = scores + exp(-0.2) .* sqrt(x(i, j).^2 + x(i, j + 1).^2) ...
                + 3 * (cos(2 * x(i, j)) + sin(2 * x(i, j + 1)));
        end

    end

end

function scores = adjimanfcn(x)

    n = size(x, 2);
    assert(n == 2, 'Adjiman function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (cos(X) .* sin(Y)) - (X ./ ((Y.^2) + 1));
end

function scores = alpinen1fcn(x)
    scores = sum(abs(x .* sin(x) + 0.1 * x), 2);
end

function scores = alpinen2fcn(x)
    scores = prod(sqrt(x) .* sin(x), 2);
end

function scores = bartelsconnfcn(x)

    n = size(x, 2);
    assert(n == 2, 'Bartels Conn function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = abs((X.^2) + (Y.^2) + (X .* Y)) + abs(sin(X)) + abs(cos(Y));
end

function scores = bealefcn(x)
    n = size(x, 2);
    assert(n == 2, 'Beale''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (1.5 - X + (X .* Y)).^2 + ...
        (2.25 - X + (X .* (Y.^2))).^2 + ...
        (2.625 - X + (X .* (Y.^3))).^2;
end

function scores = birdfcn(x)

    n = size(x, 2);
    assert(n == 2, 'Bird function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = sin(X) .* exp((1 - cos(Y)).^2) + ...
        cos(Y) .* exp((1 - sin(X)).^2) + ...
        (X - Y).^2;
end

function scores = bohachevskyn1fcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Bohachevsky N. 1 function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (X.^2) + (2 * Y.^2) - (0.3 * cos(3 * pi * X)) - (0.4 * cos(4 * pi * Y)) + 0.7;
end

function scores = bohachevskyn2fcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Bohachevsky N. 2 function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (X.^2) + (2 * Y.^2) - (0.3 * cos(3 * pi * X)) .* (cos(4 * pi * Y)) + 0.3;
end

function scores = boothfcn(x)

    n = size(x, 2);
    assert(n == 2, 'Booth''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (X + (2 * Y) - 7).^2 + ((2 * X) + Y - 5).^2;
end

function scores = brentfcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Brent function is defined only on the 2-D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (X + 10).^2 + (Y + 10).^2 + exp(-X.^2 - Y.^2);
end

function scores = brownfcn(x)

    n = size(x, 2);
    scores = 0;

    x = x.^2;

    for i = 1:(n - 1)
        scores = scores + x(:, i).^(x(:, i + 1) + 1) + x(:, i + 1).^(x(:, i) + 1);
    end

end

function scores = bukinn6fcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Bukin N. 6 functions is only defined on a 2D space.')

    X = x(:, 1);
    X2 = X.^2;
    Y = x(:, 2);

    scores = 100 * sqrt(abs(Y - 0.01 * X2)) + 0.01 * abs(X + 10);
end

function scores = crossintrayfcn(x)

    n = size(x, 2);
    assert(n == 2, 'The Cross-in-tray function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    expcomponent = abs(100 - (sqrt(X.^2 + Y.^2) / pi));

    scores = -0.0001 * ((abs(sin(X) .* sin(Y) .* exp(expcomponent)) + 1).^0.1);
end

function scores = deckkersaartsfcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Deckkers-Aarts function is defined only on the 2-D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (100000 * X.^2) + Y.^2 +- (X.^2 + Y.^2).^2 + (10^ - 5) * (X.^2 + Y.^2).^4;
end

function scores = dropwavefcn(x)
    n = size(x, 2);
    assert(n == 2, 'Drop-Wave function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    numeratorcomp = 1 + cos(12 * sqrt(X.^2 + Y.^2));
    denumeratorcom = (0.5 * (X.^2 + Y.^2)) + 2;
    scores =- numeratorcomp ./ denumeratorcom;
end

function scores = easomfcn(x)

    n = size(x, 2);
    assert(n == 2, 'The Easom''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = -cos(X) .* cos(Y) .* exp(-(((X - pi).^2) + ((Y - pi).^2)));
end

function scores = eggcratefcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Egg Crate function is defined only on the 2-D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = X.^2 + Y.^2 + (25 * (sin(X).^2 + sin(Y).^2));
end

function scores = exponentialfcn(x)
    x2 = x.^2;

    scores = -exp(-0.5 * sum(x2, 2));
end

function scores = goldsteinpricefcn(x)
    n = size(x, 2);
    assert(n == 2, 'The Goldstein-Price function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (1 + ((X + Y + 1).^2) * (19 - (14 * X) + (3 * (X.^2)) - 14 * Y + (6 .* X .* Y) + (3 * (Y.^2)))) .* ...
        (30 + ((2 * X - 3 * Y).^2) .* (18 - 32 * X + 12 * (X.^2) + 48 * Y - (36 .* X .* Y) + (27 * (Y.^2))));
end

function scores = gramacyleefcn(x)
    n = size(x, 2);
    assert(n == 1, 'Gramacy & Lee function is only defined on a 1-D space.')

    scores = (sin(10 .* pi .* x) ./ (2 * x)) + ((x - 1).^4);
end

function scores = griewankfcn(x)

    n = size(x, 2);

    sumcomp = 0;
    prodcomp = 1;

    for i = 1:n
        sumcomp = sumcomp + (x(:, i).^2);
        prodcomp = prodcomp .* (cos(x(:, i) / sqrt(i)));
    end

    scores = (sumcomp / 4000) - prodcomp + 1;
end

function scores = happycatfcn(x, alpha)

    if nargin < 2
        alpha = 0.5;
    end

    n = size(x, 2);
    x2 = sum(x .* x, 2);
    scores = ((x2 - n).^2).^(alpha) + (0.5 * x2 + sum(x, 2)) / n + 0.5;
end

function scores = himmelblaufcn(x)
    n = size(x, 2);
    assert(n == 2, 'Himmelblau''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = ((X.^2 + Y - 11).^2) + ((X + (Y.^2) - 7).^2);
end

function scores = holdertablefcn(x)

    n = size(x, 2);
    assert(n == 2, 'The Holder-table function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    expcomponent = exp(abs(1 - (sqrt(X.^2 + Y.^2) / pi)));

    scores = -abs(sin(X) .* cos(Y) .* expcomponent);
end

function scores = keanefcn(x)
    n = size(x, 2);
    assert(n == 2, 'Keane function is defined only on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    numeratorcomp = (sin(X - Y).^2) .* (sin(X + Y).^2);
    denominatorcomp = sqrt(X.^2 + Y.^2);
    scores =- numeratorcomp ./ denominatorcomp;
end

function scores = leonfcn(x)
    n = size(x, 2);
    assert(n == 2, 'Leon function is defined only on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = 100 * ((Y - X.^3).^2) + ((1 - X).^2);
end

function scores = levin13fcn(x)
    n = size(x, 2);
    assert(n == 2, 'Levi''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);
    scores = sin(3 * pi * X).^2 + ...
        ((X - 1).^2) .* (1 + sin(3 * pi * Y).^2) + ...
        ((Y - 1).^2) .* (1 + sin(2 * pi * Y).^2);
end

function scores = matyasfcn(x)
    n = size(x, 2);
    assert(n == 2, 'Matyas''s function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = 0.26 * (X.^2 + Y.^2) - 0.48 * X .* Y;
end

function scores = mccormickfcn(x)

    n = size(x, 2);
    assert(n == 2, 'The McCormick function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = sin(X + Y) + ((X - Y).^2) - 1.5 * X + 2.5 * Y + 1;
end

function scores = periodicfcn(x)

    sin2x = sin(x).^2;
    sumx2 = sum(x.^2, 2);
    scores = 1 + sum(sin2x, 2) -0.1 * exp(-sumx2);

end

function scores = powellsumfcn(x)
    n = size(x, 2);
    absx = abs(x);

    scores = 0;

    for i = 1:n
        scores = scores + (absx(:, i).^(i + 1));
    end

end

function scores = qingfcn(x)
    n = size(x, 2);
    x2 = x.^2;

    scores = 0;

    for i = 1:n
        scores = scores + (x2(:, i) - i).^2;
    end

end

function scores = quarticfcn(x)

    n = size(x, 2);

    scores = 0;

    for i = 1:n
        scores = scores + i * (x(:, i).^4);
    end

    scores = scores + rand;
end

function f = rastriginfcn(x)
    n = size(x, 2);
    A = 10;
    f = (A * n) + (sum(x.^2 - A * cos(2 * pi * x), 2));
end

function scores = ridgefcn(x, d, alpha)

    if nargin < 3
        alpha = 0.5;
    end

    if nargin < 2
        d = 1;
    end

    x1 = x(:, 1);
    scores = x1 + d * (sum(x(:, 2:end).^2, 2).^alpha);
end

function scores = rosenbrockfcn(x)
    scores = 0;
    n = size(x, 2);
    assert(n >= 1, 'Given input X cannot be empty');
    a = 1;
    b = 100;

    for i = 1:(n - 1)
        scores = scores + (b * ((x(:, i + 1) - (x(:, i).^2)).^2)) + ((a - x(:, i)).^2);
    end

end

function scores = salomonfcn(x)
    x2 = x.^2;
    sumx2 = sum(x2, 2);
    sqrtsx2 = sqrt(sumx2);

    scores = 1 - cos(2 .* pi .* sqrtsx2) + (0.1 * sqrtsx2);
end

function scores = schaffern1fcn(x)
    n = size(x, 2);
    assert(n == 2, 'Schaffer function N. 1 is defined only on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    numeratorcomp = (sin((X.^2 + Y.^2).^2).^2) - 0.5;
    denominatorcomp = (1 + 0.001 * (X.^2 + Y.^2)).^2;
    scores = 0.5 + numeratorcomp ./ denominatorcomp;
end

function scores = schaffern2fcn(x)

    n = size(x, 2);
    assert(n == 2, 'The Schaffer N. 2 function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    sincomponent = sin((X.^2) - (Y.^2)).^2;

    scores = 0.5 + ((sincomponent - 0.5) ./ (1 + 0.001 * (X.^2 + Y.^2)).^2);
end

function scores = schaffern3fcn(x)
    n = size(x, 2);
    assert(n == 2, 'Schaffer function N. 3 is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    numeratorcomp = (sin(cos(abs(X.^2 - Y.^2))).^2) - 0.5;
    denominatorcomp = (1 + 0.001 * (X.^2 + Y.^2)).^2;
    scores = 0.5 + numeratorcomp ./ denominatorcomp;
end

function scores = schaffern4fcn(x)
    n = size(x, 2);
    assert(n == 2, 'Schaffer function N. 4 is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    numeratorcomp = (cos(sin(abs(X.^2 - Y.^2))).^2) - 0.5;
    denominatorcomp = (1 + 0.001 * (X.^2 + Y.^2)).^2;
    scores = 0.5 + numeratorcomp ./ denominatorcomp;
end

function scores = schwefel220fcn(x)
    scores = sum(abs(x), 2);
end

function scores = schwefel221fcn(x)
    scores = max(abs(x), [], 2);
end

function scores = schwefel222fcn(x)

    absx = abs(x);
    scores = sum(absx, 2) + prod(absx, 2);
end

function scores = schwefel223fcn(x)
    scores = sum(x.^10, 2);
end

function scores = schwefelfcn(x)
    n = size(x, 2);
    scores = 418.9829 * n - (sum(x .* sin(sqrt(abs(x))), 2));
end

function scores = shubert3fcn(x)
    n = size(x, 2);

    scores = 0;

    for i = 1:n

        for j = 1:5
            scores = scores + j * sin(((j + 1) * x(:, i)) + j);
        end

    end

end

function scores = shubert4fcn(x)
    n = size(x, 2);

    scores = 0;

    for i = 1:n

        for j = 1:5
            scores = scores + j * cos(((j + 1) * x(:, i)) + j);
        end

    end

end

function scores = shubertfcn(x)
    n = size(x, 2);

    scores = 1;

    for i = 1:n
        inner_sum = 0;

        for j = 1:5
            inner_sum = inner_sum + j * cos(((j + 1) * x(:, i)) + j);
        end

        scores = inner_sum .* scores;
    end

end

function f = spherefcn(x)
    f = sum(x.^2);
end

function scores = styblinskitankfcn(x)
    n = size(x, 2);
    scores = 0;

    for i = 1:n
        scores = scores + ((x(:, i).^4) - (16 * x(:, i).^2) + (5 * x(:, i)));
    end

    scores = 0.5 * scores;
end

function scores = sumsquaresfcn(x)

    [m, n] = size(x);
    x2 = x.^2;
    I = repmat(1:n, m, 1);
    scores = sum(I .* x2, 2);

end

function scores = threehumpcamelfcn(x)

    n = size(x, 2);
    assert(n == 2, 'The Three-hump camel function is only defined on a 2D space.')
    X = x(:, 1);
    Y = x(:, 2);

    scores = (2 * X.^2) - (1.05 * (X.^4)) + ((X.^6) / 6) + X .* Y + Y.^2;
end

function scores = wolfefcn(x)
    n = size(x, 2);
    assert(n == 3, 'The Wolfe function is defined only on the 3-D space.')
    X = x(:, 1);
    Y = x(:, 2);
    Z = x(:, 3);

    scores = (4/3) * (((X.^2 + Y.^2) - (X .* Y)).^(0.75)) + Z;
end

function scores = xinsheyangn1fcn(x)
    n = size(x, 2);

    scores = 0;

    for i = 1:n
        scores = scores + rand * (abs(x(:, i)).^i);
    end

end

function scores = xinsheyangn2fcn(x)
    scores = sum(abs(x), 2) .* exp(-sum(sin(x.^2), 2));
end

function scores = xinsheyangn3fcn(x)
    beta = 15;
    m = 5;

    scores = exp(-sum((x / beta).^(2 * m), 2)) - (2 * exp(-sum(x.^2, 2)) .* prod(cos(x).^2, 2));
end

function scores = xinsheyangn4fcn(x)
    scores = (sum(sin(x).^2, 2) - exp(-sum(x.^2, 2))) .* exp(-sum(sin(sqrt(abs(x))).^2, 2));
end

function scores = zakharovfcn(x)

    n = size(x, 2);
    comp1 = 0;
    comp2 = 0;

    for i = 1:n
        comp1 = comp1 + (x(:, i).^2);
        comp2 = comp2 + (0.5 * i * x(:, i));
    end

    scores = comp1 + (comp2.^2) + (comp2.^4);
end
