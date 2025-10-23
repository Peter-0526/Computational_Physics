% 二维网格随机游走
N = 1000;
ra = Random_generator_16807(N);     % 随机数

X = zeros(N+1,1);
Y = zeros(N+1,1);

for i = 2:N+1
    if ra(i-1) < 0.25
        X(i) = X(i-1) + 1;
        Y(i) = Y(i-1);
    elseif  ra(i-1) < 0.5
        X(i) = X(i-1);
        Y(i) = Y(i-1) + 1;
    elseif ra(i-1) < 0.75
        X(i) = X(i-1) - 1;
        Y(i) = Y(i-1);
    else
        X(i) = X(i-1);
        Y(i) = Y(i-1) - 1;
    end
end

plot(X,Y)

