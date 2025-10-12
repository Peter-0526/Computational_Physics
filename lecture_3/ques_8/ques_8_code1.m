% 计算定积分\int_{0}^{5}\sqrt{x^2 + 2\sqrt{x}}dx
N0 = 10000000;
X = Random_generator_16807(2*N0);
fun = @(x) sqrt(x.^2 + 2.*sqrt(x));
gfun = @(x) sqrt(x.^(3/2) + 2).*4*5^(1/4);
N1 = [100,1000,10000,100000,1000000,10000000];

for i = 1 : 6
    N = N1(i);
    disp(['N = ',num2str(N)]);

    % 掷石法
    M1 = sqrt(25 + 2*sqrt(5));
    V = M1*5;
    X1 = X(1:N).*5;
    X2 = X(N+1:2*N).*M1;
    Index1 = find(X2 < fun(X1));
    ans1 = V.*length(Index1)/N;
    % 计算有效数字
    var1 = ans1*(V-ans1)/N;
    k1 = ceil(log10(ans1)) - floor(log10(sqrt(var1)));
    disp(['掷石法:   ','I1 = ',num2str(ans1),'   标准差',num2str(sqrt(var1)),...
        '   有效数字位数: ', num2str(k1)]);

    % 平均值法
    ans2 = 5.*sum(fun(X1))/N;
    % 计算有效数字
    ave_f2 = sum(fun(X1).^2)/N;     % <f^2>
    var2 = 25*(ave_f2-(ans2/5)^2)/N;
    k2 = ceil(log10(ans2)) - floor(log10(sqrt(var2)));
    disp(['平均值法:  ','I2 = ',num2str(ans2),'   标准差',num2str(sqrt(var2)),...
        '   有效数字位数: ', num2str(k2)]);

    % 重要抽样法
    Y1 = X(1:N).^(4/5)*5;
    ans3 = sum(gfun(Y1))./N;
    % 计算有效数字
    ave_fg2 = sum(gfun(Y1).^2)./N;  % <(f/g)^2>
    var3 = (ave_fg2 - ans3.^2)/N;
    k3 = ceil(log10(ans3)) - floor(log10(sqrt(var3)));
    disp(['重要抽样法: ','I3 = ',num2str(ans3),'   标准差',num2str(sqrt(var3)),...
        '   有效数字位数: ', num2str(k3)]);
end