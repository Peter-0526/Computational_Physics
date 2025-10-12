% 计算定积分\int_{0}^{5}\sqrt{x^2 + 2\sqrt{x}}dx
N0 = 10000000;
N1 = [100,1000,10000,100000,1000000,10000000];
X = Random_generator_16807(5*N0);

fun = @(x,y,z,u,v) 5+x.^2-y.^2+3.*x.*y-z.^2+u.^3-v.^3;
V = 7/10 * 4/7 * 9/10 * 2 * 13/11;

% 平均值法
tic
for i = 1:6
    N = N1(i);
    X1 = 7/10*X(1:N);
    Y1 = 4/7*X(N+1:2*N);
    Z1 = 9/10*X(2*N+1:3*N);
    U1 = 2*X(3*N+1:4*N);
    V1 = 13/11*X(4*N+1:5*N);
    I = V * sum(fun(X1,Y1,Z1,U1,V1))./N;
% 计算有效位数
    ave_f2 = sum(fun(X1,Y1,Z1,U1,V1).^2)/N;     % <f^2>
    var = V^2*(ave_f2-(I./V)^2)/N;
    k = ceil(log10(I)) - floor(log10(sqrt(var)));
    disp(['N : ',num2str(N),'   I = ',num2str(I),'   标准差',num2str(sqrt(var)),...
        '   有效数字位数: ', num2str(k)])
end
toc
