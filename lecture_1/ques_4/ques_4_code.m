tic

% 输入a,b,c,初始参数定义
a = -0.5;b = 0.2;c = -1;
pfun = @(x)(x >= -1 & x < 0).*...
    (1/2 - a/2 + b /(c^2) * sinh(c) - b/c*exp(- c*x))+...
    (x >= 0 & x <= 1).*...
    (1/2 + a/2 + b /(c^2) * sinh(c) - b/c*exp(- c*x));
N1 = 100;
x = linspace(-1,1,N1);
y = pfun(x);
plot(x,y,'r','LineWidth',2);
hold on;

M_1 = max(y(1:N1/2)) + 0.05;
M_2 = max(y(N1/2+1:N1)) + 0.05;
Fun1 = @(x) (x >= 0 & x < M_1/(M_1+M_2)).*((1+M_2/M_1)*x - 1)+...
    (x >= M_1/(M_1+M_2) & x < 1).*((1 + M_1/M_2) * x - M_1/M_2);
Fun2 = @(x)(x >= -1 & x < 0).*M_1 + (x >= 0 & x < 1).*M_2;
% 绘制p(x) 判断p(x)是否满足条件


if(min(y) < 0)
    disp('不符合条件')
else    % 符合条件则继续

    % 产生随机数
    a1 = 16807;m = 2147483647;q = 127773;r = 2836;
    N = 20000000;
    I = zeros(1,N);
    cl = fix(clock());
    I(1) = cl(1) + 70*(cl(2) + 12 * (cl(3) +...
    31 * (cl(4) + 23 * (cl(5) + 59 * cl(6)))));

    for i = 2 : N+2
        I(i) = a1 * mod(I(i-1), q) - r * floor(I(i-1)/q);
        if I(i) < 0
            I(i) = I(i) + m;
        end
    end
    I = I./m;

    % xi_1 : I(1:N/2);xi_2 : I(N/2+1:N)
    xi_x(1,:) = Fun1(I(1:N/2));
    xi_x(2,:) = I(N/2+1:N).*Fun2(xi_x(1,:));
    scatter(xi_x(1,1:1000),xi_x(2,1:1000),'filled','b')
    legend('p(x)','散点分布')

    % 判断 若xi_x(2,i) <= pfun(xi_x(1,i)),保留xi_x
    Ifind = find(xi_x(2,:) <= pfun(xi_x(1,:)));
    xx = xi_x(1,Ifind);
    figure;
    histogram(xx,100,'Normalization','pdf');
    hold on;
    plot(x,y,'r','LineWidth',2);
    legend('拟合曲线','p(x)')
end

toc