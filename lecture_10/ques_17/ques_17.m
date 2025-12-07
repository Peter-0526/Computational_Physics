N = 1000000;
DN = 200;
fun = @(x,lambda) lambda*sin(pi*x);
min = 0;
max = 1;

bifurcation_plot(min, max, DN, N/100, fun);     % 绘制倍周期分叉图

% 求解分叉临界\lambda值及第一个Feigenbaum数
tol = 1e-10;
target = 2;     % 由1过渡到2的临界点
lambda = zeros(1,8);
for i = 1:length(lambda)
lambda(i) = find_Feigenbaum1(0.88, 0.5, tol, N, fun, 2^i);
end
Delta = zeros(1,length(lambda)-2);
for i = 1:length(Delta)
    Delta(i) = (lambda(i+1)-lambda(i))/(lambda(i+2)-lambda(i+1));
end
fprintf("lambda:\n")
for i = 1:length(lambda)
fprintf("%9f\n",lambda(i))
end
fprintf("Delta: \n")
for i = 1:length(Delta)
fprintf("%9f\n",Delta(i))
end

% 求解第二个Feigenbaum数
x_a = zeros(1,length(lambda)-1);
d = zeros(1,length(x_a));
alpha = zeros(1,length(lambda)-2);
for i = 1:length(lambda)-1
    [x_a(i),d(i)] = find_Feigenbaum2(lambda(i+1),lambda(i),tol,N,fun,i);
end
for i = 1:length(alpha)
    alpha(i) = d(i)/d(i+1);
end
fprintf("d: \n")
for i = 1:length(d)
fprintf("%9f\n",d(i))
end
fprintf("alpha: \n")
for i = 1:length(alpha)
fprintf("%9f\n",alpha(i))
end

function bifurcation_plot(max, min, n_plot, N, fun) % 绘制倍周期分叉图
dx = (max - min)/1000;
lambda = min:dx:max;         L = length(lambda);
X_Plot = repelem(lambda,n_plot);
Y_Plot = zeros(1,L*n_plot);
Counts = zeros(1,L);

figure
hold on
for i = 1:L
    x = zeros(1,N);
    x(1) = 0.5;
    for j = 2:N
        x(j) = fun(x(j-1),lambda(i));
    end
    for k = 1:n_plot
        if abs(x(1,N+k-n_plot) - x(1,N-n_plot)) < 1e-5
            break;
        end
    end
    Counts(i) = k;
    Y_Plot(1,n_plot*(i-1)+1 : n_plot*i) = x(1,N-n_plot+1:N);
end

% 绘制倍周期分叉图
scatter(X_Plot,Y_Plot,6,'blue','filled')
xlabel('\lambda','FontSize',18)
ylabel('X','FontSize',18)
title('$x_{n+1} = \lambda \sin(\pi x_n)$', 'Interpreter', 'latex','FontSize',18,'FontName', 'Arial')
end

function [X_inf,count] = x_cycle(lambda,N,fun,n_plot)% 计算周期数T及迭代结果X_inf
x = zeros(1,N+n_plot);
    x(1) = 0.5;
    for j = 2:N+n_plot
        x(j) = fun(x(j-1),lambda);
    end
    for k = 1:n_plot
        if abs(x(1,N+k) - x(1,N)) < 1e-6
            break;
        end
    end
    X_inf = sort(x(1,N:N+k-1));
    count = k;
end

function mid = find_Feigenbaum1(high, low, tol, N, fun, target) % 第一个Feigenbaum数计算
while high - low > tol
    mid = (low + high)/2;
    [~,count] = x_cycle(mid,N,fun,target);
    if count >= target
        high = mid;
    else
        low = mid;
    end
end
end

function [alpha,d] = find_Feigenbaum2(high, low, tol, N, fun, n) % 第二个Feigenbaum数计算
k = round((2^n + 2)/3).*(mod(n,2) == 0) + round((2^n + 4)/3).*(mod(n,2) == 1);
k1 = k - 1.*(mod(n,2) == 1);
X_inf = zeros(1,k1);
while abs(X_inf(k1)-0.5) > tol
    mid = (high + low)/2;
    alpha = mid;
    [X_inf,~] = x_cycle(mid,N,fun,2^n);
    if X_inf(k1) > 0.5 && mod(n,2) == 1 || X_inf(k1) < 0.5 && mod(n,2) == 0
        low = mid;
    else
        high = mid;
    end
end
d = X_inf(k) - X_inf(k-1);
end