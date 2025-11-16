% Metropolis
N = 100000;
N_heat = 5000;
dr = 0.9;
beta = 5;
count = 0;
xi = Random_generator_16807(3*N);
xi_x = xi(1:N);     xi_y = xi(N+1:2*N);     xi_r = xi(2*N+1:3*N);
x = zeros(1,N);
y = zeros(1,N);
pfun = @(X,Y) exp(-beta.*((-(X.^2 + Y.^2) + (X.^4 + Y.^4)/2 + (X - Y).^4/3)));
x2_pfun = @(X,Y) X.^2.*exp(-beta.*((-(X.^2 + Y.^2) + (X.^4 + Y.^4)/2 + (X - Y).^4/3)));

C0 = integral2(pfun,-Inf,Inf,-Inf,Inf,'AbsTol',1e-12);
exact = integral2(x2_pfun,-Inf,Inf,-Inf,Inf,'AbsTol',1e-12)/C0;

surf_plot(pfun,beta); % p(x,y)平面分布图

for i = 2:N
    x_try = x(i-1) + dr*(xi_x(i) - 0.5);
    y_try = y(i-1) + dr*(xi_y(i) - 0.5);
    r = pfun(x_try,y_try)/pfun(x(i-1),y(i-1));
    if min(1,r) > xi_r(i)
        x(i) = x_try;
        y(i) = y_try;
        count = count + 1;
    else
        x(i) = x(i-1);
        y(i) = y(i-1);
    end
end

accept_rate = count/(N-1);
% burn_period(N, x, y);       % <x^2>计算结果与N关系图，用于估算热化长度

x2_mean = sum(x(N_heat+1:N).^2)/(N-N_heat);
y2_mean = sum(y(N_heat+1:N).^2)/(N-N_heat);
x2y2_mean = x2_mean + y2_mean;
err_x2 = abs(x2_mean - exact)/exact;
err_y2 = abs(y2_mean - exact)/exact;
err_x2y2 = abs(x2y2_mean - 2*exact)/(2*exact);

disp(['<x^2> = ',num2str(x2_mean),' <y^2> = ',num2str(y2_mean),' <x^2 + y^2> = ',num2str(x2y2_mean)])
disp(['exact = ',num2str(exact),' err<x^2> = ',num2str(err_x2),' err<y^2> = ',num2str(err_y2),' err<x^2+y^2> = ',num2str(err_x2y2)])
scatter_plot(x, y, beta);       % 散点分布图
chain_plot(x, y, beta);         % Markov Chain图

function surf_plot(fun,beta)     % p(x,y)平面分布图
[X,Y] = meshgrid(-2.5:0.1:2.5);
Z = fun(X,Y);
surf(X,Y,Z)
xlabel('X','FontSize',18)
ylabel('Y','FontSize',18)
title(['p（x,y）平面分布图,\beta = ',num2str(beta)],'FontSize',18)
figure
end

function burn_period(N, x, y)   % <x^2>计算结果与N关系图，用于估算热化长度
x_test = 1:1000:N;
l = length(x_test);
x2_mean = zeros(1,l);
y2_mean = zeros(1,l);
x2y2_mean = zeros(1,l);
for i = 1:l
    x2_mean(i) = sum(x(1:x_test(i)).^2)/x_test(i);     % <x^2>
    y2_mean(i) = sum(y(1:x_test(i)).^2)/x_test(i);     % <y^2>
    x2y2_mean(i) = x2_mean(i) + y2_mean(i);   % <x^2 + y^2>
end
plot(x_test,x2_mean,x_test,y2_mean,x_test,x2y2_mean,'LineWidth',2)
xlabel('N','FontSize',18)
title('<x^2>,<y^2>,<x^2+y^2>计算结果随序列长度关系图（未舍去热化长度）','FontSize',18)
legend('<x^2>','<y^2>','<x^2+y^2>','FontSize',18)
figure
end

function scatter_plot(x, y, beta)   % 散点分布图
scatter(x(10001:20000),y(10001:20000),10,'filled')
xlabel('X','FontSize',18)
ylabel('Y','FontSize',18)
title(['Markov链散点分布图（N = 10000, \beta = ',num2str(beta), '）' ],'FontSize',18)
figure
end

function chain_plot(x, y, beta)     % Markov Chain图
plot(x(1000:5000),y(1000:5000),'LineWidth',0.2)
xlabel('X','FontSize',18)
ylabel('Y','FontSize',18)
title(['Markov链（N = 4000, \beta = ',num2str(beta), '）' ],'FontSize',18)
end
