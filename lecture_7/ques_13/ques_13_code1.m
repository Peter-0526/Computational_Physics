% 初始参量定义
N_heat = 2000;  % 热化长度
N = 500000;     % 有效长度

alpha = 1;
beta = 1;
Gamma_1 = [0.5:0.5:5,6:1:15,20:5:100];
nbins = 100;
err_1 = zeros(1,length(Gamma_1));
accept_1 = zeros(1,length(Gamma_1));
err_2 = zeros(1,length(Gamma_1));
accept_2 = zeros(1,length(Gamma_1));

for k = 1:length(Gamma_1)
    gamma_1 = Gamma_1(k);
    disp(['alpha = ',num2str(alpha),' beta = ',num2str(beta),' gamma = ',num2str(gamma_1)])
    fun1 = @(x) 1/(beta*gamma(alpha)).*(x/beta).^(alpha - 1).*exp(-x./beta);     % 权重函数1
    fun2 = @(x) 1/(beta*gamma(alpha)).*(x/beta).^(alpha - 1).*exp(-x./beta)...  % 权重函数2
        .*(x - alpha*beta).^2;
    tfun = @(x) exp(-x./gamma_1)/gamma_1;   % 提议分布
    i_tfun = @(x) -gamma_1.*log(x);         % 根据提议分布采样
    x_Max = 4*beta*alpha;   dx = x_Max/nbins;
    x_plot = dx/2:dx:x_Max-dx/2;

    % p(x) = f(x)抽样
    [x1,accept_1(k)] = metropolis_Sample(fun1,tfun,i_tfun,N,N_heat);
    s1 = histcounts(x1,nbins,'Normalization','pdf');
    sum1 = sum((x1 - alpha*beta).^2)/N;
    err_1(k) = abs(sum1/(alpha*beta^2) - 1);
    disp(['N = ',num2str(N),' sum1 = ',num2str(sum1),' relative error = ',num2str(err_1(k)), ...
    ' accept_rate = ', num2str(accept_1(k))])
    
    % p(x) = f(x) * (x - \alpha\beta) ^ 2
    [x2,accept_2(k)] = metropolis_Sample(fun2,tfun,i_tfun,N,N_heat);
    [s2,edges] = histcounts(x2,nbins,'Normalization','pdf','BinLimits',[0,x_Max]);
    S2 = zeros(1,nbins);
    for i = 1:nbins
        if s2(i) ~= 0
            S2(i) = fun2(x_plot(i))^2/s2(i);
        end
    end
    sum2 = sum(S2)/(sum(fun2(x_plot)));
    err_2(k) = abs(sum2/(alpha*beta^2) - 1);
    disp(['N = ',num2str(N),' sum2 = ',num2str(sum2),' relative error = ', ...
    num2str(err_2(k)),' accept_rate = ', num2str(accept_2(k))])
end

plot(Gamma_1,err_1,Gamma_1,accept_1,'LineWidth',2);
legend('Relative error1','Accept rate1','FontSize',18)
title(['p(x) = f(x) 抽样结果图(\alpha = ',num2str(alpha),' \beta = ',num2str(beta),')'],'FontSize',20)
figure;
plot(Gamma_1,err_2,Gamma_1,accept_2,'LineWidth',2);
legend('Relative error2','Accept rate2','FontSize',18)
title(['p(x) = (x - \alpha\beta)^2f(x) 抽样结果图(\alpha = ',num2str(alpha),' \beta = ',num2str(beta),')'],'FontSize',20)

filename = 'testdata.xlsx';
writematrix(Gamma_1',filename,'Sheet',1,'Range','A3:A39')
writematrix(err_1',filename,'Sheet',1,'Range','B3:B39')
writematrix(accept_1',filename,'Sheet',1,'Range','C3:C39')
writematrix(err_2',filename,'Sheet',1,'Range','D3:D39')
writematrix(accept_2',filename,'Sheet',1,'Range','E3:E39')

function [x,accept_r] = metropolis_Sample(func,tfun,i_tfun,N,N_heat)
    x = zeros(1,N+N_heat);
    xi = Random_generator_16807(2*N+2*N_heat);
    x(1) = 0;
    count = 0;
    
    for i = 2:N+N_heat
        x_try = i_tfun(xi(i));
        r = func(x_try).*tfun(x(i-1))/(func(x(i-1))*tfun(x_try));
        if r >= min(1,xi(i+N+N_heat))
            x(i) = x_try;
            count = count + 1;
        else
            x(i) = x(i-1);
        end
    end
    accept_r = count/(N+N_heat-1);
    x = x(N_heat+1:N_heat+N);
end

