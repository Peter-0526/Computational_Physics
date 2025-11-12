% 初始参量定义
N_heat = 2000;  % 热化长度
N = 500000;     % 有效长度

Alpha = 2:2:10;
Beta = 2:2:10;
nbins = 100;
err_1 = zeros(length(Alpha),length(Beta));
accept_1 = zeros(length(Alpha),length(Beta));
err_2 = zeros(length(Alpha),length(Beta));
accept_2 = zeros(length(Alpha),length(Beta));

for j = 1:length(Beta)
    for i = 1:length(Alpha)
        alpha = Alpha(i);               beta = Beta(j);
        gamma_1 = beta*(alpha + 6);     gamma_2 = beta*(alpha + 14);
        disp(['alpha = ',num2str(alpha),' beta = ',num2str(beta),' gamma_1 = ', ...
            num2str(gamma_1),' gamma_2 = ',num2str(gamma_2)])
        fun1 = @(x) 1/(beta*gamma(alpha)).*(x/beta).^(alpha - 1).*exp(-x./beta);     % 权重函数1
        fun2 = @(x) 1/(beta*gamma(alpha)).*(x/beta).^(alpha - 1).*exp(-x./beta)...  % 权重函数2
            .*(x - alpha*beta).^2;
        tfun1 = @(x) exp(-x./gamma_1)/gamma_1;   % 提议分布
        i_tfun1 = @(x) -gamma_1.*log(x);         % 根据提议分布采样
        tfun2 = @(x) exp(-x./gamma_2)/gamma_2;
        i_tfun2 = @(x) -gamma_2.*log(x);
        x_Max = 4*beta*alpha;   dx = x_Max/nbins;
        x_plot = dx/2:dx:x_Max-dx/2;

        % p(x) = f(x)抽样
        [x1,accept_1(i,j)] = metropolis_Sample(fun1,tfun1,i_tfun1,N,N_heat);
        s1 = histcounts(x1,nbins,'Normalization','pdf');
        sum1 = sum((x1 - alpha*beta).^2)/N;
        err_1(i,j) = abs(sum1/(alpha*beta^2) - 1);
        disp(['N = ',num2str(N),' sum1 = ',num2str(sum1),' relative error = ',num2str(err_1(i,j)), ...
            ' accept_rate = ', num2str(accept_1(i,j))])

        % p(x) = f(x) * (x - \alpha\beta) ^ 2
        [x2,accept_2(i,j)] = metropolis_Sample(fun2,tfun2,i_tfun2,N,N_heat);
        [s2,edges] = histcounts(x2,nbins,'Normalization','pdf','BinLimits',[0,x_Max]);
        S2 = zeros(1,nbins);
        for k = 1:nbins
            if s2(k) ~= 0
                S2(k) = fun2(x_plot(k))^2/s2(k);
            end
        end
        sum2 = sum(S2)/(sum(fun2(x_plot)));
        err_2(i,j) = abs(sum2/(alpha*beta^2) - 1);
        disp(['N = ',num2str(N),' sum2 = ',num2str(sum2),' relative error = ', ...
            num2str(err_2(i,j)),' accept_rate = ', num2str(accept_2(i,j))])
    end
end

filename = 'testdata.xlsx';
writematrix(Alpha',filename,'Sheet',2,'Range','A3:A7')
writematrix(Beta,filename,'Sheet',2,'Range','C1:G1')
writematrix(err_1,filename,'Sheet',2,'Range','C3:G7')
writematrix(Alpha',filename,'Sheet',2,'Range','A13:A17')
writematrix(Beta,filename,'Sheet',2,'Range','C11:G11')
writematrix(err_2,filename,'Sheet',2,'Range','C13:G17')

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

