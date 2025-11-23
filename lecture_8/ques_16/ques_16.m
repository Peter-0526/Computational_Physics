N =100000;
m = 5000;
T = [0.5:0.5:2,2.05:0.05:2.5,2.6,2.8,3:0.5:5];
H = -5:0.5:5;
t = length(T);
h = length(H);
xi = Random_generator_16807(3*N*t);
l = 20;
M_1 = zeros(t,1);
M_2 = zeros(t,1);
chi = zeros(t,1);
E_1 = zeros(t,1);
E_2 = zeros(t,1);
C_V = zeros(t,1);

% 外加磁场强度为0，改变温度，研究相变
for i = 1:t
    [sigma, M, E] = simulate_Ising(l, N, m, xi(1+3*N*(i-1):3*N*i), T(i),0 ,0);
    M_1(i) = sum(M)/(N-m);
    M_2(i) = sum(M.^2)/(N-m);
    chi(i) = (M_2(i) - M_1(i)^2)/T(i);
    E_1(i) = sum(E)/(N-m);
    E_2(i) = sum(E.^2)/(N-m);
    C_V(i) = (E_2(i) - E_1(i)^2)/T(i)^2*(l^2);
end
chi_plot(T,chi,l)
CV_plot(T,C_V,l)
magne_plot(T,M_1,l)
energy_plot(T,E_1,l)

% 不同温度下的磁滞回线
T2 = [1.5,2,2.25,3];
t2 = length(T2);
for j = 1:t2
    M_1_ = zeros(h,1);
    for i = 1:h
        [sigma, M, E] = simulate_Ising(l, N, m, xi(1+3*N*(i-1):3*N*i), T2(j), H(i),-1);
        M_1_(i) = sum(M)/(N-m);
    end
    mag_cur_plot(H,M_1_,l,T2(j))
    hold on
    for i = 1:h
        [sigma, M, E] = simulate_Ising(l, N, m, xi(1+3*N*(i-1):3*N*i), T2(j), H(i),1);
        M_1_(i) = sum(M)/(N-m);
    end
    mag_cur_plot(H,M_1_,l,T2(j))
    legend('-H \rightarrow +H','+H \rightarrow -H','FontSize',14)
    figure
end

function [sigma, M, E] = simulate_Ising(l, N, m, xi, T, H, init) % Matropolis抽样方法，
% 生成 末态自旋位形sigma, 磁场强度M, 平均能量E   
if init >= 0       
    sigma = ones(l,l);
else
    sigma = -ones(l,l);
end
M = zeros(N-m,1);
E = zeros(N-m,1);
E_now = -2-H;
for k = 1:N
    i = ceil(xi(k)*l);
    j = ceil(xi(k+N)*l);
    dE = sigma(i,j) * H;
    if i ~= 1 && i ~= l
        dE = dE + 2*sigma(i,j)*(sigma(i-1,j) + sigma(i+1,j));
    elseif i == 1
        dE = dE + 2*sigma(i,j)*(sigma(l,j) + sigma(i+1,j));
    else
        dE = dE + 2*sigma(i,j)*(sigma(i-1,j) + sigma(1,j));
    end
    if j ~= 1 && j ~= l
        dE = dE + 2*sigma(i,j)*(sigma(i,j-1) + sigma(i,j+1));
    elseif j == 1
        dE = dE + 2*sigma(i,j)*(sigma(i,j+1) + sigma(i,l));
    else
        dE = dE + 2*sigma(i,j)*(sigma(i,j-1) + sigma(i,1));
    end
    r = xi(k+2*N);
    if(r < min(exp(-dE/T),1))
        sigma(i,j) = -sigma(i,j);
        E_now = E_now + dE/(l^2);
        if k == m+1
            M(1) = sum(sum(sigma))/(l^2);
            E(1) = E_now;
        elseif k > m
            M(k-m) = M(k-m-1) + 2*sigma(i,j)/(l^2) ;
            E(k-m) = E_now;
        end
    else
        if k == m+1
            M(1) = sum(sum(sigma))/(l^2);
            E(1) = E_now;
        elseif k > m
            M(k-m) = M(k-m-1);
            E(k-m) = E_now;
        end
    end
end
end

function mag_plot(M,T)      % 绘制磁场强度M随操作次数变化图像
plot(M)
xlabel('Time','FontSize',18,'FontName', 'Arial')
ylabel('Magnetization','FontSize',18,'FontName', 'Arial')
title(['磁场强度M随操作次数变化(去除热化长度,T = ',num2str(T),')'],'FontSize',18,'FontName', 'Arial')
figure
end

function elec_plot(E,T)     % 绘制能量E随操作次数变化图像
plot(E)
xlabel('Time','FontSize',18,'FontName', 'Arial')
ylabel('energy per spin','FontSize',18,'FontName', 'Arial')
title(['能量E随操作次数变化(去除热化长度,T = ',num2str(T),')'],'FontSize',18,'FontName', 'Arial')
figure
end

function lattice_plot(l,sigma)   % 绘制格点自旋构型
[po_x,po_y] = find(sigma - 1);
[ne_x,ne_y] = find(sigma + 1);
scatter(po_x,po_y,14,'filled','Color','blue')
hold on
scatter(ne_x,ne_y,14,'filled','Color','red')
xlabel('X','FontSize',18)
ylabel('Y','FontSize',18)
title([num2str(l),'\times',num2str(l),'格点2维Ising模型自旋构型' ],'FontSize',18,'FontName', 'Arial')
figure
end

function chi_plot(T,chi,l)      % 绘制磁化率\chi随温度变化关系图
plot(T,chi,'o-','LineWidth',1)
xlabel('T','FontSize',18,'FontName', 'Arial')
ylabel('\chi','FontSize',18,'FontName', 'Arial')
title(['磁化率\chi与温度关系示意图(格点边长',num2str(l),'\times',num2str(l),')'],'FontSize',18,'FontName', 'Arial')
figure
end

function CV_plot(T,C_V,l)       % 绘制C_V随温度变化关系图
plot(T,C_V,'o-','Color','r','LineWidth',1)
xlabel('T','FontSize',18,'FontName', 'Arial')
ylabel('C_V','FontSize',18,'FontName', 'Arial')
title(['C_V与温度关系示意图(格点边长',num2str(l),'\times',num2str(l),')'],'FontSize',18,'FontName', 'Arial')
figure
end

function magne_plot(T,M_1,l)    % 绘制磁场强度随温度变化关系图
plot(T,M_1,'o-','Color','m','LineWidth',1)
xlabel('T','FontSize',18,'FontName', 'Arial')
ylabel('M','FontSize',18,'FontName', 'Arial')
title(['<M>与温度关系示意图(格点边长',num2str(l),'\times',num2str(l),')'],'FontSize',18,'FontName', 'Arial')
figure
end

function energy_plot(T,E_1,l)   % 绘制平均能量随温度变化关系图
plot(T,E_1,'o-','Color','g','LineWidth',1)
xlabel('T','FontSize',18,'FontName', 'Arial')
ylabel('<E>','FontSize',18,'FontName', 'Arial')
title(['<E>与温度关系示意图(格点边长',num2str(l),'\times',num2str(l),')'],'FontSize',18,'FontName', 'Arial')
figure
end

function mag_cur_plot(H,M_1,l,T)    % 绘制磁滞回线
plot(H,M_1,'o-','LineWidth',1)
xlabel('H','FontSize',18,'FontName', 'Arial')
ylabel('<M>','FontSize',18,'FontName', 'Arial')
title(['<M>与温度关系示意图(格点边长',num2str(l),'\times',num2str(l),',T = ',num2str(T),')'], ...
    'FontSize',18,'FontName', 'Arial')
end