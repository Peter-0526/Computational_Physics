tic
% 初始参量定义
tau = 1;    % 弛豫时间
f_E = 10;    % 电场力
w_E = 5*pi;    % 电磁场变化角频率

k = 10000;   % 粒子的总数
T = 5;     % 运动总时间
nT = max(100,round(50*w_E/pi));  % 单位时间分格数  
sigma = 1;
t = 0:1/nT:T;

V_x = zeros(T*nT+1,k);
V_y = zeros(T*nT+1,k);
V_self_x = zeros(T*nT+1,1);
V_self_y = zeros(T*nT+1,1);
[V_x(1,:),V_y(1,:)] = Random_Gauss_generator(k,sigma/sqrt(2));  % 粒子初始速度

% x-y方向作用力
funx = @(v,t,A) -v./tau + f_E.*sin(w_E*t) + A;
funy = @(v,t,A) -v./tau + A;

pause(2)
A = Random_Gauss_generator(2*k*T*nT,sqrt(nT)*sigma); % 随机力
% fourth-order Runge-Kutta
for i = 2:T*nT+1
    h = 1/nT;
    k_1 = h*funx(V_x(i-1,:),(i-1)/nT,A((2*i-4)*k+1:(2*i-3)*k));
    k_2 = h*funx(V_x(i-1,:)+k_1/2,(i-1)/nT+h/2,A((2*i-4)*k+1:(2*i-3)*k));
    k_3 = h*funx(V_x(i-1,:)+k_2/2,(i-1)/nT+h/2,A((2*i-4)*k+1:(2*i-3)*k));
    k_4 = h*funx(V_x(i-1,:)+k_3,(i-1)/nT+h,A((2*i-4)*k+1:(2*i-3)*k));
    V_x(i,:) = V_x(i-1,:) + (k_1+2*k_2+2*k_3+k_4)/6;
    k_1 = h*funy(V_y(i-1,:),(i-1)/nT,A((2*i-3)*k+1:(2*i-2)*k));
    k_2 = h*funy(V_y(i-1,:)+k_1/2,(i-1)/nT+h/2,A((2*i-3)*k+1:(2*i-2)*k));
    k_3 = h*funy(V_y(i-1,:)+k_2/2,(i-1)/nT+h/2,A((2*i-3)*k+1:(2*i-2)*k));
    k_4 = h*funy(V_y(i-1,:)+k_3,(i-1)/nT+h,A((2*i-3)*k+1:(2*i-2)*k));
    V_y(i,:) = V_y(i-1,:) + (k_1+2*k_2+2*k_3+k_4)/6;
end

V_1_ave = sum(V_x(1,:))./k;
V_2_ave = sum(V_x(1,:).^2 + V_y(1,:).^2)/k;

% x,y方向的自相关函数\
tfun = @(t) V_2_ave.*exp(-t/tau) + V_1_ave*f_E.*...
    (sin(w_E.*t)/tau - w_E.*cos(w_E.*t) + w_E.*exp(-t./tau))./(w_E^2 + 1/tau^2);
for i = 1:T*nT+1
    V_self_x(i) = sum(V_x(i,:).*V_x(1,:))/k;
    V_self_y(i) = sum(V_y(i,:).*V_y(1,:))/k;
end
figure
plot(t,V_self_x,'LineWidth',1);
hold on
plot(t,V_self_y,'LineWidth',1);
plot(t,V_self_x + V_self_y,'LineWidth',2);
plot(t,tfun(t),'LineWidth',2)
legend('C_x','C_y','C_{total}','C_{theory}','FontSize',16)
title(['速度自相关函数C(t)图像','(f_E = ',num2str(f_E),',w_E = ',num2str(w_E),')'],'FontSize',18)
xlabel('t/s','FontSize',18)

toc