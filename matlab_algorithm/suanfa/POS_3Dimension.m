
% 主调用函数
clc;clear;close all;
[ZuoBiao_x, ZuoBiao_y] = meshgrid(-10:0.1:10,-5:0.1:5);   % 产生二维坐标
ZuoBiao_z = f_xy(ZuoBiao_x, ZuoBiao_y);
figure(1);
s = mesh(ZuoBiao_x, ZuoBiao_y, ZuoBiao_z);      % 画网格曲面图
s.FaceColor = 'flat';    % 修改网格图的属性


% 初始化种群 
N = 100;                        % 初始种群个数  
D = 2;                          % 空间维数  
iter = 50;                      % 迭代次数       
x_limit = [-10, 10; -5, 5];     % 位置限制  
v_limit = [ -2,  2; -1, 1];     % 速度限制  

x = zeros(N, D);
for i = 1:D 
    x(:,i) = x_limit(i, 1) + (x_limit(i, 2) - x_limit(i, 1)) * rand(N, 1);%初始种群的位置  
end  
v(:,1) = rands(N, 1) * 1;       % 初始种群x方向的速度 
v(:,2) = rands(N, 1) * 2;       % 初始种群y方向的速度 

p_best = x;                     % (初始化)每个个体的历史最佳位置  
f_best = zeros(1, D);           % (初始化)种群的历史最佳位置  

fp_best = zeros(N, 1) - inf;    % (初始化)每个个体的历史最佳适应度为负无穷  
fg_best = -inf;                 % (初始化)种群历史最佳适应度为负无穷

w = 0.8;                        % 惯性权重
c1 = 0.5;                       % 自我学习因子  
c2 = 0.5;                       % 群体学习因子 

figure(2);
s = mesh(ZuoBiao_x, ZuoBiao_y, ZuoBiao_z);    % 画网格曲面图
s.FaceColor = 'flat';                         % 修改网格图的属性
hold on;
plot3(x(:,1), x(:,2),f_xy(x(:,1), x(:,2)), 'ro');
title('初始状态图');  


figure(3); 
i = 1;  
record = zeros(iter, 1);                 % 记录器  
while i <= iter  
    fx = f_xy(x(:,1), x(:,2));       % 个体当前适应度     
    for j = 1:N        
        if fp_best(j) < fx(j)        % 记录最大值
            fp_best(j) = fx(j);      % 更新个体历史最佳适应度  
            p_best(j,:) = x(j,:);    % 更新个体历史最佳位置  
        end   
    end  
    if fg_best < max(fp_best)  
        [fg_best, ind_max] = max(fp_best);    % 更新群体历史最佳适应度  
        f_best = p_best(ind_max, :);          % 更新群体历史最佳位置  
    end  
    
    v = v * w + c1 * rand() * (p_best - x) + c2 * rand() * (repmat(f_best, N, 1) - x); % 速度更新
    
    % 速度处理
    for t=1:N
        for k=1:D
            if v(t,k) > v_limit(k,2)      % 超速处理
                v(t,k) = v_limit(k,2);
            elseif v(t,k) < v_limit(k,1)  % 慢速处理
                v(t,k) = v_limit(k,1);
            end   
        end
    end

    x = x + v;            % 位置更新

    % 边界处理
    for t=1:N
        for k=1:D
            if x(t,k) > x_limit(k,2)      % 超过边界上限
                x(t,k) = x_limit(k,2);
            elseif x(t,k) < x_limit(k,1)  % 超过边界下限
                x(t,k) = x_limit(k,1);
            end
        end
    end
    
    record(i) = fg_best;   % 最大值记录  
    i = i + 1;
    if mod(i,10)  == 0
        i                  % 收敛进度输出
    end
    pause(0.1) 
    
end

figure(4)
plot(record);
title('收敛过程');

% 画最终状态图
figure(5);
[ZuoBiao_x, ZuoBiao_y] = meshgrid(-10:0.1:10,-5:0.1:5);   % 产生二维坐标 
s = mesh(ZuoBiao_x, ZuoBiao_y, ZuoBiao_z);      % 画网格曲面图
s.FaceColor = 'flat'; % 修改网格图的属性
hold on;
scatter3(x(:,1), x(:,2), f_xy(x(:,1), x(:,2)), 'ro');
title('最终状态图')

disp(['最大值：',num2str(fg_best)]);
disp(['变量取值：',num2str(f_best)]);


% 下面一段用来输出 粒子群最佳适应度的排行榜
[socre,ind] = sort(fp_best, 'descend')  % 降序排序
DaAn = zeros(3,2);
for i=1:N
    row = ind(i);
    DaAn(i,:) = p_best(row,:);
end

