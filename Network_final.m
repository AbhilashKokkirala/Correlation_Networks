clc; clear; close all;
%% Correlation Networks

Lqs = (load('w_final1').CL_w);     %.CL_w;
%Lqs = Lqs(3,:);
%Lqs = Lqs(1:2501);

rc = 0.8;
L = 40;

N = length(Lqs);
corr_matrix = zeros(N,N);
D = zeros(N,N);

for i = 1:N
    if i+L > N
        V_i = [Lqs(i:N), Lqs(1:i+L-N-1)];
    else
        V_i = Lqs(i:i+L-1);
    end
    
    for j = 1:N
        if j+L > N
            V_j = [Lqs(j:N), Lqs(1:j+L-N-1)];
        else
            V_j = Lqs(j:j+L-1);
        end

        corr = corrcoef(V_i, V_j);
        corr_matrix(i,j) = corr(1,2);

        if abs(corr(1,2)) > rc
            D(i,j) = 1;
        end

    end
end

for i = 1:N
    Degree(i) = sum(D(i,:));
end

for i = 1:N
    count = sum(Degree == i);
    num_nodes(i) = count;
end

alpha = power_law_coef(num_nodes)

[xx, yy] = meshgrid(1:N, 1:N);
figure(1)
pcolor(xx, yy, D)
shading interp    


figure(2)
plot(1:N, num_nodes)
xlabel('Degree')
ylabel('No of Nodes')

%%

% [DDF, D1] = histcounts(Degree, N);
% 
% figure(2)
% histogram(Degree, N)

%% FING POWER LAW COEFFICIENT ALPHA


function alpha = power_law_coef(DDF)
    %have to do log of it
    X = log(1:length(DDF));
    DDF1 = DDF;
    DDF1(DDF1 == 0) = 1;
    Y = log(DDF1);

    coef = polyfit(X, Y, 1);
    alpha = coef(1);
end

% % % function alpha = power_law_coef(DDF, Degree_vec)
% % %     %have to do log of it
% % %     X = log(Degree_vec);
% % %     DDF1 = DDF;
% % %     DDF1(DDF1 == 0) = 1;
% % %     Y = log(DDF1);
% % % 
% % %     coef = polyfit(X, Y, 1);
% % %     alpha = coef(1);
% % % end

% function DDF = DDF_func(Lqs, rc, L)
% N = length(Lqs);
% corr_matrix = zeros(N,N);
% D = zeros(N,N);
% 
% 
% for i = 1:N
%     if i+L > N
%         V_i = [Lqs(i:N), Lqs(1:i+L-N-1)];
%     else
%         V_i = Lqs(i:i+L-1);
%     end
%     
%     for j = 1:N
%         if j+L > N
%             V_j = [Lqs(j:N), Lqs(1:j+L-N-1)];
%         else
%             V_j = Lqs(j:j+L-1);
%         end
% 
%         corr = corrcoef(V_i, V_j);
%         corr_matrix(i,j) = corr(1,2);
% 
%         if corr(1,2) > rc
%             D(i,j) = 1;
%         end
% 
%     end
% end
% 
% Degree = zeros(N);
% for i = 1:N
%     Degree(i) = sum(D(i,:));
% end
% 
% [DDF, ~] = histcounts(Degree, N);
% 
% end
