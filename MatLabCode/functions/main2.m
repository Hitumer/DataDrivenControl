clear 
close all
D_x = 2;
D_u = 2;

a = ones(D_x,1);
b = ones(D_u,1);
sysFull = drss(2,2,2);
A = sysFull.A;
B = sysFull.B;
L = 21;
iter_noise =20;
N_min = 101;
N_max = 2001;
N_iter = (N_max - N_min)/(L-1)+1;
% x_direct = DirectMethod(N,x_data{1,1},x_data{end},u_data{end},D_x,D_u,L);
% x_indirect = IndirectMethod(N,x_data{1,1},x_data{end},u_data{end},D_x,D_u);
x_direct = cell(N_iter,iter_noise);
x_indirect = cell(N_iter,iter_noise);

cost_array_direct = zeros(N_iter,iter_noise);
cost_array_indirect = zeros(N_iter,iter_noise);
cost_direct = cell(1,3);
cost_indirect = cell(1,3);
for k = 1:4
for j = 1:N_iter
    N = (j - 1)*(L-1)+N_min;
    data = DataGenerate(A,B,N,iter_noise,k);
    x_data = data{1,1};
    u_data = data{1,2};
    x_original = x_data{end};
    u_original = u_data{end};

    for i = 1:iter_noise
        x_direct{j,i} = DirectMethod(N,x_data{i,1},x_data{end},u_data{end},D_x,D_u,L);
        x_indirect{j,i} = IndirectMethod(N,x_data{i,1},x_data{end},u_data{end},D_x,D_u);
        cost_array_direct(j,i) = CostFunction(x_original,x_direct{j,i});
        cost_array_indirect(j,i) = CostFunction(x_original,x_indirect{j,i});
    end
    cost_direct{k} = cost_array_direct;
    cost_indirect{k} = cost_array_indirect;

end
end

% [m,n] = size(cost_array_direct);
% for i= 1:iter_noise
%     X(:,i) = i;
%     for j = 1:N_iter
%         X(j,:) = j;
%     end
% 
% end
% for j = 1:N_iter
%     Y(j,:) = j;
%     for i = 1:iter_noise
%         Y(:,i) = i;
%     end
% end
%cost_plot = reshape(cost_array_direct,[N*M,1]);

% surf(X,Y,cost_array_direct),shading flat;
% hold on
% surf(X,Y,cost_array_indirect);
% hold off

%************************************************
%Add noise to u input
% cost_array_direct2 = zeros(N_iter,iter_noise);
% cost_array_indirect2 = zeros(N_iter,iter_noise);
% for j = 1:N_iter
%     N = (j - 1)*(L-1)+N_min;
%     data = DataGenerate(A,B,N,iter_noise,2);
%     x_data = data{1,1};
%     u_data = data{1,2};
%     x_original = x_data{end};
%     u_original = u_data{end};
% 
%     for i = 1:iter_noise
%         x_direct{j,i} = DirectMethod(N,x_data{i,1},x_data{end},u_data{end},D_x,D_u,L);
%         x_indirect{j,i} = IndirectMethod(N,x_data{i,1},x_data{end},u_data{end},D_x,D_u);
%         cost_array_direct2(j,i) = CostFunction(x_original,x_direct{j,i});
%         cost_array_indirect2(j,i) = CostFunction(x_original,x_indirect{j,i});
%     end
% 
% end
% 
% cost_array_direct3 = zeros(N_iter,iter_noise);
% cost_array_indirect3 = zeros(N_iter,iter_noise);
% for j = 1:N_iter
%     N = (j - 1)*(L-1)+N_min;
%     data = DataGenerate(A,B,N,iter_noise,3);
%     x_data = data{1,1};
%     u_data = data{1,2};
%     x_original = x_data{end};
%     u_original = u_data{end};
% 
%     for i = 1:iter_noise
%         x_direct{j,i} = DirectMethod(N,x_data{i,1},x_data{end},u_data{end},D_x,D_u,L);
%         x_indirect{j,i} = IndirectMethod(N,x_data{i,1},x_data{end},u_data{end},D_x,D_u);
%         cost_array_direct3(j,i) = CostFunction(x_original,x_direct{j,i});
%         cost_array_indirect3(j,i) = CostFunction(x_original,x_indirect{j,i});
%     end
% 
% end
% figure
% surf(X,Y,cost_array_direct2),shading flat;
% hold on
% surf(X,Y,cost_array_indirect2);
% hold off


% for i = 1:N_iter
%     figure
%     hold on
%     plot(1:1:iter_noise,cost_array_direct(i,:),'r');
%     plot(1:1:iter_noise,cost_array_direct2(i,:),'g');
%     plot(1:1:iter_noise,cost_array_indirect(i,:),'b');
%     plot(1:1:iter_noise,cost_array_indirect2(i,:),'-black')
%     hold off
%     legen('Direct','Direct with input noise','Indirect','Indirect with input noise')
% end

for j = 5
    figure
    subplot(3,1,1)
    %title(['The added noise power is ',num2str(0.01*j),'power of signal'])
    hold on
    plot(1:1:N_iter,cost_direct{1}(:,j),'r','LineWidth',3);
    plot(1:1:N_iter,cost_indirect{1}(:,j),'b','LineWidth',3);
    xlabel('Datasize');
    ylabel('Error');
    legend('Direct with additive x disturbance','Indirect with x disturbance')

    subplot(3,1,2)
    hold on
    plot(1:1:N_iter,cost_direct{2}(:,j),'r','LineWidth',3);
    plot(1:1:N_iter,cost_indirect{2}(:,j),'b','LineWidth',3); 
    xlabel('Datasize');
    ylabel('Error');
    legend('Direct with additive u disturbance','Indirect with additive u disturbance')
  
    
%     subplot(2,2,3)
%     hold on
%     plot(1:1:N_iter,cost_direct{3}(:,j),'r','LineWidth',3);
%     plot(1:1:N_iter,cost_indirect{3}(:,j),'b','LineWidth',3);
%     xlabel('Datasize');
%     ylabel('Error');
%     hold off
%     legend('Direct with u noise','Indirect with u noise');

   
    subplot(3,1,3)
    hold on
    plot(1:1:N_iter,cost_direct{4}(:,j),'r','LineWidth',3);
    plot(1:1:N_iter,cost_indirect{4}(:,j),'b','LineWidth',3);
    xlabel('Datasize');
    ylabel('Error');
    hold off
    legend('Direct with measurement noise','Indirect with measurement noise');
end
for j = 1: iter_noise
    for i = 1:4
    array_noise_power_dir(j,i) = sum(cost_direct{i}(:,j))/length(cost_direct{1}(:,j));
    array_noise_power_indir(j,i) = sum(cost_indirect{i}(:,j))/length(cost_indirect{1}(:,j));
    end 
end

    figure
    subplot(3,1,1)
    %title(['The added noise power is ',num2str(0.01*j),'power of signal'])
    hold on
    plot(1:1:iter_noise,array_noise_power_dir(:,1),'r','LineWidth',3);
    plot(1:1:iter_noise,array_noise_power_indir(:,1),'b','LineWidth',3);
    xlabel('Noise Power procent of the Signal Power');
    ylabel('Final Error');
    legend('Direct with additive x disturbance','Indirect with x disturbance')

    subplot(3,1,2)
    hold on
    plot(1:1:iter_noise,array_noise_power_dir(:,2),'r','LineWidth',3);
    plot(1:1:iter_noise,array_noise_power_indir(:,2),'b','LineWidth',3);
    xlabel('Noise Power procent of the Signal Power');
    ylabel('Final Error');
    legend('Direct with additive u disturbance','Indirect with u disturbance')
  
    
%     subplot(2,2,3)
%     hold on
%     plot(1:1:iter_noise,array_noise_power_dir(:,3),'r','LineWidth',3);
%     plot(1:1:iter_noise,array_noise_power_indir(:,3),'b','LineWidth',3);
%     xlabel('Noise Power procent of the Signal Power');
%     ylabel('Final Error');
%     hold off
%     legend('Direct with u noise','Indirect with u noise');

   
    subplot(3,1,3)
    hold on
    plot(1:1:iter_noise,array_noise_power_dir(:,4),'r','LineWidth',3);
    plot(1:1:iter_noise,array_noise_power_indir(:,4),'b','LineWidth',3);
    xlabel('Noise Power procent of the Signal Power');
    ylabel('Final Error');
    hold off
    legend('Direct with measurement noise','Indirect with measurement noise');
save("data_4_20_3.mat");
