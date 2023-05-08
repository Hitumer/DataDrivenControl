load("data_test.mat");
for k = 1:3
for j = 1:N_iter
    N = (j - 1)*(L-1)+N_min;
    data = DataGenerate(A,B,N,iter_noise,k);
    x_data = data{1,1};
    u_data = data{1,2};
    x_original = x_data{end};
    u_original = u_data{end};

    for i = 1:iter_noise
  
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

for j = 1:iter_noise
    figure
    title(['The added noise power is ',num2str(0.01*j),'power of signal'])
    hold on
    plot(1:1:N_iter,cost_direct{1}(:,j),'r');
    plot(1:1:N_iter,cost_direct{2}(:,j),'--r');
    plot(1:1:N_iter,cost_direct{3}(:,j),'-ro');
    plot(1:1:N_iter,cost_indirect{1}(:,j),'b');
    plot(1:1:N_iter,cost_indirect{2}(:,j),'--b');
    plot(1:1:N_iter,cost_indirect{3}(:,j),'-bo');
    hold off
    legend('Direct with x noise','Direct with input disturbance','Direct with u noise','Indirect with x noise','Indirect with input disturbance','Indirect with u noise');
end
hold off
save('data_test.mat','x_direct','x_indirect');
