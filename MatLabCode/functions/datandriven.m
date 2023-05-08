clear all
close all

%%%%%%%%% generate state x and input u
N = 100;
T = 1:N;
D_x = 2;
D_u = 2;
x = zeros(D_x,N);
u = zeros(D_u,N);
a = ones(D_x,1);
b = ones(D_u,1);
A = [1,1/N;1.01,0];
B = [1,0;0,1+1/N];
x_noise_rand = randn(D_x,N);
u_noise_rand = randn(D_u,N);
iter_noise = 100;
predict_step = N/2;
x(:,1) = [2;3];
for i = 1:N
    u(:,i) = [sin(pi/5*i);cos(i*pi/5)];
end
for j = 1:N-1
    x(:,j+1) = A*x(:,j)+B*u(:,j);
end
power_x = norm(x.^2);
power_u = norm(u.^2);
%% indirect method 


M_con_rec =zeros(D_x,D_x+D_u,iter_noise);

noise_power_total = ([1:iter_noise]*0.0001+0.0001)*sqrt(power_x);
noise_power_u = ([1:iter_noise]*0.001+0.0001)*power_u*0;
u_noise = noise_power_u(1);
for m = 1:iter_noise

    noise_power = noise_power_total(m);
    x_noise = x_noise_rand/norm(x_noise_rand)*noise_power;
    x_noise_add(:,:,m) = x_noise_rand/norm(x_noise_rand)*noise_power;
    x_bat=x+x_noise;
    x_bat_con(:,:,m) = x_bat;
    L = length(x(:,1));
    x_input1 = [x_bat(:,1:end-1);u(:,1:end-1)];
    x_input2 = x_bat(:,2:end);
    M_con = x_input2*pinv(x_input1);
    M_con_rec(:,:,m) = M_con;
    A_con = M_con(:,1:D_x);
    B_con = M_con(:,D_x+1:end);
    x_con(:,m,1) = x(:,1);
    for i=1:N-1
        x_con(:,m,i+1) = A_con* x(:,i) + B_con*u(:,i); 
    end
    x_rcon = reshape(x_con(:,m,:),2,N);
    SNR(m) = norm((x_bat(:,predict_step) - x_rcon(:,predict_step)).^2)/noise_power;
end
  
    
[SNR_opt,Idx] =min(SNR);
noise_power_opt = noise_power_total(Idx);
x_noise_opt = x_noise_rand/norm(x_noise_rand)*noise_power_opt;
x_bat_opt =x+x_noise_opt;
L = length(x(:,1));
x_input1 = [x_bat_opt(:,1:end-1);u(:,1:end-1)];
x_input2 = x_bat_opt(:,2:end);
M_con_opt = x_input2*pinv(x_input1);
A_con_opt = M_con_opt(:,1:D_x);
B_con_opt = M_con_opt(:,D_x+1:end);
x_con_opt = reshape(x_con(:,Idx,:),2,N);
%% direct method 


% 
% H_u = [[u(:,1),u(:,2),u(:,3),u(:,4),u(:,5)];
%         [u(:,2),u(:,3),u(:,4),u(:,5),u(:,6)];
%         [u(:,3),u(:,4),u(:,5),u(:,6),u(:,7)];
%         [u(:,4),u(:,5),u(:,6),u(:,7),u(:,8)];
%         [u(:,5),u(:,6),u(:,7),u(:,8),u(:,9)];
%         [u(:,6),u(:,7),u(:,8),u(:,9),u(:,10)]];
 x_0 = zeros(2,1);
x_hankel = zeros(D_x,N);

for i = 1:N
    x_hankel(:,i) =x(:,i);
    u_hankel(:,i) = u(:,i);
end

H_x = hankelmatrix(x_hankel,predict_step);
H_u = hankelmatrix(u_hankel,predict_step);
% H_x = [[x(:,1),x(:,2),x(:,3),x(:,4),x(:,5)];
%         [x(:,2),x(:,3),x(:,4),x(:,5),x(:,6)];
%         [x(:,3),x(:,4),x(:,5),x(:,6),x(:,7)];
%         [x(:,4),x(:,5),x(:,6),x(:,7),x(:,8)];
%         [x(:,5),x(:,6),x(:,7),x(:,8),x(:,9)];
%         [x(:,6),x(:,7),x(:,8),x(:,9),x(:,10)]];
% %         [x(:,7),x(:,8),x(:,9),x(:,10),x_0];
% %              [x(:,8),x(:,9),x(:,10),x_0,x_0];
% %              [x(:,9),x(:,10),x_0,x_0,x_0];
% %              [x(:,10),x_0,x_0,x_0,x_0]];
 for m = 1:iter_noise
    noise_power = noise_power_total(m);
    x_noise = x_noise_rand/norm(x_noise_rand)*noise_power;
    x_bat_direct =x+x_noise;

    state_init = zeros(D_u*(predict_step)+D_x,1);
for i = 1:predict_step
    state_init((i-1)*D_x+1:D_u*i,:) = u(:,i);
end
state_init(end-D_x+1:end,1) = x_bat_direct(:,1);
%     state_init = [u(:,1);u(:,2);u(:,3);u(:,4);u(:,5);u(:,6);x_bat_direct(:,1)];
    alpha = pinv([H_u;H_x(1:D_x,:)])*state_init;
    H_x_new = H_x(3:end,:);
   
    x_new = H_x_new*alpha;
    x_new = reshape(x_new,[2,predict_step-1]);
    x_predict = [x(:,1),x_new];
    x_predict_con(:,:,m) = x_predict;

    SNR_direct(m) = norm((x_bat_direct(:,length(x_predict)) - x_predict).^2)/noise_power;
    
 end
    
 [SNR_direct_opt, Idx_direct] = min(SNR_direct);
 x_direct_opt = x_predict_con(:,:,Idx_direct);


%% comparesion

Idx_same = find(SNR > SNR_direct);
if isempty(Idx_same)
    Idx_same =5;

end

%% 

subplot(1,2,1);
plot(1:iter_noise,SNR);
title('Changed noise power in Indirect method ')
hold on
plot(Idx,SNR_opt,'or');
hold off
subplot(1,2,2);
plot(1:iter_noise,SNR_direct,1:iter_noise,SNR);
hold on
plot(Idx_same(1),SNR(Idx_same(1)),'or');
title(['Nosie power same at ',num2str(noise_power_total(Idx_same(end))/power_x),'signal power']);
hold off
figure()
hold on
x_con_plot = reshape(x_con(:,Idx,:),[D_x,length(x_con(:,Idx,:))]);
x_predict = reshape(x_predict_con(:,:,Idx),[D_x,length(x_predict_con(:,:,Idx))]);
plot(T,x,'-*blue', T,x_bat_con(:,:,Idx),'-or',T,x_con_plot,'-xblack',T(1:length(x_predict)),x_predict,'-.green');
title(['With',num2str(noise_power_total(Idx)/power_x),'*signal power, optimization in indirect']);
h =legend({'orignal1','orginal2','AddNoise1','AddNoise2','indrect1','indirect2','direct1','direct2'});
set(h,'Fontsize',7);
hold off
figure()
hold on
x_bat_con_plot = reshape(x_bat_con(:,:,Idx_same(end)),[D_x,length(x_bat_con(:,:,Idx_same(end)))]);
Idx_same_plot = reshape(x_con(:,Idx_same(end),:),[D_x,length(x_con(:,Idx_same(end),:))]);
x_predict_con_polt = reshape(x_predict_con(:,:,Idx_same(end)),[D_x,length(x_predict_con(:,:,Idx_same(end)))]);
plot(T,x,'-*blue', T,x_bat_con_plot, '-or', T,Idx_same_plot,'-xblack',T(1:length(x_predict)),x_predict_con_polt,'-.green');
title(['With',num2str(noise_power_total(Idx_same(end))/power_x),'*signal power, same effect']);
h = legend({'orignal1','orginal2','AddNoise1','AddNoise2','indrect1','indirect2','direct1','direct2'});
set(h,'Fontsize',7);
hold off
figure()
plot(T,x,'-*blue', T,reshape(x_bat_con(:,:,Idx_same(1)),[D_x,length(x_bat_con(:,:,Idx_same(1)))]), '-or', T,reshape(x_con(:,Idx_same(1),:),[D_x,length(x_con(:,Idx_same(1),:))]),'-xblack',T(1:length(x_predict_con(:,:,Idx_same(1)))),reshape(x_predict_con(:,:,Idx_same(1)),[D_x,length(x_predict_con(:,:,Idx_same(1)))]),'-.green');
title(['With',num2str(noise_power_total(Idx_same(1))/power_x),'*signal power, direct better than indirect']);
h = legend({'orignal1','orginal2','AddNoise1','AddNoise2','indrect1','indirect2','direct1','direct2'});
set(h,'Fontsize',7);
