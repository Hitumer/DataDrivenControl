clear all
close all

%%%%%%%%% generate state x and input u
N = 10;
T = 1:N;
D_x = 2;
D_u = 2;
x = zeros(D_x,N);
u = zeros(D_u,N);
a = ones(D_x,1);
b = ones(D_u,1);
A = [1,0.1;1,0];
B = [1,0;0,2];
x(:,1) = [2;3];

for i = 1:N
    u(:,i) = [sin(pi/5*i);cos(i*pi/5)];
end
for j = 1:N-1
    x(:,j+1) = A*x(:,j)+B*u(:,j);
end

%% indirect method 

power_x = norm(x.^2);
x_noise_rand = randn(D_x,N);
iter_noise = 1000;
M_con_rec =zeros(D_x,D_x+D_u,iter_noise);
for m = 1:iter_noise
    noise_power = (m*0.0001+0.001)*power_x;
    x_noise = x_noise_rand/norm(x_noise_rand)*noise_power;
    x_bat =x+x_noise;
    L = length(x(:,1));
    x_input1 = [x_bat(:,1:end-1);u(:,1:end-1)];
    x_input2 = x_bat(:,2:end);
    M_con = x_input2*pinv(x_input1);
    M_con_rec(:,:,m) = M_con;
    A_con = M_con(:,1:D_x);
    B_con = M_con(:,D_x+1:end);
    x_con(:,m,1) = x(:,1);
    for i=1:N-1
        x_con(:,m,i+1) = A_con* x_con(:,m,i) + B_con*u(:,i); 
    end
    x_rcon = reshape(x_con(:,m,:),2,10);
    SNR(m) = norm((x_bat - x_rcon).^2)/noise_power; 
end
    

[SNR_opt,Idx] =min(SNR);

noise_power_opt = (Idx*0.0001+0.001)*power_x;
x_noise_opt = x_noise_rand/norm(x_noise_rand)*noise_power_opt;
x_bat_opt =x+x_noise_opt;
L = length(x(:,1));
x_input1 = [x_bat_opt(:,1:end-1);u(:,1:end-1)];
x_input2 = x_bat_opt(:,2:end);
M_con_opt = x_input2*pinv(x_input1);
A_con_opt = M_con_opt(:,1:D_x);
B_con_opt = M_con_opt(:,D_x+1:end);
x_con_opt = reshape(x_con(:,Idx,:),2,10);
%% direct method 
% state_init = zeros(1,22);
% H_u = zeros(2,2*N);
% H_x = hankel(x(:,1));
% for i = 1:N
%     H_u(:,(i-1)*2+1:i*2) = hankel(u(:,i));
%     state_init(1,i) = u(1,i);
%     state_init(1,i+N+1) = u(2,i);
% end
% state_init(1,11) = x(1,1);
% state_init(1,22) = x(1,2);
% 
% alpha = state_init*pinv([H_u,H_x]);
% H_xstate = zeros(2,2*(N-1));
% for i = 1:N-1
%     H_xstate(:,(i-1)*2+1:i*2) = hankel(x(:,i+1));
% 
%     end
% x_state = alpha*H_xstate;

% figure();
% plot(x_recon,u);
% c = x_input2*pinv(x_input1);
subplot(1,2,1);
plot(1:m,SNR);
hold on
plot(Idx,SNR_opt,'ro');
hold off
subplot(1,2,2);
hold on
plot(T,x,'*blue');
plot(T,x_bat_opt,'or');
plot(T,x_con_opt,'xblack');
legend({'orignal1','orginal2','noise1','noise2','reconstucture1','reconstucture2'});
hold off
