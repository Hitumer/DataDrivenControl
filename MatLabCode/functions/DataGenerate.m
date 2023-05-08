function data = DataGenerate(A,B,N,iter_noise,flag)

        D_x = 2;
        D_u = 2;
        [V_a, D_a] = eig(A);
        [V_b, D_b] = eig(B);
        x = zeros(D_x,N);
        u = zeros(D_u,N);
        x_org(:,1) = [1;1];
        x_data = cell(iter_noise+1,1);
        u_data = cell(iter_noise+1,1);
        data = cell(1,2);
%         for i = 1:N
%             u(:,i) = [sin(pi/5*i);cos(i*pi/5)];
%         end
        u = 20*(rand(2,N)-0.5);
        u_data{end} = u;
        for j = 1:N-1
            x_org(:,j+1) = A*x_org(:,j)+B*u(:,j);
        end
        x_data{end} = x_org;
        power_x = norm(x_org.^2);
        x_noise_rand = randn(D_x,N); 
        power_x_arry = ones(D_x,N);
    switch flag
        % x noise influence
        case 1
            for m = 1:iter_noise
                noise_power = ((m-1)*0.0001)*power_x;
                %x_noise = x_noise_rand/norm(x_noise_rand)*noise_power;
                x_noise = ((m-1)*0.0001)*power_x*power_x_arry;
                x(:,1) = [1,1];
                for j = 1:N-1
                    x(:,j+1) = A*(x(:,j)+x_noise(:,j))+B*u(:,j);
                end            
                x_data{m} =x;
            end                
            data{1,1} = x_data;
            data{1,2} = u_data;
        case 2
            % input disturbance
            for m = 1:iter_noise
                noise_power = ((m-1)*0.002)*power_x;
               
                %x_k_noise = x_noise_rand/norm(x_noise_rand)*noise_power;
                x_k_noise = ((m-1)*0.0001)*power_x*power_x_arry;
                x(:,1) = [1,1];
                for j = 1:N-1
                    x(:,j+1) = A*x(:,j)+B*u(:,j) + x_k_noise(:,j);
                end            
                x_data{m} =x;
            end
            data{1,1} = x_data;
            data{1,2} = u_data;
        case 3 
            %u noise influence
            power_u = norm(u.^2);
            u_noise_rand = randn(D_u,N);
            for m = 1:iter_noise
                noise_power_u = ((m-1)*0.01)*power_u;
                u_noise = u_noise_rand/norm(u_noise_rand)*noise_power_u;
                u_data{m} = u_noise+u;
            end
            for i = 1:iter_noise
                u = u_data{i};
                x(:,1) = [1;1];
                for j = 1:N-1
                    x(:,j+1) = A*x(:,j)+B*u(:,j);
                end
                 x_data{i} = x;
            end
            x_data{end} = x_org;
            data{1,1} = x_data;
            data{1,2} = u_data;
            %measurement noise
        case 4 
            for m = 1:iter_noise
                noise_power = ((m-1)*0.002)*power_x;
                %x_noise = x_noise_rand/norm(x_noise_rand)*noise_power;
                x_noise = ((m-1)*0.0001)*power_x*power_x_arry;
                x(:,1) = [1,1];
                for j = 1:N-1
                    x(:,j+1) = A*(x(:,j))+B*u(:,j);
                end   
                x = x_noise+x;
                x_data{m} =x;
            end                
            data{1,1} = x_data;
            data{1,2} = u_data;
    end       
end
       