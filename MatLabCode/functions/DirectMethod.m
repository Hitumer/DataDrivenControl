function x_predict = DirectMethod(N,x_bat,x,u,D_x,D_u,L) 
    x_predict = zeros(2,N-1);
    power_x = norm(x.^2);
    for j = 1:(N-1)/(L-1)
    
        if j ==1
            for i = 1:N
                 x_hankel1(:,i) =x_bat(:,i);
                 u_hankel1(:,i) = u(:,i);
                 H_x = hankelmatrix(x_hankel1,L);
                 H_u = hankelmatrix(u_hankel1,L);
                 state_init = zeros(1,(L*D_u+D_x));
                 for n = 1:L
                    state_init(1,n*2-1) = u(1,n);
                    state_init(1,n*2) = u(2,n);
                 end
            end
                 state_init(1,end-D_x+1:end) = x(:,1);
            
                 inves = pinv([H_u;H_x(1:D_x,:)]);
                 alpha = state_init*inves';
    
                 x_state = H_x(D_x+1:end,:)*alpha';
                 x_predict(:,1) = x(:,1);
                 x_predict(:,2:L) = reshape(x_state,D_x,length(x_state)/D_x);
        else
    
            for i = (L-1)*(j-1)+1:N
                x_hankel(:,i-(L-1)*(j-1)) =x_bat(:,i);
                u_hankel(:,i-(L-1)*(j-1)) = u(:,i);
            end
                H_x = hankelmatrix(x_hankel,L);
                H_u = hankelmatrix(u_hankel,L);
                state_init = zeros(1,(L*D_u+D_x));
                for k = 1:L
                    state_init(1,k*D_u-1) = u(1,k+(j-1)*(L-1));
                    state_init(1,k*D_u) = u(2,k+(j-1)*(L-1));
                end
                state_init(1,end-D_x+1:end) = x_predict(:,(L-1)*(j-1)+1);
        
        inves = pinv([H_u;H_x(1:D_x,:)]);
        alpha = state_init*inves';
    
        x_state = H_x(D_x+1:end,:)*alpha';
        x_predict(:,(L-1)*(j-1)+2:(L-1)*j+1) = reshape(x_state,D_x,length(x_state)/D_x);
        x_hankel = zeros(D_x,N - (L-1)*(j));
        u_hankel = zeros(D_u,N - (L-1)*(j));
        end      
    end
end