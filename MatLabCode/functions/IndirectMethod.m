function x_indirect = IndirectMethod(N,x_bat,x,u,D_x,D_u)
    x_input1 = [x_bat(:,1:end-1);u(:,1:end-1)];
    x_input2 = x_bat(:,2:end);
    M_con = x_input2*pinv(x_input1);
    A_con = M_con(:,1:D_x);
    B_con = M_con(:,D_x+1:end);
    x_con(:,1) = x(:,1);

    for i=1:N-1
        x_con(:,i+1) = A_con* x_con(:,i) + B_con*u(:,i); 
    end
    x_indirect(:,1) = x(:,1);
   for i=1:N-1
        x_indirect(:,i+1) = A_con* x_indirect(:,i) + B_con*u(:,i); 
   end
end