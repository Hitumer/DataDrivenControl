function cost = CostFunction(x_reference,x)
         [~,N] = size(x);
        n = floor(N/4);
        Q = [1,0,;0,1];
        n = n-1;
%        cost = sum(sum((x - x_reference)'*Q*(x - x_reference)))/(m*n);
         power_x = sum(sum(x_reference(:,1:n).^2));
         error = (x(:,1:n)-x_reference(:,1:n)).^2;
         error = error'*Q;
         cost = sum(sum(error))/(power_x);
end