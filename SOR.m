% SOR Method

function[x_SOR, comp_time_SOR] = SOR(A,b,n,omega)

x0 = zeros(1,n);
xt = x0;
tic
for t = 1:500
for i = 1:n
    xt(i) = 0;
    x(i) = (1-omega)*x0(i) + omega*(b(i) - sum(A(i,:)*xt'))/A(i,i);
    xt = x0;
    xt(1:i) = x(1:i);
end
comp_time = toc;
    x0 = x;
end
comp_time_SOR = comp_time;
x_SOR = x';
end