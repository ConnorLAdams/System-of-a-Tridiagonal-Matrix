% Jacobi Method

function[x_Jacobi] = Jacobi(A,b,n)

x0 = zeros(1,n);
for t = 1:500
for i = 1:n
    xt = x0;
    xt(i) = 0;
    x(i) = (b(i)-sum(A(i,:)*xt'))/A(i,i);
end
    x0 = x;
end
x_Jacobi = x';
end