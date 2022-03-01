% Gauss-Seidel Method

function[x_Gauss_Seidel] = Gauss_Seidel(A,b,n)

x0 = zeros(1,n);
for j = 1:500
for i = 1:n
    xt(i) = 0;
    x(i) = (b(i)-sum(A(i,:)*xt'))/A(i,i);
    xt = x0;
    xt(1:i) = x(1:i);
end
x0 = x;
end
x_Gauss_Seidel = x';

end