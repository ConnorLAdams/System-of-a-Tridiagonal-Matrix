% Gauss-elimination without pivoting

function[x_Gauss_Wo_Pivot] = Gauss_WO_Pivot(H,b_solve,L,U)

% Forward substitution - Ly = b
b = b_solve;
[m,n] = size(H);
y(1,1)=b(1)/L(1,1);


for i = 2:m
    y(i,1) = (b(i)- L(i,i-1:-1:1)*y(i-1:-1:1))/L(i,i);
    
end


% Backward substitution - Ux = y
b = y; % letting Ux = y = a new b.
x = zeros(n,1);
x(n,1)=b(n)/U(n,n);
for i = n-1:-1:1
    x(i,1) = (b(i)- U(i,i+1:n)*x(i+1:n))/U(i,i);
end
x_Gauss_Wo_Pivot = x;

end