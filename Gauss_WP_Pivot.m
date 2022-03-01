% Gauss with partial pivoting

function[x_Gauss_WP_Pivot] = Gauss_WP_Pivot(H,b,n)

if n==1
        x_Gauss_WP_Pivot = .5;
else
    % Creating the augmented matrix
Augmen = [H,b];
    % Partial Pivoting has Begun
for i = 1:n-1
    [Max, Row] = max(abs(Augmen(i:n,i)));
    store = Augmen(i,:);
    Augmen(i,:) = Augmen(Row+i-1,:);
    Augmen(Row+i-1,:) = store;
        % Time for echlon form
    for j = i+1:n
       multiplier = Augmen(j,i)/Augmen(i,i);
       Augmen(j,:) = Augmen(j,:) - multiplier * Augmen(i,:);
    end
    
end
for i = 2:n
    multiplier = Augmen(1:i-1,i)/Augmen(i,i);
    Augmen(1:i-1,:) = Augmen(1:i-1,:) - multiplier * Augmen(i,:);
end
x_Gauss_WP_Pivot = Augmen(:,n+1)./diag(Augmen);
end
end