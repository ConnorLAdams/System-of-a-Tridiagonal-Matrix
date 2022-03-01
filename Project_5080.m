% Project MAT 5080

clear
close all
clc

% Parameters
n = 4;
a = sym('a',[1,n-1]);
b = sym('b',[1,n]);
c = sym('c',[1,n-1]);
d = sym('d',[1,n]);
c_over_d = c(1:n-1)./d(1:n-1);

% Matrices setup
T = diag(b) + diag(a,1) + diag(c,-1);
[L,U] = LUFactorization(a,d,c_over_d,n);


% Showing L*U = T
LUverification(L,U,T,a,b,c,d,n)

%Notice that T_new is the exact same as T so this means L and U are the LU
%factorizations of T.


% Part 2 ii
n = 5;

H = diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);


% Part 3
    % Solve for eigenvalues
eig_H = round(double((eig(H))),5);

    % Solve for the cosine equation
for k = 1:n
    eig_cosine(k) = round(2-2*cos(k*pi/(n+1)),5);
end
    % Checking to see if the 2 arrays are equivalent
eigCheck(eig_H,eig_cosine,n)

spectral_radius1 = max(abs(double(eig(H))));

    % Set up for problem 4
n = 500;
x_exact = zeros(n,n);
x_Gauss_Seidel = x_exact; x_Gauss_Wo_Pivot = x_exact;
x_Gauss_WP_Pivot = x_exact; x_Jacobi = x_exact; x_LUFactorization = x_exact;
x_SOR_1 = x_exact; x_SOR_1_half = x_exact; x_SOR_half = x_exact;
% New Parameters for 4
for n = 1:500
H = diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);

clear a b c d
a = -1*ones(1,n-1);
c = a;
b = 2*ones(1,n);
d(1) = b(1);
for i = 2:n
    d(i) = b(i) - (a(i-1)*c(i-1)/d(i-1));
end
c_over_d = c(1:n-1)./d(1:n-1);

% Calculating a new L and U
[L,U] = LUFactorization(a,d,c_over_d,n);

format long
  
b_solve = zeros(n,1);
b_solve(1) = 1; b_solve(n) = 1;

    % Exact solutions
    x_exact(1:n,n) = ones(n,1);

    % Part 4a - Gaussian Elimination Without Pivoting
        % n comes form part 2ii so n=4  

tic
x_Gauss_Wo_Pivot(1:n,n) = Gauss_WO_Pivot(H,b_solve,L,U);
comp_time_Wo(n) = toc;

tic
x_Gauss_WP_Pivot(1:n,n) = Gauss_WP_Pivot(H,b_solve,n);
comp_time_WP(n) = toc;
    % Part 4b LU Factorization method page 5 of week 4

tic
x_LUFactorization(1:n,n) = LUFactorization_Solution(L,U,b_solve);
comp_time_LU(n) = toc;
    % Part 4c - Jacobi
        % For fun make a while loop to see how long it takes to get to the
        % correct points.

tic
x_Jacobi(1:n,n) = Jacobi(H,b_solve,n);
comp_time_J(n) = toc;

    % Part 4d - Gauss Seidel

tic    
x_Gauss_Seidel(1:n,n) = Gauss_Seidel(H,b_solve,n);
comp_time_GS(n) = toc;

    % Part 4e - SOR
     i = 1;       
for omega = .5:.5:1.5
    [x_SOR(:,i), comp_time_SOR(i)] = SOR(H,b_solve,n,omega);
    i = i + 1;
end
x_SOR_half(1:n,n) = x_SOR(:,1);
x_SOR_1(1:n,n) = x_SOR(:,2);
x_SOR_1_half(1:n,n) = x_SOR(:,3);
comp_time_SOR_half(n) = comp_time_SOR(1);
comp_time_SOR_1(n) = comp_time_SOR(2);
comp_time_SOR_1_half(n) = comp_time_SOR(3);
clear x_SOR
end

% Plot the CPU time of each method
    iterations = 1:n;
    sample_size = round((n - 1)/2);
    
    % finding slope and intercepts of each method
cpu_times = [comp_time_GS; comp_time_J; comp_time_LU; comp_time_SOR_1;...
    comp_time_SOR_1_half; comp_time_SOR_half; comp_time_Wo; comp_time_WP];
for i = 1:8
[slope(i),intercept(i)] = slope_CPU(cpu_times(i,:),iterations,sample_size,n);
end
close all
color = [0,0,1; 0,0,0; 1,0,0; 0,1,0; 1,.5,0; 0,1,1; 1,0,1; .25,.75,.5];
for i = 1:8
    figure(1);
loglog(1:500,cpu_times(i,:),'color',color(i,:),'Linewidth',1.25)
hold on
    figure(2);
loglog(1:500,intercept(i)*iterations.^(slope(i)),'--','color',color(i,:),'Linewidth',1.25)
hold on
end
figure(1)
legend('Gauss-Seidel','Jacobi','LU','SOR, omega = 1','SOR, omega = 1.5',...
    'SOR, omega = .5','Gauss-Elimination w/o pivot',...
    'Gauss-Elimination w/ partial pivot','Location','SouthEast')
title('CPU Time Versus H Matrix Size')
xlabel('Matrix Size (nxn)')
ylabel('CPU Time (sec)')
figure(2)
legend('line fit for Gauss-Seidel','line fit for Jacobi','line fit for LU',...
    'line fit for SOR, omega = 1','line fit for SOR, omega = 1.5',...
    'line fit for SOR, omega = .5',...
    'line fit for Gauss-Elimination w/o pivot',...
    'line fit for Gauss-Elimination w/ partial pivot','Location',...
    'SouthEast')
xlabel('Matrix Size (nxn)')
ylabel('Average CPU Time (sec)')
title('Average CPU Time Versus H Matrix Size')
hold off


% plot the spectral radius for each matrix size.

% Setting up for B_SOR
for n = 1:200
H = diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);
D_H = diag(diag(H));
I = eye(n);

Lower = tril(H,-1);
Upper = triu(H,1);
        
i = 1;

for omega = 0:.01:2
    B_SOR = inv(I-omega*inv(D_H)*Lower) * ((1-omega)*I + omega*inv(D_H)*Upper);
    spectral_radius(i) = max(abs(eig(B_SOR)));
    i = 1+i;
end
omeg = find(spectral_radius==min(spectral_radius));
w(n) = (omeg-1)/100;
end

figure(3)
plot(1:n,w,'b')
hold on
yline(2,'--r','Horizontal Asymptote')
hold off
ylim([0,2.2])
xlabel('Matrix Size (nxn)')
ylabel('Relaxation Parameter (w)')
legend('Optimal Relaxation Parameter (w)','Location','East')
title('Relaxation Parameter Versus Matrix Size for H')

figure(4)
plot(0:.01:2,spectral_radius)
xlabel('Relaxation Parameter (w)')
ylabel('Spectral Radius')
title('Spectral Radius Versus Relaxation Parameter')
legend('Spectral Radius','Location','West')


% Time to compute the error.

for i = 1:500
   norm_exact(i) = norm(x_exact(:,i));
   norm_GS(i) = norm(x_Gauss_Seidel(:,i));
   norm_GWOP(i) = norm(x_Gauss_Wo_Pivot(:,i));
   norm_GWPP(i) = norm(x_Gauss_WP_Pivot(:,i));
   norm_J(i) = norm(x_Jacobi(:,i));
   norm_LU(i) = norm(x_LUFactorization(:,i));
   norm_SOR_1(i) = norm(x_SOR_1(:,i));
   norm_SOR_1_half(i) = norm(x_SOR_1_half(:,i));
   norm_SOR_half(i) = norm(x_SOR_half(:,i));
end
norms = [norm_exact;norm_GS;norm_GWOP;norm_GWPP;norm_J;norm_LU;norm_SOR_1;...
    norm_SOR_1_half;norm_SOR_half];

% Errors

for i = 1:9
    errors(i,:) = norms(1,:) - norms(i,:);
end

color = [0,0,1; 0,0,0; 1,0,0; 0,1,0; 1,.5,0; 0,1,1; 1,0,1; .25,.75,.5];
figure(5)
for i = 2:9
   plot(1:500,errors(i,:),'color',color(i-1,:),'LineWidth',1.25)
   hold on
end
xlabel('Matrix Size (n)')
ylabel('Errors')
legend('error_{GS}','error_{GWOP}','error_{GWPP}','error_{J}','error_{LU}','error_{SOR 1}',...
    'error_{SOR 1.5}','error_{SOR .5}','Location','SouthEast')
title('The errors For Finding The Solutions of HX=b Versus Matrix Size of H')
grid on
