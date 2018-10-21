xl=0; xr=1;
J = 100;
dx = (xr-xl) / J; 
tf = 0.1;
Nt = 10000;
dt = tf/Nt;

mu = dt/(dx)^2;

if mu > 0.5
    error('mu should < 0.5!');
end

% Evaluate the initial conditions
x = xl : dx : xr; % generate the grid point % f(1:J+1) since array index starts from 1
%f = sin(pi*x) + sin(2*pi*x);
Cb = 0.1;
f = Cb * ones(length(x));
f_psi = zeros(length(x));
% store the solution at all grid points for all time steps 
Cp = zeros(J+1,Nt);
Cm = zeros(J+1,Nt);
Psi = zeros(J+1,Nt);
k = 0.1;
A = full(gallery('tridiag',length(x),1,-2,1));
Ainv = inv(A);

for n = 1:Nt
    t = n*dt; % current time
    % boundary condition at left side
    %Neuman Boundary conditions
    Psi0 = -11;
    if n==1 % first time step
        for j=2:J % interior nodes
            Cp(j,n) = f(j) + mu*(f(j+1)-2*f(j)+f(j-1))-(mu/dx^2)*(f_psi(j+1)-2*f_psi(j)+f_psi(j-1))*(f(j+1)-f(j))-(mu/dx^2)*(f_psi(j+1)-f_psi(j))*(f(j+1)-f(j)); 
            Cm(j,n) = f(j) + mu*(f(j+1)-2*f(j)+f(j-1))+(mu/dx^2)*(f_psi(j+1)-2*f_psi(j)+f_psi(j-1))*(f(j+1)-f(j))+(mu/dx^2)*(f_psi(j+1)-f_psi(j))*(f(j+1)-f(j)); 
            sum_eta = 0;
            for eta=2:J
                 sum_eta = sum_eta - k * Ainv(j,eta) * f_psi(eta);
            end
            Psi(j,n) = sum_eta;
        end
        
        Cp(1,n) = Cp(2,n);
        Cp(J+1,n) = Cb;
        Cm(1,n) = Cm(2,n);
        Cm(J+1,n) = Cb;
        Psi(1,n) = Psi0;
        Psi(J+1,n) = 0;
    else
        for j=2:J %
            % the left-end point % the right-end point interior nodes
            Cp(j,n)=Cp(j,n-1)+mu*(Cp(j+1,n-1)-2*Cp(j,n-1)+Cp(j-1,n-1)) - (mu/dx^2) * (Cp(j+1,n-1)-Cp(j,n-1)*(Psi(j+1,n-1)-Psi(j,n-1))) - (mu/dx^2) * (Cp(j,n-1)*(Psi(j+1,n-1)-2*Psi(j,n-1)+Psi(j-1,n-1)));
            Cm(j,n)=Cm(j,n-1)+mu*(Cm(j+1,n-1)-2*Cm(j,n-1)+Cm(j-1,n-1))+ (mu/dx^2) * (Cm(j+1,n-1)-Cm(j,n-1)*(Psi(j+1,n-1)-Psi(j,n-1))) + (mu/dx^2) * (Cm(j,n-1)*(Psi(j+1,n-1)-2*Psi(j,n-1)+Psi(j-1,n-1)));
            sum_eta = 0;
            for eta=2:J
                 sum_eta = sum_eta - k * Ainv(j,eta) * f_psi(eta);
            end
            Psi(j,n) = sum_eta;
        end
        
        Cp(1,n) = Cp(2,n);
        Cp(J+1,n) = Cb;
        Cm(1,n) = Cm(2,n);
        Cm(J+1,n) = Cb;
        Psi(1,n) = Psi0;
        Psi(J+1,n) = 0;
    end
    
end



% Plot the results
tt = dt : dt : Nt*dt;
figure(1)
colormap(gray);
surf(x,tt, Cp');
xlabel('x')
ylabel('t')
zlabel('u')
title('Cp')
figure(2)
surf(x,tt, Cm'); % 3-D surface plot xlabel('x')
ylabel('t')
zlabel('u')
title('Cm')

figure(3)
surf(x,tt, Psi'); % 3-D surface plot xlabel('x')
xlabel('x')
ylabel('t')
zlabel('u')
title('Psi')