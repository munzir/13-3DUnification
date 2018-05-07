A=csvread('../../build/initial_traj.csv');
x=A(:,1); psii=A(:,2); th=A(:,3);
dx=A(:,4); dpsi=A(:,5); dth=A(:,6);
X0=A(:,7); Y0=A(:,8);
ddth=A(:,9); tau0=A(:,10);
dt=0.01;
t=0:dt:dt*(length(A)-1);

figure; 
subplot(2,2,1); plot(t,x,t,dx,t,ddth); title('$x$, $\dot{x}$ and $\ddot{\theta}$', 'Interpreter', 'latex'); legend({'$x$', '$\dot{x}$', '$\ddot{\theta}$'}, 'Interpreter', 'latex'); grid on
subplot(2,2,2); plot(t,th,t,dth,t,ddth); title('$\theta$, $\dot{\theta}$ and $\ddot{\theta}$', 'Interpreter', 'latex'); legend({'$\theta$', '$\dot\theta$', '$\ddot\theta$'}, 'Interpreter', 'latex'); grid on
subplot(2,2,3); plot(t,psii,t,dpsi,t,tau0); title('$\psi$, $\dot{\psi}$ and $\tau_0$', 'Interpreter', 'latex'); legend({'$\psi$', '$\dot\psi$', '$\tau_0$'}, 'Interpreter', 'latex'); grid on
subplot(2,2,4); plot(X0,Y0); xlabel('X_0'); ylabel('Y_0'); grid on
