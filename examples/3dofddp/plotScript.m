Aref=csvread('../../build/initial_traj.csv');
A=csvread('../../build/mpc_traj.csv');
Aref = Aref(1:length(A), :);

xref=Aref(:,1); psiiref=Aref(:,2); thref=Aref(:,3);
dxref=Aref(:,4); dpsiref=Aref(:,5); dthref=Aref(:,6);
X0ref=Aref(:,7); Y0ref=Aref(:,8);
ddthref=Aref(:,9); tau0ref=Aref(:,10);
x=A(:,1); psii=A(:,2); th=A(:,3);
dx=A(:,4); dpsi=A(:,5); dth=A(:,6);
X0=A(:,7); Y0=A(:,8);
ddth=A(:,9); tau0=A(:,10);

dt=0.01;
t=0:dt:dt*(length(A)-1);

figure; 
subplot(3,3,1); plot(t,x,t,xref); title('$x$ and ${x}_{ref}$', 'Interpreter', 'latex'); legend({'$x$', '${x}_{ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,2); plot(t,dx,t,dxref); title('$\dot{x}$ and $\dot{x}_{ref}$', 'Interpreter', 'latex'); legend({'$\dot{x}$', '$\dot{x}_{ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,3); plot(t,psii,t,psiiref); title('$\psi$ and $\psi_{ref}$', 'Interpreter', 'latex'); legend({'$\psi$', '$\psi_{ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,4); plot(t,dpsi,t,dpsiref); title('$\dot\psi$ and $\dot\psi_{ref}$', 'Interpreter', 'latex'); legend({'$\dot\psi$', '$\dot\psi_{ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,5); plot(t,th,t,thref); title('$\theta$ and $\theta_{ref}$', 'Interpreter', 'latex'); legend({'$\theta$', '$\theta_{ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,6); plot(t,tau0,t,tau0ref); title('$\tau_0$ and $\tau_{0ref}$', 'Interpreter', 'latex'); legend({'$\tau_0$', '$\tau_{0ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,7); plot(t,ddth,t,ddthref); title('$\ddot\theta$ and $\ddot\theta_{ref}$', 'Interpreter', 'latex'); legend({'$\ddot\theta$', '$\ddot\theta_{ref}$'}, 'Interpreter', 'latex'); grid on
subplot(3,3,8); plot(X0,Y0,X0ref,Y0ref); legend({'state', 'ref'}); grid on

