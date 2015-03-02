xdel(winsid());
clear ;
clc;

exec('mpcfunctions.sci'); //loading functions in mpcfunctions.sci

//Airplane Model From Maciejowski [pg 71]

//state space Model
Ac=[-1.2822 0 0.98 0;0 0 1 0;-5.4293 0 -1.8366 0;-128.2 128.2 0 0]
Bc=[-0.3;0;-17;0]
Cc=[0 1 0 0;0 0 0 1;-128.2 128.2 0 0]
Dc=[0;0;0]

//size of matrices
Ns=size(Ac,1)
Nu=size(Bc,2)
Ny=size(Cc,1)

//discrete model
Ts=0.5 //sampling time
A=eye(Ns,Ns)+Ts*Ac
B=Bc*Ts
C=Cc
D=Dc

//Delta formulation
[Adelta,Bdelta,Cdelta,Ddelta]=mpcdelta(A,B,C,D)

///setting of mpc problem
Np=10; //prediction horizon
Nc=3; //control horizon

Qx=eye(size(Adelta,1),size(Adelta,1))
Rx=eye(size(Bdelta,2),size(Bdelta,2))

Qy=eye(size(Cdelta,1),size(Cdelta,1))
Ry=eye(size(Bdelta,2),size(Bdelta,2))
[H,F,G,Su,Sx]=predmat(Np,Nc,Qy,Ry,Adelta,Bdelta,Cdelta)

//states=[xxx pitch angle xxx altitude]
//outputs=[pitch angle altitude altitude rate]

//Constraints
elevator_angle_max=0.262
elevator_angle_min=-0.262

elevator_slew_max=0.524
elevator_slew_min=-0.524

pitch_angle_max=0.349
pitch_angle_min=-0.349

//Constraintes on output
altitude_min=-1000000
altitude_max=1000000
altitude_rate_min=-100000
altitude_rate_max=100000


//Constraints on input rate
lbu=[elevator_slew_min]
ubu=[elevator_slew_max]

//constraints on input
umax=elevator_angle_max
umin=elevator_angle_min
//Constraints on states
lbx=[-100000;pitch_angle_min;-10000;-10000;umin]
ubx=[100000;pitch_angle_max;10000;10000;umax]




lby=[pitch_angle_min;altitude_min;altitude_rate_min]
uby=[pitch_angle_max;altitude_max;altitude_rate_max]

//Constraints matrices
// A.x>=b


[Acon,bcon,Sxcon]=mpcconstraints(Su,Sx,lbu,ubu,lby,uby,Np,Nc)

//

//Simulation

sim_time=10
dt=0.5
time_vec=[0:dt:sim_time]
Npoints=length(time_vec)
xdata=zeros(size(Adelta,1),Npoints)
ydata=zeros(size(Cdelta,1),Npoints)
udata=zeros(size(Bdelta,2),Npoints)
x0=[0;0.0;0;0.5;0]

xref=[0;0;0;40;0]
xdata(:,1)=x0
ydata(:,1)=C*xdata(:,1)

for i=1:Npoints-1
    //soln=qp_solve(H,F*(xdata(:,i)-xref),Acon',bcon+Sxcon*(xdata(:,i)-xref),0)
    soln=qp_solve(H,F*(ydata(:,i)),Acon',bcon,0)
    u(i)=soln(1)
    
    xdata(:,i+1)=Adelta*xdata(:,i)+Bdelta*(u(i))
    ydata(:,i+1)=Cdelta*ydata(:,i)
end

scf(3)
clf(3)
subplot(321)
plot(time_vec,xdata(1,:)')
subplot(322)
plot(time_vec,xdata(2,:)')
subplot(323)
plot(time_vec,xdata(3,:)')
subplot(324)
plot(time_vec,xdata(4,:)')
subplot(325)
plot(time_vec,xdata(5,:)')
