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
Ts=0.5//sampling time
temp=(eye(Ns,Ns)+0.5*Ts*Ac)*inv((eye(Ns,Ns)-0.5*Ts*Ac))
temp=inv(eye(Ns,Ns)-Ac*Ts)
//A=temp
A=eye(Ns,Ns)+Ts*Ac 
B=Bc*Ts
C=Cc
D=Dc
//[A,B,C,D]=mpcdiscretize(Ac,Bc,Cc,Dc,Ts)
//steady state matrix
Cref=[0,0,0,1];
Ssmat=[A-eye(size(A,1),size(A,1)) B;Cref zeros(size(Cref,1),size(B,2))]

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
altitude_rate_min=-30
altitude_rate_max=30


//Constraints on input rate
lbu=[elevator_slew_min]
ubu=[elevator_slew_max]

//constraints on input
umax=elevator_angle_max
umin=elevator_angle_min
//Constraints on states
lbx=[-100000;pitch_angle_min;-10000;-10000;umin]
ubx=[100000;pitch_angle_max;10000;10000;umax]



///for delta formulation we add the input as a state in the output to be able to add constraints
lby=[pitch_angle_min;altitude_min;altitude_rate_min;elevator_angle_min]
uby=[pitch_angle_max;altitude_max;altitude_rate_max;elevator_angle_max]

//Constraints matrices
// A.x>=b


[Acon,bcon,Sxcon,bxcon,Sxxcon,Axcon,bucon]=mpcconstraints(Su,Sx,lbu,ubu,lby,uby,Np,Nc)

//

//Simulation

sim_time=1.5
dt=Ts
time_vec=[0:dt:sim_time]
Npoints=length(time_vec)
xdata=zeros(size(Adelta,1),Npoints)
ydata=zeros(size(Cdelta,1),Npoints)
udata=zeros(size(Bdelta,2),Npoints)
x0=[0;0.0;0;0.0;0]

xref=[0;0.0;0;40;0]
xdata(:,1)=x0
ydata(:,1)=Cdelta*xdata(:,1)



//bucon=[repmat(eye(Nu,Nu),Nc,1)*lbu;repmat(eye(Nu,Nu),Nc,1)*-ubu]


for i=1:Npoints-1
    //soln=qp_solve(H,F*(xdata(:,i)-xref),Acon',bcon+Sxcon*(xdata(:,i)-xref),0)
   // [soln,iact,iter,f]=qp_solve(H,(F*(xdata(:,i)-xref))',Acon',bcon+Sxcon*(xdata(:,i)-xref),0)
    cxdata=xdata(:,i)-xref
    //if(pmodulo(i*dt,Ts)==0)
   soln=qld(H,F*cxdata,-1*Axcon,-1*(bxcon+Sxxcon*cxdata),bucon(1:Nc),-1*bucon(Nc+1:$),0)
    //soln=qp_solve(H,F*cxdata,Acon',bcon+Sxcon*cxdata,0)
    udata(:,i)=soln(1)
    //udata(:,i)=1;
    //else
    //    udata(:,i+1)=udata(:,i);
    //end
    xdata(:,i+1)=Adelta*xdata(:,i)+Bdelta*(udata(:,i))
    ydata(:,i+1)=Cdelta*xdata(:,i+1)
end

scf(3)
clf(3)
subplot(321)
plot(time_vec,xdata(1,:)')
subplot(322)
plot(time_vec,xdata(2,:)'*(180/3.142))
subplot(323)
plot(time_vec,xdata(3,:)')
subplot(324)
plot(time_vec,xdata(4,:)')
subplot(325)
plot2d2(time_vec,xdata(5,:)'*(180/3.142))
subplot(326)
plot2d2(time_vec,udata(1,:)')


scf(4)
clf(4)
subplot(221)
plot(time_vec,ydata(1,:)'*(180/3.142))
plot([0 sim_time],[pitch_angle_max pitch_angle_max]*(180/3.142),'r')
plot([0 sim_time],[pitch_angle_min pitch_angle_min]*(180/3.142),'g')
subplot(222)
plot(time_vec,ydata(2,:)')
subplot(223)
plot(time_vec,ydata(3,:)')
plot([0 sim_time],[altitude_rate_max altitude_rate_max],'r')
plot([0 sim_time],[altitude_rate_min altitude_rate_min],'r')
subplot(224)
plot2d2(time_vec,xdata(5,:)'*(180/3.142))
plot([0 sim_time],[elevator_angle_max elevator_angle_max]*(180/3.142),'r')
plot([0 sim_time],[elevator_angle_min elevator_angle_min]*(180/3.142),'r')

