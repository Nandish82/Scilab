xdel(winsid());
clear ;
clc;



exec('mpcfunctions.sci'); //loading functions in mpcfunctions.sci
Ac=[-1.2822 0 0.98 0;0 0 1 0;-5.4293 0 -1.8366 0;-128.2 128.2 0 0]
Bc=[-0.3;0;-17;0]
Cc=[0 1 0 0;0 0 0 1;-128.2 128.2 0 0]
Dc=[0;0;0]

Ns=size(Ac,1)
Nu=size(Bc,2)
Ny=size(Cc,1)

Ts=0.25
A=eye(Ns,Ns)+Ts*Ac
B=Bc*Ts
C=Cc
D=Dc

SSmat=[A-eye(Ns,Ns) B;C(2,:) zeros(1,Nu)]
r=40
ss=inv(SSmat)*[zeros(Ns,1);r]
xss=ss(1:4)
uss=ss(5)


Np=4
Nc=4
Qx=eye(4,4)
R=10


// Constraints Matrices
//C*x>=b
//Number of inputs =1
// lb<u<ub
lb=-0.262-uss
ub=0.262-uss
bcon1=[lb*ones(Nc,1);-ub*ones(Nc,1)]
//[-1;1]*u>=[-ub;lb]
Cons=[eye(Nc,Nc);-eye(Nc,Nc)]
//qp_solve(Q, p, C, b, me)
//me=number of equality constraints



[H,F,G,Su,Sx]=predmat(Np,Nc,Qx,R,A,B)


///Simulation simple euler
sim_time=20
dt=0.25
N=int(sim_time/dt)
time_vec=[0:dt:(N-1)*dt]

[Acon,bcon,Sxcon,bxcon,Sxxcon,Axcon,bucon]=mpcconstraints(Su,Sx,lbu,ubu,lbx,ubx,Np,Nc)


x0=[0;0.0;0;0]

xdata=zeros(size(A,1),N)
xdata(:,1)=x0
u=zeros(1,length(time_vec))
//lqr controller
lbx=[-100000 -10000 -100000 -100000]';
ubx=[100000 10000 100000 100000]';

Big=sysdiag(Qx,R);    //Now we calculate C1 and D12
[w,wp]=fullrf(Big);C1=wp(:,1:4);D12=wp(:,5:$);   //[C1,D12]'*[C1,D12]=Big
P=syslin('d',A,B,C1,D12);    //The plant (continuous-time)
[klqr,xlqr]=lqr(P)

for i=1:N-1
   [cxdata=xdata(:,i)-xref
    soln=qld(H,F*cxdata,-1*Axcon,-1*(bxcon+Sxxcon*cxdata),bucon(1:Nc),-1*bucon(Nc+1:$),0)
    //u(i)=klqr*xdata(:,i)//soln(1)
    u(i)=soln(1)
    xdata(:,i+1)=A*xdata(:,i)+B*(u(i)+uss)
end
scf(1)
clf(1)
subplot(221)
plot(time_vec,xdata(1,:)')
subplot(222)
plot(time_vec,xdata(2,:)')
subplot(223)
plot(time_vec,xdata(3,:)')
subplot(224)
plot(time_vec,xdata(4,:)')
scf(2)
clf(2)
plot(time_vec,u)

//[Adelta,Bdelta,Cdelta,Ddelta]=mpcdelta(A,B,C,D)
//
//Qx=eye(5,5)
//R=10
//
//lbu=-0.524
//ubu=0.524
//
//lbx=[-100000;-0.349;-100000;-100000;-0.28]
//
//ubx=[100000;0.349;100000;100000;0.28]
//
//Np=20
//Nc=20
//
//[H,F,G,Su,Sx]=predmat(Np,Nc,Qx,R,Adelta,Bdelta)
//
//sim_time=4
//dt=0.25
//N=int(sim_time/dt)
//time_vec=[0:dt:(N-1)*dt]
//
//
//
//x0=[0;0.0;0;10;0]
//
//xdata=zeros(size(Adelta,1),N)
//xdata(:,1)=x0
//u=zeros(1,length(time_vec))
//
//for i=1:N-1
//    [Acon,bcon,Sxcon]=mpcconstraints(Su,Sx,lbu,ubu,lbx,ubx,Np,Nc)
//    soln=qp_solve(H,F*xdata(:,i),Acon',bcon+Sxcon*xdata(:,i),0)
//    //soln=-inv(H)*F*xdata(:,i)
//    u(i)=soln(1)
//    xdata(:,i+1)=Adelta*xdata(:,i)+Bdelta*(u(i))
//end
//
//
//scf(3)
//clf(3)
//subplot(321)
//plot(time_vec,xdata(1,:)')
//subplot(322)
//plot(time_vec,xdata(2,:)')
//subplot(323)
//plot(time_vec,xdata(3,:)')
//subplot(324)
//plot(time_vec,xdata(4,:)')
//subplot(325)
//plot(time_vec,xdata(5,:)')
//
//scf(4)
//clf(4)
//plot(time_vec,u,'o')
