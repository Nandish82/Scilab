clear all
close all

%Servo model
A=[-2.00000000005751e+001    -1.04719755089633e+000     0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000;
    1.90985931710274e+002    -2.00000000000000e-001    -1.06683333333333e+000     0.00000000000000e+000     2.13366666666666e+001;
    0.00000000000000e+000     6.00000000000001e+000     0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000;
    0.00000000000000e+000     0.00000000000000e+000     1.06683333333333e+000    -2.50000000000000e+000    -2.13366666666666e+001;
    0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000     6.00000000000001e+000     0.00000000000000e+000]
B=[
    1.00000000050712e+000;
    0.00000000000000e+000;
    0.00000000000000e+000;
    0.00000000000000e+000;
    0.00000000000000e+000;]
C=[0.00000000000000e+000     0.00000000000000e+000     1.11718525420157e+000     0.00000000000000e+000    -2.23437050840314e+001;
    0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000     1.00000000000000e+000]
    
D=[0.00000000000000e+000;
    0.00000000000000e+000;
   ]

C1=[0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000     1.00000000000000e+000;
    0.00000000000000e+000     0.00000000000000e+000     0.00000000000000e+000     1.00000000000000e+000     0.00000000000000e+000]

B1=[1.00000000050712e+000     0.00000000000000e+000;
    0.00000000000000e+000     0.00000000000000e+000;
    0.00000000000000e+000     0.00000000000000e+000;
    0.00000000000000e+000     -9.54929658551370e-001;
    0.00000000000000e+000     0.00000000000000e+000;]
D1=[    0.00000000000000e+000     0.00000000000000e+000;
    0.00000000000000e+000     0.00000000000000e+000;]


C=C1
%Discretization Approximate
Ts=0.1
Ak=(eye(size(A,1),size(A,1))+Ts*A)
Bk=Ts*B1
Bd=Bk(:,2)
Bk=Bk(:,1)
Ck=C1
Dk=D1

%Delta formulation
Gc=ss(A,B(:,1),C,D)
Gd=c2d(Gc,Ts,'zoh')
Ak=Gd.a
Bk=Gd.b
Ck=Gd.c
Dk=Gd.d
[Adelta,Bdelta,Cdelta,Ddelta]=mpcdelta(Ak,Bk(:,1),Ck,Dk(:,1))


%%observer with disturbance model


%%observer with delta formulation

poles_obs=eig(Adelta)/10

L=place(Adelta',Cdelta',poles_obs)
L=L' 

%%MPC Control Setup
Qy=eye(size(Ck,1),size(Ck,1))
Qx=eye(size(Ak,1),size(Ak,1))
R=1*eye(size(Bk,2),size(Bk,2))
Rrate=eye(size(Bk,2),size(Bk,2))

Qydelta=Qy
Qxdelta=blkdiag(Qx,R)
Rdelta=Rrate

Np=10
Nc=10
[H,F,G,Su,Sx]=mpcpredmat(Np,Nc,Qx,R,Ak,Bk)

%%constraints
INPUT_RATE_MIN=-1000000
INPUT_RATE_MAX=1000000


INPUT_VOLTAGE_MIN=-8000
INPUT_VOLTAGE_MAX=8000

%%state constraints
CURRENT_MIN=-100000
CURRENT_MAX=100000

POSITION_1_MIN=-1000000
POSITION_1_MAX=1000000

SPEED_1_MIN=-1000000
SPEED_1_MAX=1000000


POSITION_2_MIN=-99999900
POSITION_2_MAX=999999

SPEED_2_MIN=-10000
SPEED_2_MAX=100000

%%
lbrate=[INPUT_RATE_MIN];
ubrate=[INPUT_RATE_MAX]
lbu=[INPUT_VOLTAGE_MIN];
ubu=[INPUT_VOLTAGE_MAX]
%%state constraints
lbx=[CURRENT_MIN SPEED_1_MIN  POSITION_1_MIN  SPEED_2_MIN POSITION_2_MIN]'
ubx=[CURRENT_MAX SPEED_1_MAX  POSITION_1_MAX  SPEED_2_MAX POSITION_2_MAX]'
%%output constraints
lby=[SPEED_2_MIN POSITION_2_MIN]'
uby=[SPEED_2_MAX POSITION_2_MAX]'

%%constraints calculation
[Acon,bcon,Sxcon]=mpcconstraints(Su,Sx,lbu,ubu,lbx,ubx,Np,Nc)

%%Simulation
sim_time=20
dt=Ts
time_vec=[0:dt:sim_time]
N=length(time_vec)
xdata_p=zeros(size(Ak,1),N);
xdata_e=zeros(size(Adelta,1),N);
ydata_e=zeros(size(Cdelta,1),N);

udata=zeros(size(Bk,2),N);
udata_rate=zeros(size(Bk,2),N);
xdata_p(:,1)=[0,0,0,0,0]' 

xdata_e(:,1)=[0,0,0,0,0,0]'
ydata_e(:,1)=Cdelta*[0,0,0,0,0,0]'

ydata_p(:,1)=Ck*xdata_p(:,1)


u=100

xdata_p(:,1)=[0;0;0;0;0]
xdata_e(:,1)=[0;0;0;0;0;0]
xref=[0;0;0;0;0;0]
u=0;
td=0
Ssmat=[Ak-eye(5,5) Bk;0 0 0 0 1 0]
ss=inv(Ssmat)*[0 0 0 0 0 10]'
xss=ss(1:5)
uss=ss(6)
for i=1:N-1
    cxdata=xdata_p(:,i)-xss
   %%if(pmodulo(i*dt,Ts)==0)
    %%soln=qld(H,F*cxdata,-1*Axcon,-1*(bxcon+Sxxcon*cxdata),1*bucon(1:Nc),-1*bucon(Nc+1:$),0)
    %%soln=qpsolve(H,F*cxdata,-1*Axcon,(bxcon+Sxxcon*cxdata),1*bucon(1:Nc),-1*bucon(Nc+1:$),0)
    soln=quadprog(H,F*cxdata,Acon,bcon+Sxcon*cxdata);
    u=soln(1)+uss;
    udata(:,i)=u;
    
    xdata_p(:,i+1)=Ak*xdata_p(:,i)+Bk*(udata(:,i))%+Bd*td;
    ydata_p(:,i+1)=Ck*xdata_p(:,i);
    %%xdata_e(:,i+1)=Adelta*xdata_e(:,i)+Bdelta*deltaU(i)+L*(ydata_p(:,i)-Cdelta*xdata_e(:,i))
    %xdata_e(:,i+1)=Adelta*xdata_e(:,i)+Bdelta*soln(1)
    %ydata_e(:,i+1)=Cdelta*xdata_e(:,i+1)
end


subplot(221)
plot(time_vec,ydata_p(1,:)')
subplot(222)
plot(time_vec,ydata_p(2,:)')
subplot(223)
plot(time_vec,udata(1,:)')
subplot(224)
%plot(time_vec,xdata(5,:)'*(180/3.142))
%plot(time_vec,xdata_e(6,:))
