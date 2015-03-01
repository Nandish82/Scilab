xdel(winsid());
clear ;
clc;


function [H,F,G,Su,Sx]=predmat(Np,Nc,Q,R,A,B,varargin)
    
    //Np--prediction Horizon
    //Nc--control Horizon
    
    //Test if C matrix is given. if given C=varargin(1) else C=identity(Ns,Ns);
    
    //matrix sizes
    // Sux=(Np x Ns)*(Nc x Nu)
    // Suy=(Np x Ny) *(Nc*Nu)
    // Sxx=(Np*Ns) * (Ns)
    // Sxy=(Np*Ny) *(Ns)
    // Qtildex=(Np*Ns) * (Np*Ns)
    // Qtildey=(Np*Ny) * (Np*Ny)
    // Rtildex=(Nc*Nu) * (Nc*Nu)
    // Rtildex=Rtildey
    // H=Cm'*Qtildexy*Cm+R
    // F=
    
    testcmatrix=argn(2)-6
    //
            //Cm(i,j)=1
    
    if testcmatrix==1 then
        C=varargin(1)
    else
        C=eye(size(A,1),size(A,1))
    end
    Ns=size(A,1)
    Nu=size(B,2)
    Ny=size(C,1)
    Ns=Ny //to cater for both state and output cost function
    a=zeros(Ns,Nu)
    Su=zeros(Np*Ns,Nc*Nu)
    Sx=zeros(Np*Ns,Ns)
    Qtilde=zeros(Np*Ns,Np*Ns)
    Rtilde=zeros(Nc*Nu,Nc*Nu)
    //Creating Cm matrix and Mm matrix
    for i=1:Np
        Sx((i-1)*Ns+1:i*Ns,1:size(C,2))=C*(A^i)
        for j=1:Nc
            if i>=j & i<=Nc then
                Su((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=C*(A^(i-j))*B
            end
                
            if i>=j & i>Nc then
                if j==Nc then
                    for k=0:(i-Nc)
                        a=a+(C*(A^k)*B)
                    end
                     Su((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=a
                     a=zeros(Ns,Nu)
                 else
                     Su((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=C*(A^(i-j))*B
                end
                
            end
            
    end
end

//Qtilde and Rtilde Matrix

P=Q //for now. P should be the solution to the algebraic ricatti equation

for i=1:Np
    if i<=Np-1 then
        Qtilde((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=Q
    else
       Qtilde((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=P
    end
    
end

for i=1:Nc
        Rtilde((i-1)*Nu+1:i*Nu,(i-1)*Nu+1:i*Nu)=R
    
end
    
   
    H=Su'*Qtilde*Su+Rtilde
    F=Su'*Qtilde*Sx
    G=Sx'*Qtilde*Sx+C'*Q*C
    
endfunction

function [Adelta,Bdelta,Cdelta,Ddelta]=mpcdelta(A,B,C,D)
    Ns=size(A,1)
    Nu=size(B,2)
    Ny=size(C,1)
    Adelta=[A B;zeros(Nu,Ns) eye(Nu,Nu)]
    Bdelta=[B;eye(Nu,Nu)]
    Cdelta=[C D]
    Ddelta=[D]
endfunction

function [Acon,bcon,Sxcon]=mpcconstraints(Su,Sx,lbu,ubu,lbx,ubx,Np,Nc)
    // We write the constraints in the form
    // Acon.u=>bcon+Scon*x
    
    
    Ns=size(Sx,1) /Np//extract size from constraints 
    Nu=size(lbu,1)
    
    if Nu==0 then
        lb=-10000000;
        ub=10000000
        Nu=1
    end
    
        Aucon=[eye(Nc*Nu,Nc*Nu);-eye(Nc*Nu,Nc*Nu)]
        bucon=[repmat(eye(Nu,Nu),Nc,1)*lbu;repmat(eye(Nu,Nu),Nc,1)*-ubu]
        Sucon=zeros(size(bucon,1),Ns)
    
    if lbx~=[] then
        Axcon=[Su;-Su]
        bxcon=[repmat(eye(Ns,Ns),Np,1)*lbx;repmat(eye(Ns,Ns),Np,1)*-ubx]
        Sxxcon=[-Sx;Sx]
    else
        Axcon=[]
        bxcon=[]
        Sxxcon=[]
        
    end
    
    //since in scilab we have to put constraints in the form
    Acon=[Aucon;Axcon]
    bcon=[bucon;bxcon]
    Sxcon=[Sucon;Sxxcon]
    
endfunction





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
   [Acon,bcon,Sxcon]=mpcconstraints(Su,Sx,lb,ub,lbx-xss,ubx-xss,Np,Nc)
    soln=qp_solve(H,F*(xdata(:,i)-xss),Acon',bcon+Sxcon*(xdata(:,i)-xss),0)
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
