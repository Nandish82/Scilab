xdel(winsid());
clear ;
clc;

exec('mpcfunctions.sci'); //loading functions in mpcfunctions.sci

[A,B,C,D,x0]=Readjac("ServoMech_Model_.jac0")


//size of matrices
Ns=size(A,1)
Nu=size(B,2)
Ny=size(C,1)

//discrete model
Ts=0.1//sampling time
Ad=eye(Ns,Ns)+Ts*A
Bd=B*Ts
Cd=C
Dd=D

//Delta formulation
[Adelta,Bdelta,Cdelta,Ddelta]=mpcdelta(Ad,Bd,Cd,Dd)

Qx=eye(size(A,1),size(A,1)) //weights on state
Ri=eye(size(B,2),size(B,2))  //weights on input
Rd=eye(size(B,2),size(B,2))  //weights on delta

Qy=eye(size(C,1),size(C,1)) //weights on output

Np=10; //prediction
Nc=3;  // control 

///constraints

//input
input_max=1000
input_min=-1000

//state
state_max=[1000 1000 1000 1000 1000]
state_min=[-1000 -1000 -1000 -1000 -1000]

//outpu
output_max=[1000 1000]
output_min=[-1000 -1000]

//delta
delta_max=[1000]
delta_min=[-1000]


f_type=1; //0=NORMAL  1=DELTA
p_type=0; //0=state 1=output

//// p_type | f_type | s_type
////   0        0        00   -State and Normal-0
///    0        1        01   -State and Delta-1
///    1        0        10   -Ouput and Normal-2
////   1        1        11   -Output and Delta-3

s_type=f_type+(p_type*2) //sum of types

Np=10 //predicition horizon
Nc=3 //control horizon

select (s_type)
case (0) then
    Q=Qx
    R=Ri
    [H,F,G,Su,Sx]=predmat(Np,Nc,Q,R,Ad,Bd) 
    lbu=[input_min]
    ubu=[input_max]
    lbA=[state_min]
    ubA=[state_max]
case(1) then
    Q=sysdiag(Qx,Ri)
    R=Rd
    [H,F,G,Su,Sx]=predmat(Np,Nc,Q,R,Adelta,Bdelta) 
    lbu=[delta_min]
    ubu=[delta_max]
    lbA=[state_min input_min]
    ubA=[state_max input_max]
case(2) then
    Q=Qy
    R=Ri
    [H,F,G,Su,Sx]=predmat(Np,Nc,Q,R,Ad,Bd,Cd)
    lbu=[input_min]
    ubu=[input_max]
    lbA=[output_min]
    ubA=[output_max]
case(3) then
    Q=sysdiag(Qy,Ri)
    R=Rd   
    [H,F,G,Su,Sx]=predmat(Np,Nc,Q,R,Adelta,Bdelta,Cdelta)
    lbu=[delta_min]
    ubu=[delta_max]
    lbA=[output_min input_min]
    ubA=[output_max input_max]
end

mpc.Q=Q
mpc.R=R

if f_type==0 then
    mpc.A=Ad
    mpc.B=Bd
    mpc.C=Cd
    mpc.D=Dd
else
    mpc.A=Adelta
    mpc.B=Bdelta
    mpc.C=Cdelta
    mpc.D=Ddelta
end

//constraints matrices

[Acon,bcon,Sxcon,bxcon,Sxxcon,Axcon,bucon]=mpcconstraints(Su,Sx,lbu',ubu' ,lbA' ,ubA',Np,Nc)






    


