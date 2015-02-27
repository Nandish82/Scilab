clear 

function [H,F,G,Cm,Mm]=predmat(Np,Nc,Q,R,A,B,varargin)
    
    testcmatrix=argn(2)-6 //test if matrix C has been given or no. 
    if testcmatrix==1 then
        C=varargin(1)
    else
        C=eye(2,2)
    end
    H=C
    F=0
    G=0
    Cm=0
    Mm=0
    
endfunction



Ac=[-1.2822 0 0.98 0;0 0 1 0;-5.4293 0 -1.8366 0;-128.2 128.2 0 0]
Bc=[-0.3;0;-17;0]
Cc=[0 1 0 0;0 0 0 1;-128.2 128.2 0 0]
Dc=[0;0;0;0]

Ns=size(Ac,1)
Nu=size(Bc,2)
Ny=size(Cc,1)

Ts=0.1
A=eye(Ns,Ns)+Ts*Ac
B=Bc*Ts
C=Cc
D=Dc

ref=2
matSS=[A-eye(size(A,1),size(A,1)) B;C(3,:) zeros(size(C(3,:),1),size(B,2))]
refD=[zeros(Ns,1) ;ref]



