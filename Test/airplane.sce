clear 

function [H,F,G,Cm,Mm]=predmat(Np,Nc,Q,R,A,B,varargin)
    
    //Np--prediction Horizon
    //Nc--control Horizon
    
    //Test if C matrix is given. if given C=varargin(1) else C=identity(Ns,Ns);
    testcmatrix=argn(2)-6
    //
            //Cm(i,j)=1
    
    if testcmatrix==1 then
        C=varargin(1)
    else
        C=eye(Ns,Ns)
    end
    Ns=size(A,1)
    Nu=size(B,2)
    Ny=size(C,1)
    Ns=Ny //to cater for both state and output cost function
    a=zeros(Ns,Nu)
    Cm=zeros(Np*Ns,Nc*Nu)
    Mm=zeros(Np*Ns,Ns)
    //Creating Cm matrix
    for i=1:Np
        Mm((i-1)*Ns+1:i*Ns,1:Ns)=A^i
        for j=1:Nc
            if i>=j & i<=Nc then
                Cm((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=C*(A^(i-j))*B
            end
            if i>=j & i>Nc then
                if j==Nc then
                    for k=0:(i-Nc)
                        a=a+(C*(A^k)*B)
                    end
                     Cm((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=a
                     a=zeros(Ns,Nu)
                 else
                     Cm((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=C*(A^(i-j))*B
                end
                
            end
            
    end
end

    
   
    H=Cm
    F=Mm
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

Q=eye(4,4)
R=1

//[H,F]=predmat(5,3,Q,R,A,B)
[H,F]=predmat(5,3,Q,R,A,B,C)

