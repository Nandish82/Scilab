function [Acon,bcon,Sxcon]=mpcconstraints(Su,Sx,lbu,ubu,lbx,ubx,Np,Nc)
    % We write the constraints in the form
    % Acon.u=>bcon+Scon*x
    
    
    Nsy=size(Sx,1)/Np%extract size from constraints 
    Ns=size(Sx,2)
    Nu=size(lbu,1)
    
    if Nu==0
        lb=-10000000;
        ub=10000000
        Nu=1
    end
    
        Aucon=[eye(Nc*Nu,Nc*Nu);-eye(Nc*Nu,Nc*Nu)]
        bucon=[repmat(eye(Nu,Nu),Nc,1)*lbu;repmat(eye(Nu,Nu),Nc,1)*-ubu]
        Sucon=zeros(size(bucon,1),Ns)
    
    if isempty(lbx)
        Axcon=[]
        bxcon=[]
        Sxxcon=[]
       
    else
        Axcon=[Su;-Su]
        bxcon=[repmat(eye(Nsy,Nsy),Np,1)*lbx;repmat(eye(Nsy,Nsy),Np,1)*-ubx]
        Sxxcon=[-Sx;Sx]
        
    end
    
    %since in MATLAB  we have to put constraints in the form
    %Au<=b+m*x
    Acon=[Aucon;Axcon]/-1
    bcon=[bucon;bxcon]/-1
    Sxcon=[Sucon;Sxxcon]/-1
    