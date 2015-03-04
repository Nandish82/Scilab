
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
    // F=Su'*Qtilde*Sx [Nc x Nu]*[Ns ]
    
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

function [Acon,bcon,Sxcon,bxcon,Sxxcon,Axcon,bucon]=mpcconstraints(Su,Sx,lbu,ubu,lbx,ubx,Np,Nc)
    // We write the constraints in the form
    // Acon.u=>bcon+Scon*x
    
    
    Nsy=size(Sx,1) /Np//extract size from constraints 
    Ns=size(Sx,2)
    Nu=size(lbu,1)
    
    disp(Nsy)
    
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
        bxcon=[repmat(eye(Nsy,Nsy),Np,1)*lbx;repmat(eye(Nsy,Nsy),Np,1)*-ubx]
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
