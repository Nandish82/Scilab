function [H,F,G,Su,Sx]=mpcpredmat(Np,Nc,Q,R,A,B,varargin)
%Np--prediction Horizon
    %Nc--control Horizon
    
    %Test if C matrix is given. if given C=varargin(1) else C=identity(Ns,Ns);
    
    %matrix sizes
    % Sux=(Np x Ns)*(Nc x Nu)
    % Suy=(Np x Ny) *(Nc*Nu)
    % Sxx=(Np*Ns) * (Ns)
    % Sxy=(Np*Ny) *(Ns)
    % Qtildex=(Np*Ns) * (Np*Ns)
    % Qtildey=(Np*Ny) * (Np*Ny)
    % Rtildex=(Nc*Nu) * (Nc*Nu)
    % Rtildex=Rtildey
    % H=Cm'*Qtildexy*Cm+R
    % F=
    nargin
    testcmatrix=nargin-6
    %
            %Cm(i,j)=1
    %varargin is a CELL in matlab to extract the elements of the first cell
    %use { }
    if testcmatrix==1 
        C=varargin{1,:}
    else
        C=eye(size(A,1),size(A,1))
    end
    Ns=size(A,1);
    Nu=size(B,2);
    Ny=size(C,1);
    Ns=Ny %to cater for both state and output cost function;
    a=zeros(Ns,Nu);
    Su=zeros(Np*Ns,Nc*Nu);
    Sx=zeros(Np*Ns,Ns);
    Qtilde=zeros(Np*Ns,Np*Ns);
    Rtilde=zeros(Nc*Nu,Nc*Nu);
    %Creating Cm matrix and Mm matrix
    for i=1:Np
        Sx((i-1)*Ns+1:i*Ns,1:size(C,2))=C*(A^i);
        for j=1:Nc
            if (i>=j & i<=Nc)
                Su((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=C*(A^(i-j))*B;
            end
                
            if (i>=j & i>Nc)
                if j==Nc
                    for k=0:(i-Nc)
                        a=a+(C*(A^k)*B);
                    end
                     Su((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=a;
                     a=zeros(Ns,Nu);
                 else
                     Su((i-1)*Ns+1:i*Ns,(j-1)*Nu+1:j*Nu)=C*(A^(i-j))*B;
                end
                
            end
            
    end
end

%Qtilde and Rtilde Matrix

P=Q ;%for now. P should be the solution to the algebraic ricatti equation

for i=1:Np
    if i<=Np-1 
        Qtilde((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=Q;
    else
       Qtilde((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns)=P;
    end
    
end

for i=1:Nc
        Rtilde((i-1)*Nu+1:i*Nu,(i-1)*Nu+1:i*Nu)=R;
    
end
    
   
    H=Su'*Qtilde*Su+Rtilde;
    F=Su'*Qtilde*Sx;
    G=Sx'*Qtilde*Sx+C'*Q*C;
    
    