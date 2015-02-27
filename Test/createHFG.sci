
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//function [H,F,G,cm3,Mm,Qbar,RN]=createHFG(N,A,B,Q,P,R)

//the size of Q is the a square matrix of size Ns*Ns..number of states
//the size of P is the same except it is the Terminal weight
//the size of R is usually the same as the above



// ! L.9: mtlb(B) can be replaced by B() or B whether B is an M-file or not.
sizeB = size(mtlb_double(mtlb(B)));
cosizeB = sizeB(2);//number of columns of B i.e the number of inputs
// ! L.11: mtlb(A) can be replaced by A() or A whether A is an M-file or not.
sizeA = size(mtlb_double(mtlb(A)));
sizeA = sizeA(1);

//sizeC=size(C);
//sizeC=sizeC(1); %number of outputs


// ! L.18: mtlb(N) can be replaced by N() or N whether N is an M-file or not.
// ! L.18: mtlb(A) can be replaced by A() or A whether A is an M-file or not.
// ! L.18: mtlb(N) can be replaced by N() or N whether N is an M-file or not.
cm3 = zeros(mtlb_double(mtlb(N))*size(size(mtlb_double(mtlb(A))),"*"),mtlb_double(mtlb(N))*cosizeB);
Mm = [];
// ! L.20: mtlb(N) can be replaced by N() or N whether N is an M-file or not.

for i = mtlb_imp(1,mtlb_double(mtlb(N)))
  // ! L.21: mtlb(A) can be replaced by A() or A whether A is an M-file or not.
  Mm = [Mm;mtlb_double(mtlb(A))^i];
  // ! L.22: mtlb(N) can be replaced by N() or N whether N is an M-file or not.

  for j = mtlb_imp(1,mtlb_double(mtlb(N)))
    if i>=j then
      //cm2((i-1)*ndims(B)+1:(i-1)*2+2,j)=(A^(i-j))*B
      //cm2((i-1)*ndims(A)+1:(i-1)*ndims(A)+1,)=
      // ! L.26: mtlb(A) can be replaced by A() or A whether A is an M-file or not.
      // ! L.26: mtlb(B) can be replaced by B() or B whether B is an M-file or not.
      cm3((i-1)*sizeA+1:i*sizeA,(j-1)*cosizeB+1:j*cosizeB) = (mtlb_double(mtlb(A))^(i-j))*mtlb_double(mtlb(B));
    end;
  end;
end;

Q_N = [];
// ! L.32: mtlb(N) can be replaced by N() or N whether N is an M-file or not.

for i = mtlb_imp(1,mtlb_s(mtlb_double(mtlb(N)),1))
  // ! L.33: mtlb(Q) can be replaced by Q() or Q whether Q is an M-file or not.
  // !! L.33: Matlab function blkdiag not yet converted, original calling sequence used.
  Q_N = 1*mtlb_double(blkdiag(Q_N,mtlb(Q)));
end;

// ! L.36: mtlb(P) can be replaced by P() or P whether P is an M-file or not.
// !! L.36: Matlab function blkdiag not yet converted, original calling sequence used.
Qbar = blkdiag(Q_N,mtlb(P));

RN = []
// ! L.39: mtlb(N) can be replaced by N() or N whether N is an M-file or not.

for i = mtlb_imp(1,mtlb_double(mtlb(N)))
  // ! L.40: mtlb(R) can be replaced by R() or R whether R is an M-file or not.
  // !! L.40: Matlab function blkdiag not yet converted, original calling sequence used.
  RN = blkdiag(RN,mtlb(R));
end;


H = mtlb_a((cm3'*mtlb_double(Qbar))*cm3,mtlb_double(RN));
F = (cm3'*mtlb_double(Qbar))*Mm;
// ! L.46: mtlb(Q) can be replaced by Q() or Q whether Q is an M-file or not.
G = mtlb_a((Mm'*mtlb_double(Qbar))*Mm,mtlb_double(mtlb(Q)));
//end
