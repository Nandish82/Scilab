
function [Adelta,Bdelta,Cdelta,Ddelta]=mpcdelta(A,B,C,D)
    Ns=size(A,1);
    Nu=size(B,2);
    Ny=size(C,1);
    Adelta=[A B;zeros(Nu,Ns) eye(Nu,Nu)];
    Bdelta=[B;eye(Nu,Nu)];
    Cdelta=[C D];
    Ddelta=[D];