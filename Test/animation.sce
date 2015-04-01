
// animation_point.sce

clear; xdel(winsid());

// Create data
t = 0:0.005:1;    // Time data
x = sin(2*%pi*t); // Position data

//Model Parameters
m=0.5
A=[0 1;0 0]
B=[0;1/m]
C=[1 0]
D=0
Ts=0.01

 s1=syslin('c',A,B,C,D)  //creates continous LTI system
 sld=dscr(s1,Ts)         //discretize the system
 //Assign them to matrices.
 Ad=sld.A
 Bd=sld.B
 Cd=sld.C
 Dd=sld.D
 
 sim_time=5
 
time_vec=0:Ts:sim_time
xdata=length
 

// Draw initial figure
figure(1);
plot(x(1),0,'o');
h_compound = gce();
h_compound.children.mark_size = 20;
h_compound.children.mark_background = 2;
h_axes = gca();
h_axes.data_bounds = [-1.5,-1.5;1.5,1.5];

// Animation Loop
i = 1;
while i<=length(x)
    drawlater();
    h_compound.children.data = [t(i),x(i)];
    //h_compound.children.data = [t(1:i),x(1:i)];
    drawnow();
    i = i+1;
end
