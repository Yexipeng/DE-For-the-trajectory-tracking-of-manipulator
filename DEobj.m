function Object=object(path1,path2,path3,delta1,delta2,delta3,flag)  %path,delta是2000维
global TE G ts

q1_1=0;q2_1=0;q3_1=0;
tol1_1=0;tol2_1=0;tol3_1=0;
e1_1=0.1;e2_1=0.1;e3_1=0.1;
Tol=[];


tmax=3*TE; %目标函数积分上限为3TE
q1d=1.0;q2d=1.0;q3d=1;
q1op_1=0;dq1op_1=0;
q2op_1=0;dq2op_1=0;
q3op_1=0;dq3op_1=0;


x1_1=2;x2_1=0;
x3_1=0;x4_1=0;
x5_1=0;x6_1=0;
for k=1:1:G 

t(k)=k*ts;

    if t(k)<=TE
       q1op(k)=path1(k); %要逼近的最优轨迹
       dq1op(k)=(q1op(k)-q1op_1)/ts;
       ddq1op(k)=(dq1op(k)-dq1op_1)/ts;
    else
        q1op(k)=q1d;
        dq1op(k)=0;
        ddq1op(k)=0;        
    end    

    if t(k)<=TE
       q2op(k)=path2(k); %要逼近的最优轨迹
       dq2op(k)=(q2op(k)-q2op_1)/ts;
       ddq2op(k)=(dq2op(k)-dq2op_1)/ts;
    else
        q2op(k)=q2d;
        dq2op(k)=0;
        ddq2op(k)=0;        
    end    
    
    if t(k)<=TE
       q3op(k)=path3(k); %要逼近的最优轨迹
       dq3op(k)=(q3op(k)-q3op_1)/ts;
       ddq3op(k)=(dq3op(k)-dq3op_1)/ts;
    else
        q3op(k)=q3d;
        dq3op(k)=0;
        ddq3op(k)=0;        
    end    

%离散模型

g=9.8;
m1=1;m2=0.5;m3=0.2; %质量
L1=1;L2=1.5;L3=0.8;   %连杆长度
r1=0.1;r2=0.2;r3=0.15; %连杆质心长度

I1=4/3*m1^2;        %转动惯量
I2=4/3*m2^2;
I3=4/3*m3^2;

p1=m2*r2^2+m3*L2^2;
p2=m3*r3^2;
p3=m3*r3*L2;

b1=(m2*r2+m3*L2)*g;
b2=m3*r3*g;



M=[I1+p1*cos(x3_1)*cos(x3_1)+p2*cos(x3_1+x5_1)+2*p3*cos(x3_1)*cos(x3_1+x5_1) 0 0;
    0 I2+p1+p2+2*p3*cos(x5_1) p2+p3*cos(x5_1);
    0 p2+p3*cos(x5_1) I3+p2];

C=[-1/2*p1*x4_1*sin(2*x3_1)-1/2*p2*(x4_1+x6_1)*sin(2*x3_1+2*x5_1)-p3*x4_1*sin(2*x3_1+x5_1)-p3*x6_1*cos(x3_1)*sin(x3_1+x5_1) -1/2*p1*x2_1*sin(2*x3_1)-1/2*p2*x2_1*sin(2*x3_1+x5_1)-p3*x2_1*sin(2*x3_1+x5_1) -1/2*p1*x2_1*sin(2*x3_1+x5_1)-p3*x2_1*cos(x3_1)*sin(x3_1+x5_1);
   1/2*p1*x2_1*sin(2*x3_1)+1/2*p2*x2_1*sin(2*x3_1+x5_1)+p3*x2_1*sin(2*x3_1+x5_1) -p3*x6_1*sin(x5_1) -p3*(x4_1+x6_1)*sin(x5_1);
   1/2*p1*x2_1*sin(2*x3_1+x5_1)+p3*x2_1*cos(2*x3_1)*sin(x3_1+x5_1) -p3*x4_1*sin(x5_1) 0];

g=[0;b1*cos(x3_1)+b2*cos(x3_1+x5_1);b2*cos(x3_1+x5_1)];

tol=[tol1_1 tol2_1 tol3_1]';

dq=[x2_1;x4_1;x6_1];

ddq=inv(M)*(tol-C*dq-g);

x2(k)=x2_1+ts*ddq(1);
x1(k)=x1_1+ts*x2(k);

x4(k)=x4_1+ts*ddq(2);
x3(k)=x3_1+ts*x4(k);

x6(k)=x6_1+ts*ddq(3);
x5(k)=x5_1+ts*x6(k);

q1(k)=x1(k); 
dq1(k)=x2(k); 
e1(k)=q1op(k)-q1(k);
de1(k)=(e1(k)-e1_1)/ts;         
    
q2(k)=x3(k); 
dq2(k)=x4(k); 
e2(k)=q2op(k)-q2(k);
de2(k)=(e2(k)-e2_1)/ts;         

q3(k)=x5(k); 
dq3(k)=x6(k); 
e3(k)=q3op(k)-q3(k);
de3(k)=(e3(k)-e3_1)/ts;  

e=[e1(k);e2(k);e3(k)];
de=[de1(k);de2(k);de3(k)];

alpha=30;
epsilon=0.001;

Kp=[50 0 0;0 50 0;0 0 50];
Kd=[150 0 0;0 150 0;0 0 150];


eta=de+alpha*e;
w=alpha*M*de+alpha*C*e;
rho=2.0*tanh(2.6*norm(eta));
v=(eta*rho^2)/(norm(eta)*rho+epsilon);
tol=Kp*eta+w-v+0.006*M*[ddq1op;ddq2op;ddq3op]+0.006*C*[dq1op;dq2op;dq3op];


  
energy(k)=0.33*abs(tol(1)*dq1(k))+0.33*abs(tol(2)*dq2(k))+0.34*abs(tol(3)*dq3(k));
lj(k)=0.3*abs(dq1(k)*0.001)+0.4*abs(dq2(k)*0.001)+0.3*abs(dq3(k)*0.001);
    
    x1_1=x1(k);x2_1=x2(k);
    x3_1=x3(k);x4_1=x4(k);
    x5_1=x5(k);x6_1=x6(k);

    e1_1=e1(k);
    e2_1=e2(k);
    e3_1=e3(k);
    q1op_1=q1op(k);dq1op_1=dq1op(k);
    q2op_1=q2op(k);dq2op_1=dq2op(k);
    q3op_1=q3op(k);dq3op_1=dq3op(k);
    
    tol1_1=tol(1);
    tol2_1=tol(2);
    tol3_1=tol(3);
    
    tol1k(k)=tol(1);
    tol2k(k)=tol(2);  
    tol3k(k)=tol(3); 
    
    Tol=[Tol tol1k(k)];
    
end 
%************计算总能量******************%
energy_all=0;
lj_all=0;
for k=1:1:G
    energy_all=energy_all+0.02*energy(k);
    lj_all=lj_all+lj(k);
end
d1=1*sum(delta1);%参考轨迹的逼近误差
d2=1*sum(delta2);%参考轨迹的逼近误差
d3=1*sum(delta3);%参考轨迹的逼近误差
%********计算目标********%
delta_all=0.3*d1+0.3*d2+0.4*d3;

Object=1*delta_all+0*energy_all;  %used for main.m


if flag==1
    t(1)=0;
    q10=2;q20=0;q30=0;
    for k=1:1:G   %>TE 不包含原点
        t(k)=k*ts;
    if t(k)<=TE
        qr1(k)=(q1d-q10)*(t(k)/TE-1/(2*pi)*sin(4*pi*t(k)/TE))+q10;;   %不含原点的参考轨迹
        qr2(k)=(q2d-q20)*(t(k)/TE-1/(2*pi)*sin(4*pi*t(k)/TE))+q20;   %不含原点的参考轨迹
        qr3(k)=(q3d-q30)*(t(k)/TE-1/(2*pi)*sin(4*pi*t(k)/TE))+q30;
    
    end
  
    end
        figure(1);
        subplot(221);
        
        plot(t,qr1,'b',t,q1op,'r','linewidth',1);
        legend('Ideal trajectory','Optimal trajectory','location','northwest');
        xlabel('Time (s)');ylabel('First Joint angle tracking');
        
        
        subplot(222);
        plot(t,qr2,'b',t,q2op,'r','linewidth',1);
        legend('Ideal trajectory','Optimal trajectory','location','northwest');
        xlabel('Time (s)');ylabel('Second Joint angle tracking');
        
        subplot(223);
        plot(t,qr3,'b',t,q3op,'r','linewidth',1);
        legend('Ideal trajectory','Optimal trajectory','location','northwest');
        xlabel('Time (s)');ylabel('Third Joint angle tracking');
        
        figure(2)
        subplot(221);
        plot(t,dq1op);
        
        subplot(222);
        plot(t,dq2op);
        
        subplot(223);
        plot(t,dq3op);
        
        figure(4)
        subplot(221);
        plot(t,tol1k);
        
        subplot(222);
        plot(t,tol2k);
        
        subplot(223);
        plot(t,tol3k);
        
        figure(5)
        subplot(2,2,1);
        plot(t,q1op,'r',t,q1,'b');
        
        subplot(2,2,2);
        plot(t,q2op,'r',t,q2,'b');
        
        subplot(2,2,3);
        plot(t,q3op,'r',t,q3,'b');
        
        figure(6)
        subplot(2,2,1);
        plot(t,q1op-q1);
        
        subplot(2,2,2);
        plot(t,q2op-q2);
        
        subplot(2,2,3);
        plot(t,q3op-q3);
   
        
end
end
