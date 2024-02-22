function [eva_mse,test_mse]=KF()
%测试一维匀速直线运动,速度和观测含噪声（其实不一定匀速，近似匀速）
% 状态转移方程
% X(k)=F*X(k-1)+G*u X=[x],F=[1],G=[1],u=[velocity],噪声 velocity=w~N(5,1)
% 观测方程
% Z(k)=H*X(k)+v    Z=[x],H=[1],噪声 v~N(0,1)
N=20;
x0=0;
v=5;
vel_var=1;
meas_var=1;
[z,pos]=RobotMove1(x0,v,vel_var,meas_var,N);
Xa=[0];%初始状态
Pa=[1];%过程噪声的协方差，基于经验，可以一开始设得很大，而后收敛
Q=[1];%控制矩阵的协方差
H=[1];%观测矩阵
R=[1];%观测噪声的协方差矩阵
u=[v];
Xs=zeros(N,1);
for i=1:N
    %Predict
    [Xp,Pp]=predict(Xa,Pa,u,Q);
    %Calculate K
    K=calculateK(Pp,H,R);
    %Update
    [Xa,Pa]=update(z(i),Xp,Pp,K,H,R);
    %result
    Xs(i)=Xa;
end
test_error = Xs(:)-pos(:,1);
test_mse =  sum(test_error.^2) / length(test_error );
eva_error = z-pos(:,1);
eva_mse =  sum(eva_error.^2) / length(eva_error );
plot(Xs);
hold on
plot(z(:),'o','MarkerSize',3);
plot(pos(:,1));
title("Mse:",test_mse);
xlabel('Time/(s)');ylabel('Position/(m)');
legend("Filter","Measurement","Robot");
hold off
end



function [Xa,Pa]=update(z,Xp,Pp,K,H,R)
I=eye(1);
Xa=Xp+K*(z-H*Xp);
Pa = (I-K*H)*Pp*(I-K*H)' + K*R*K';
end
function K=calculateK(P,H,R)
K=P'*H'/(H*P*H'+R);
end
function [Xp,Pp]=predict(Xa,Pa,u,Q)
%Calculate F G
F=[1];
G=[1];
%1 Write down state transition equation
Xp=F*Xa+G*u;
%2 inference Covariance
Pp=F*Pa*F'+G*Q*G';
end
function [robot,position]=RobotMove1(x0,v,vel_var,meas_var,N)
%Gaussian 噪声
%initial
x=x0;
robot=zeros(N,1);
position=zeros(N,2);
vel = v;
dt=1;
vel_std=sqrt(vel_var);
for i=1:N
    %move
    velocity=vel+ randn()*vel_std;
    if(velocity<0)
        velocity=0;
    end
   
    dx=velocity;
    x=x+ dx * dt;
    % sense
    measurement(1) = x +randn()*sqrt(meas_var);
    %result
    robot(i,1)=measurement(1);
    position(i,1)=x;
    position(i,2)=velocity;
end
end

function [robot,position]=RobotMove2(x0,v,vel_var,meas_var,N)
%Levy 噪声
%initial
x=x0;
robot=zeros(N,1);
position=zeros(N,2);
vel = v;
dt=1;
vel_std=sqrt(vel_var);
for i=1:N
    %move
    velocity=vel+ randn()*vel_std;
   
   
    dx=velocity;
    x=x+ dx * dt;
    % sense
    measurement(1) = x +stblrnd(1.5,0,sqrt((meas_var)/2),0);
    %result
    robot(i,1)=measurement(1);
    position(i,1)=x;
    position(i,2)=velocity;
end
end
function [robot,position]=RobotMove3(x0,v,vel_var,meas_var,N)
%均匀噪声
%initial
x=x0;
robot=zeros(N,1);
position=zeros(N,2);
vel = v;
dt=1;
vel_std=sqrt(vel_var);
for i=1:N
    %move
    velocity=vel+ randn()*vel_std;
   
   
    dx=velocity;
    x=x+ dx * dt;
    % sense
    measurement(1) = x +rand(sqrt(12*meas_var))-sqrt(3*meas_var);
    %result
    robot(i,1)=measurement(1);
    position(i,1)=x;
    position(i,2)=velocity;
end
end