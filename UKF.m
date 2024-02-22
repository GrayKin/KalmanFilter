function [test_mse, eva_mse]=UKF()
%测试一维匀速直线运动,速度和观测含噪声（其实不一定匀速，近似匀速）
% 状态转移方程
% X(k)=F*X(k-1)+G*u X=[x],F=[1],G=[1],u=[velocity],噪声 velocity=w~N(5,1)
% 观测方程
% Z(k)=H*X(k)+v    Z=[x],H=[1],噪声 v~N(0,1)
N=20;
Xs=zeros(N,1);
v=5;
vel_var=1;
meas_var=1;
R=[1];
Q=[1];
x=0;
P=[1];
n=length(x);
alpha=0.001;
beta=2;
kappa=3-n;
[z,pos]=RobotMove1(x,v,vel_var,meas_var,N);
for i=1:N
    % Generate Sigma Points
    [sigmas_pre,Wm,Wc]=SigmaPoints(n,alpha,beta,kappa,x,P);
    % predict
    [x,P]=UT_pre(sigmas_pre,Wm,Wc,Q,v);
    [sigmas_post,Wm,Wc]=SigmaPoints(n,alpha,beta,kappa,x,P);
    [zp,S]=UT_post(sigmas_post,Wm,Wc,R);
    % Calculate Kalman Gain
    [K]=CalculateK(sigmas_pre,sigmas_post,Wc,x,zp,S);
    % update
    [x,P]=update(x,zp,z(i),K,S,P);
    % result
    Xs(i)=x;
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

function [K]=CalculateK(sigmas_pre,sigmas_post,Wc,x, zp,S)
Pxz =cross_variance(x, zp,sigmas_pre, sigmas_post,Wc);
K = Pxz/S;
end
function [Pxz]=cross_variance(x, z,sigmas_pre, sigmas_post,Wc)
Pxz=zeros(size(sigmas_pre,2),size(sigmas_post,2));
N=size(sigmas_pre,1);
for i=1:N
    dx=sigmas_pre(i,:)-x;
    dz=sigmas_post(i,:)-z;
    Pxz=Pxz+ Wc(i) *dx'*dz;
end
end
function [x,P]=update(x,zp,z,K,S,P)
x=x+K*(z-zp);
P=P-K*S*K';
end
function [x,P]=UT_pre(sigmas,Wm,Wc,Q,v)
% Get through function
sigmas=predict_transit(sigmas,v);
% calculate mean
x=Wm'*sigmas;
% calculate covariance
y=sigmas-x;
P=y'*diag(Wc)*y+Q;
end
function [sigmas]=predict_transit(sigmas,v)
sigmas=sigmas+v;
end
function [z,S]=UT_post(sigmas,Wm,Wc,R)
% Get through function
sigmas=post_transit(sigmas);
% calculate mean
z=Wm'*sigmas;
% calculate covariance
y=sigmas-z;
S=y'*diag(Wc)*y+R;
end
function [sigmas]=post_transit(sigmas)
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