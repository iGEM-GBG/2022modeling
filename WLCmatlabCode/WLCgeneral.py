# Generated with SMOP  0.41-beta
# from libsmop import *
# WLCgeneral.m
Compound='prot'
# WLCgeneral.m:3

disp('The chosen compound is ' + Compound)
N=50
# WLCgeneral.m:6

aT=10 ** 4
# WLCgeneral.m:7

if Compound == 'prot':
    L=10
# WLCgeneral.m:10
    lp=0.6
# WLCgeneral.m:11
else:
    if Compound == 'DNA':
        L=10.2
# WLCgeneral.m:13
        lp=40
# WLCgeneral.m:14
    else:
        L=1.75
# WLCgeneral.m:16
        lp=0.01
# WLCgeneral.m:17

R=zeros(1,aT)
# WLCgeneral.m:20

Rx=zeros(1,N)
# WLCgeneral.m:21

Ry=zeros(1,N)
# WLCgeneral.m:22

Rz=zeros(1,N)
# WLCgeneral.m:23

Rx_sum=zeros(1,aT)
# WLCgeneral.m:24

Ry_sum=zeros(1,aT)
# WLCgeneral.m:25

Rz_sum=zeros(1,aT)
# WLCgeneral.m:26


sigma=sqrt(L / (dot(N,lp)))
# WLCgeneral.m:29

#Create probability distribution for theta (vikta efter omkrets för given vinkel i sfär):
N_Dist=1000
# WLCgeneral.m:32
lb=0
# WLCgeneral.m:33
ub=copy(pi)
# WLCgeneral.m:34
Angle=linspace(lb,ub,N_Dist)
# WLCgeneral.m:35
probdist_angles_core=sin(Angle) / 2
# WLCgeneral.m:36
cdf_angles=cumsum(probdist_angles_core)
# WLCgeneral.m:37
cdf_angles=cdf_angles / max(cdf_angles)
# WLCgeneral.m:38

cdf_angles[1]=0
# WLCgeneral.m:39
cdf_angles[end()]=1
# WLCgeneral.m:40
ProbDist_theta=makedist('PiecewiseLinear','x',Angle,'Fx',cdf_angles)
# WLCgeneral.m:41
#     rw=random(ProbDist_angles,[1000000,1]);
#     histogram(rw)

#Skapa sannolikhetsfördelning för rho (Normalfördelning, men inom -pi:pi):
N_Dist=1000
# WLCgeneral.m:46
lb=- pi
# WLCgeneral.m:47
ub=copy(pi)
# WLCgeneral.m:48
Angle=linspace(lb,ub,N_Dist)
# WLCgeneral.m:49
probdist_angles_core=exp(- Angle ** 2 / (dot(2,sigma ** 2)))
# WLCgeneral.m:50

cdf_angles=cumsum(probdist_angles_core)
# WLCgeneral.m:51
cdf_angles=cdf_angles / max(cdf_angles)
# WLCgeneral.m:52

cdf_angles[1]=0
# WLCgeneral.m:53
cdf_angles[end()]=1
# WLCgeneral.m:54
ProbDist_rho=makedist('PiecewiseLinear','x',Angle,'Fx',cdf_angles)
# WLCgeneral.m:55
#     rw=random(ProbDist_angles,[10000000,1]);
#     histogram(rw)

for j in arange(1,aT).reshape(-1):
    # Create initial segment with random angles
#     theta_0=random(ProbDist_theta,1); # 2*pi*(randn(1)-1/2);
#     phi_0 = 2*pi*(randn(1)-1/2);
    # Create initial segment with given angle
    theta_0=1
# WLCgeneral.m:66
    phi_0=1
# WLCgeneral.m:67
    Rx[1]=multiply(dot(L / N,cos(phi_0)),sin(theta_0))
# WLCgeneral.m:70
    Ry[1]=multiply(dot(L / N,sin(phi_0)),sin(theta_0))
# WLCgeneral.m:71
    Rz[1]=dot(L / N,cos(theta_0))
# WLCgeneral.m:72
    rho=concat([normrnd(0,sigma,concat([1,N]))])
# WLCgeneral.m:74
    theta=dot(dot(2,pi),(randn(1,N) - 1 / 2))
# WLCgeneral.m:76
    for i in arange(1,N - 1).reshape(-1):
        # Copy previous vector and normalize
        u=dot(concat([[Rx(i)],[Ry(i)],[Rz(i)]]),N) / L
# WLCgeneral.m:81
        # See https://en.wikipedia.org/wiki/Rodrigues#27_rotation_formula
        # Vector n, orthogonal to u, n = (1,-a/b,0)/norm((1,-a/b,0))
        n=concat([1,- u(1) / u(2),0]) / norm(concat([1,- u(1) / u(2),0]))
# WLCgeneral.m:87
        W=concat([0,- n(3),n(2),n(3),0,- n(1),- n(2),n(1),0])
# WLCgeneral.m:90
        RotMat_rho=eye(3) + dot(sin(rho(i)),W) + dot((1 - cos(rho(i))),W ** 2)
# WLCgeneral.m:95
        u_bent=dot(RotMat_rho,u)
# WLCgeneral.m:96
        W=concat([0,- u(3),u(2),u(3),0,- u(1),- u(2),u(1),0])
# WLCgeneral.m:99
        RotMat_theta=eye(3) + dot(sin(theta(i)),W) + multiply((1 - cos(theta(i))),W ** 2)
# WLCgeneral.m:103
        u_new=dot(RotMat_theta,u_bent)
# WLCgeneral.m:104
        Rx[i + 1]=dot(L / N,u_new(1))
# WLCgeneral.m:107
        Ry[i + 1]=dot(L / N,u_new(2))
# WLCgeneral.m:108
        Rz[i + 1]=dot(L / N,u_new(3))
# WLCgeneral.m:109
    # Plot endpoints of all generated chains
# Let aT be 10^2 for short running time
#     figure(1)
#     scatter3(Rx,Ry,Rz,'k.')  # Plot endpoints   
#     axis equal; hold on
#     scatter3(0,0,0,'r') # Plot origo as a red ring
    Rx_sum[j]=sum(Rx)
# WLCgeneral.m:121
    Ry_sum[j]=sum(Ry)
# WLCgeneral.m:122
    Rz_sum[j]=sum(Rz)
# WLCgeneral.m:123
    R[j]=sqrt(Rx_sum(j) ** 2 + Ry_sum(j) ** 2 + Rz_sum(j) ** 2)
# WLCgeneral.m:125

# Plot the chain generated last
figure(2)
# cumsum is used for adding elements cumulatively for the coordinates of
# all segment chain joints
Rx_cumsum=cumsum(Rx)
# WLCgeneral.m:132

Ry_cumsum=cumsum(Ry)
# WLCgeneral.m:133
Rz_cumsum=cumsum(Rz)
# WLCgeneral.m:134
plot3(Rx_cumsum,Ry_cumsum,Rz_cumsum,'k')
axis('image')
xlabel('X')
ylabel('Y')
zlabel('Z')
# Create histogram
# ti=char(input("Title? (Title + N) "));
ti='Estimated probability density'
# WLCgeneral.m:145
figure(3)
lb=0
# WLCgeneral.m:147

ub=copy(L)
# WLCgeneral.m:148

step=dot((ub - lb),5) / (sqrt(aT))
# WLCgeneral.m:149

y,edges=histcounts(R,arange(lb,ub,step),nargout=2)
# WLCgeneral.m:150

x=arange(lb,(ub - step),step)
# WLCgeneral.m:151

k=1 / trapz(x,y)
# WLCgeneral.m:152

plot(x,dot(y,k),'b.')
xlabel('Chain reach [$nm$]','interpreter','latex','FontSize',14)
ylabel('Probability density [1/$nm$]','interpreter','latex','FontSize',14)
eval(sprintf('title('%s %s','interpreter','latex','FontSize',14)',ti,num2str(N)))