# Generated with SMOP  0.41-beta
#from libsmop import *
# TEVcompletionProb.m

    # Main file for TEV completion
# Calculate TEV completion probability by generating many paths, 
# DNA is modeled as a straight line
#    clf
#    clc
#    clear
#    close_('all')
    # dCas9 diameter. 70 å = 7 nm
import math
from turtle import dot
import numpy as np
from numpy import *


dCasLength=7
# TEVcompletionProb.m:7

N=50
# TEVcompletionProb.m:9

aT=dot(10 ** 6,2)
# TEVcompletionProb.m:10

L=11.5
# TEVcompletionProb.m:11

bp=30
# TEVcompletionProb.m:12

DL=dot(0.34,bp)
# TEVcompletionProb.m:13

rotPerbp=0.626573
# TEVcompletionProb.m:14

lp=0.6
# TEVcompletionProb.m:15

R=zeros(1,aT)
# TEVcompletionProb.m:17

Rx=zeros(1,N)
# TEVcompletionProb.m:18

Ry=zeros(1,N)
# TEVcompletionProb.m:19

Rz=zeros(1,N)
# TEVcompletionProb.m:20

Rx_sum=zeros(1,aT)
# TEVcompletionProb.m:21

Ry_sum=zeros(1,aT)
# TEVcompletionProb.m:22

Rz_sum=zeros(1,aT)
# TEVcompletionProb.m:23

Position=bp.zeros(aT,3)
# TEVcompletionProb.m:24

DNAvector=bp.concat([DL,0,0])
# TEVcompletionProb.m:26


sigma=np.sqrt(L / (dot(N,lp)))
# TEVcompletionProb.m:28

""" #Create probability distribution for theta (vikta efter omkrets för given vinkel i sfär):
N_Dist=1000
# TEVcompletionProb.m:31
lb=0
# TEVcompletionProb.m:32
ub=math.pi
# TEVcompletionProb.m:33
Angle=np.linspace(lb,ub,N_Dist)
# TEVcompletionProb.m:34
probdist_angles_core=np.sin(Angle) / 2
# TEVcompletionProb.m:35
cdf_angles=np.cumsum(probdist_angles_core)
# TEVcompletionProb.m:36
cdf_angles=cdf_angles / max(cdf_angles)
# TEVcompletionProb.m:37

cdf_angles[1]=0
# TEVcompletionProb.m:38
cdf_angles[-1]=1
# TEVcompletionProb.m:39
ProbDist_theta=makedist('PiecewiseLinear','x',Angle,'Fx',cdf_angles)
 """

# TEVcompletionProb.m:40
RxVis=zeros(2,N)
# TEVcompletionProb.m:42
RyVis=zeros(2,N)
# TEVcompletionProb.m:43
RzVis=zeros(2,N)
# TEVcompletionProb.m:44
RVis=zeros(2,N,3)
# TEVcompletionProb.m:46
nvis=0
# TEVcompletionProb.m:48
pi = math.pi
for j in arange(1,aT).reshape(-1):
    # Create initial segment with random angles
#     theta_0=random(ProbDist_theta,1); # 2*pi*(randn(1)-1/2);
#     phi_0=2*pi*(randn(1)-1/2);
    # Create initial segment with given angles
    if j <= aT / 2:
        theta_0=pi / 2 - 0.01
# TEVcompletionProb.m:56
        phi_0=0 + 0.01
# TEVcompletionProb.m:57
    else:
        theta_0=pi / 2 - 0.01
# TEVcompletionProb.m:59
        phi_0=pi - 0.01
# TEVcompletionProb.m:60
    # (Rx, Ry, Rz) är en vektor med längd som pekar från ett segment till nästa segment. Rx representerar ej segmentens koordinat i x!
    Rx[1]=multiply(dot(L / N,cos(phi_0)),sin(theta_0))
# TEVcompletionProb.m:64
    Ry[1]=multiply(dot(L / N,sin(phi_0)),sin(theta_0))
# TEVcompletionProb.m:65
    Rz[1]=dot(L / N,cos(theta_0))
# TEVcompletionProb.m:66
    rho=np.concat([np.normrnd(0,sigma,np.concat([1,N]))])
# TEVcompletionProb.m:68
    theta=dot(dot(2,pi),(np.randn(1,N) - 1 / 2))
# TEVcompletionProb.m:70
    for i in arange(1,N - 1).reshape(-1):
        # Copy previous vector u=(a,b,c)
        u=dot(np.concat([[Rx(i)],[Ry(i)],[Rz(i)]]),N) / L
# TEVcompletionProb.m:75
        # See https://en.wikipedia.org/wiki/Rodrigues#27_rotation_formula
        # Vector n, orthogonal to u, n = (1,-a/b,0)/norm((1,-a/b,0))
        n=np.concat([1,- u(1) / u(2),0]) / np.norm(np.concat([1,- u(1) / u(2),0]))
# TEVcompletionProb.m:81
        # section
        W=np.concat([0,- n(3),n(2),n(3),0,- n(1),- n(2),n(1),0])
# TEVcompletionProb.m:85
        RotMat_rho=np.eye(3) + dot(np.sin(rho(i)),W) + dot((1 - np.cos(rho(i))),W ** 2)
# TEVcompletionProb.m:90
        u_bent=dot(RotMat_rho,u)
# TEVcompletionProb.m:91
        W=np.concat([0,- u(3),u(2),u(3),0,- u(1),- u(2),u(1),0])
# TEVcompletionProb.m:94
        RotMat_theta=eye(3) + dot(sin(theta(i)),W) + multiply((1 - cos(theta(i))),W ** 2)
# TEVcompletionProb.m:98
        u_new=dot(RotMat_theta,u_bent)
# TEVcompletionProb.m:99
        Rx[i + 1]=dot(L / N,u_new(1))
# TEVcompletionProb.m:102
        Ry[i + 1]=dot(L / N,u_new(2))
# TEVcompletionProb.m:103
        Rz[i + 1]=dot(L / N,u_new(3))
# TEVcompletionProb.m:104
    # Save 2 values for visualization
    if j == 1:
        RxVis[1,arange()]=Rx
# TEVcompletionProb.m:109
        RyVis[1,arange()]=Ry
# TEVcompletionProb.m:110
        RzVis[1,arange()]=Rz
# TEVcompletionProb.m:111
    if j == (aT / 2 + 1):
        RxVis[2,arange()]=Rx
# TEVcompletionProb.m:114
        RyVis[2,arange()]=Ry
# TEVcompletionProb.m:115
        RzVis[2,arange()]=Rz
# TEVcompletionProb.m:116
    Rx_sum[j]=sum(Rx)
# TEVcompletionProb.m:118
    Ry_sum[j]=sum(Ry)
# TEVcompletionProb.m:119
    Rz_sum[j]=sum(Rz)
# TEVcompletionProb.m:120
    Position[j,arange()]=np.concat([Rx_sum(j),Ry_sum(j),Rz_sum(j)])
# TEVcompletionProb.m:121
    R[j]=sqrt(Rx_sum(j) ** 2 + Ry_sum(j) ** 2 + Rz_sum(j) ** 2)
# TEVcompletionProb.m:123

# Plot the chain generated last, with linkers and DNA
# figure(2)
start1=np.concat([0,0,dCasLength])
# TEVcompletionProb.m:129

# Calculate start point for second linker
dCas9Vector=np.concat([[0],[0],[dCasLength]])
# TEVcompletionProb.m:132

n=DNAvector / np.norm(DNAvector)
# TEVcompletionProb.m:133

W=np.concat([0,- n(3),n(2),n(3),0,- n(1),- n(2),n(1),0])
# TEVcompletionProb.m:134
rotAngle=mod(dot(rotPerbp,bp),dot(2,pi))
# TEVcompletionProb.m:137
RotMat=eye(3) + dot(sin(rotAngle),W) + dot((1 - cos(rotAngle)),W ** 2)
# TEVcompletionProb.m:138
start2=DNAvector + (dot(RotMat,dCas9Vector)).T
# TEVcompletionProb.m:139

# cumsum is used for adding elements cumulatively for the coordinates of
# all segment chain joints
Rx_cumsum1=cumsum(RxVis(1,arange())) + start1(1)
# TEVcompletionProb.m:144
Ry_cumsum1=cumsum(RyVis(1,arange())) + start1(2)
# TEVcompletionProb.m:145
Rz_cumsum1=cumsum(RzVis(1,arange())) + start1(3)
# TEVcompletionProb.m:146
Rx_cumsum2=cumsum(RxVis(2,arange())) + start2(1)
# TEVcompletionProb.m:148
Ry_cumsum2=cumsum(RyVis(2,arange())) + start2(2)
# TEVcompletionProb.m:149
Rz_cumsum2=cumsum(RzVis(2,arange())) + start2(3)
# TEVcompletionProb.m:150
# Plot the total path
# hold('on')
np.plot3(np.concat([0,start1(1),Rx_cumsum1]),np.concat([0,start1(2),Ry_cumsum1]),np.concat([0,start1(3),Rz_cumsum1]),'k')
np.plot3(np.concat([DNAvector(1),start2(1),Rx_cumsum2]),np.concat([DNAvector(2),start2(2),Ry_cumsum2]),np.concat([DNAvector(3),start2(3),Rz_cumsum2]),'k')
np.plot3(np.concat([0,DNAvector(1)]),np.concat([0,DNAvector(2)]),np.concat([0,DNAvector(3)]),'k')
np.view(np.concat([377.6,36.63]))
np.axis('image')
np.xlabel('X')
np.ylabel('Y')
np.zlabel('Z')
# Calculate distances between endpoints
distance=zeros(1,aT / 2)
# TEVcompletionProb.m:167
for j in arange(1,(aT / 2)).reshape(-1):
    Position[j,arange()]=Position(j,arange()) + start1
# TEVcompletionProb.m:169
    Position[j + aT / 2,arange()]=Position(j + aT / 2,arange()) + start2
# TEVcompletionProb.m:170
    distance[j]=sqrt(sum((Position(j,arange()) - Position(j + aT / 2,arange())) ** 2))
# TEVcompletionProb.m:171

# Make histogram
#figure(3)
ti='Estimated probability density'
# TEVcompletionProb.m:176
lb=0
# TEVcompletionProb.m:177

ub=max(distance)
# TEVcompletionProb.m:178

step=(ub - lb) / (sqrt(aT / 2))
# TEVcompletionProb.m:179

y,edges=np.histcounts(distance,arange(lb,ub,step),nargout=2)
# TEVcompletionProb.m:180

x=arange(lb,(ub - step),step)
# TEVcompletionProb.m:181

k=1 / trapz(x,y)
# TEVcompletionProb.m:182

dist=dot(y,k)
# TEVcompletionProb.m:183
np.plot(x,dist,'k')
np.xlabel('Endpoint distance (nm)','interpreter','latex','FontSize',14)
np.ylabel('Probability density','interpreter','latex','FontSize',14)
#eval(np.sprintf('title('%s %s','interpreter','latex','FontSize',14)',ti,num2str(N)))
# Model contact probability: consider points in contact if the
# distance is less than 0.5 nm
ContactProb=0
# TEVcompletionProb.m:191
for i in arange(1,50).reshape(-1):
    if dot(step,i) < 0.5:
        ContactProb=ContactProb + dist(i)
# TEVcompletionProb.m:194

ContactProb=dot(ContactProb,step)
# TEVcompletionProb.m:197
disp('Probability of contact for linker length ' + L + ' nm:')
disp(ContactProb)
# Calculate number of iterations resulting in contact
num=dot(ContactProb,aT) / 2
# TEVcompletionProb.m:203
disp('Number of iterations resulting in contact: ' + num)