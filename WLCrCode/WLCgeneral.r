library('plotly')
library('pracma')

Compound <- "DNA"# "prot" or "DNA" to model protein or DNA respectively

N <- 50# Number of segments (less effect on smoothness than aT) #Doesn't seem to affect the result very much
aT <- 10^5# Number of chains    N=20  & 10^5 : 17 sek  /  10^4: 27.8 sek, 10^5: 7.9 min

if (Compound == "prot"){
    L <- 10# Protein chain length (nm)
    lp <- 0.6# Protein persistence length (nm)
} else if (Compound == "DNA"){
    bp <- 30# Number of base pairs
    L <- 0.34*bp# DNA chain length (1 base pair: 0.34 nm, ~30 base pairs)
    lp <- 39# DNA persistence length (nm). Source: https://www.nature.com/articles/nphys2002
} else {
    L <- 10# Chain length (nm)
    lp <- 0.45# Persistence length (0.45 nm for some protein linkers, 40 nm for DNA)
}

sigma <- sqrt(L/(N*lp))# standard deviation

R <- zeros(1,aT)# Radius vector
Rx <- zeros(aT,N)# x-direction
Ry <- zeros(aT,N)# y-direction
Rz <- zeros(aT,N)# z-direction

Position <- matrix(0,aT,3)# endpoint positions

start_time <- Sys.time()
for (j in 1:aT){

    # Create initial segment with given angles
    if (j <= aT/2){# Starting direction for first linker
        theta_0 <- pi/2-0.001
        phi_0 <- 0+0.001
    } else {# Starting direction for second linker
        theta_0 <- pi/2-0.001
        phi_0 <- pi-0.001# pi
    }

    # (Rx[i], Ry[i], Rz[i]) is a vector representing a segments direction and length.
    Rx[j,1] <- L/N *cos(phi_0)*sin(theta_0)
    Ry[j,1] <- L/N *sin(phi_0)*sin(theta_0)
    Rz[j,1] <- L/N *cos(theta_0)
    

    rho <- rnorm(N,0,sigma) #Angle between consecutive segments in linker. Diffrent angles have diffrent values
    theta <- 2*pi*(rand(1,N)-1/2) #Rotation between consecutive segments in lnker. All angles have equal probability

    # Making a linker from generated angles
    for (i in 1:(N-1)){
      # Copy previous vector u=(a,b,c)
      
      u <- c(Rx[j,i], Ry[j,i], Rz[j,i])
      
      # Rotation of vectors is done with Rodrigues' formula, using the matrix notation u_rot = R*u.
      # See https://en.wikipedia.org/wiki/Rodrigues#27_rotation_formula
      
      # Vector n, orthogonal to u, n = (1,-a/b,0)/norm((1,-a/b,0))
      #disp(c(1,-u[1]/u[2],0))
      n <- c(1,-u[1]/u[2],0)/Norm(c(1,-u[1]/u[2],0))
      
      # Rotation matrice for angle rho between new section and previous
      # section
      W <- t(array(c(0, -n[3], n[2], # elements added column by column, so transpose needed
                     n[3],  0,     -n[1],
                     -n[2],  n[1],   0   ), dim=c(3,3)))
      
      # RotMat_rho=eye(3) + sin(rho(i))*W + 2*sin(rho(i)/2).^2.*W^2; # version 1
      RotMat_rho <- eye(3) + sin(rho[i])*W + (1-cos(rho[i]))* (W %*% W)
      u_bent <- RotMat_rho %*% u
      
      # Rotate the bent vector around the last vector by uniformly distributed angle theta
      u = u/Norm(u)
      W <- t(array(c( 0,    -u[3],   u[2],
                      u[3],  0,     -u[1],
                      -u[2],  u[1],   0  ), dim=c(3,3)))
      # RotMat_theta=eye(3)+sin(theta(i))*W+2*sin(theta(i)/2).^2.*W^2; # verion 1
      RotMat_theta <- eye(3)+sin(theta[i])*W+(1-cos(theta[i]))*W %*% W
      u_new <- RotMat_theta %*% u_bent
      
      #Store new vector
      Rx[j,i+1] <- u_new[1]
      Ry[j,i+1] <- u_new[2]
      Rz[j,i+1] <- u_new[3]
    }
    Position[j,] <- c(sum(Rx[j,]), sum(Ry[j,]), sum(Rz[j,]))
    
    R[j] <- Norm(Position)
}
end_time <- Sys.time()
disp(end_time - start_time)

# Plot the chain generated last
figure(2)
# cumsum is used for adding elements cumulatively for the coordinates of
# all segment chain joints
Rx_cumsum <- cumsum(Rx(},:))
Ry_cumsum <- cumsum(Ry(},:))
Rz_cumsum <- cumsum(Rz(},:))

plot3(Rx_cumsum,Ry_cumsum,Rz_cumsum,'k')#Plot the chain!
axis image
xlabel("X")
ylabel("Y")
zlabel("Z")


# Create histogram
ti <- "Estimated probability density"
figure(3)
lb <- 0# Lower boundary for the histogram. >0 to avoid noise
ub <- L# Upper boundary for the histogram
step <- (ub-lb)/(2*sqrt(aT))# Step size for the points in the histogram
[y,edges] <- histcounts(R,lb:step:ub)# Creating points for histogram
x <- lb:step:(ub-step)# Calculate the x-value for each data point in the histogram
k <- 1/trapz(x,y)# trapz: Trapezoidal numerical integration
plot(x,y*k,'b')# Plot the data points in the histogram
xlabel("Chain reach [$nm$]", 'interpreter','latex','FontSize',14)
ylabel("Density function", 'interpreter','latex','FontSize',14)
title(ti)
#eval(sprintf("title('#s #s','interpreter','latex','FontSize',14)",ti,num2str(N)))
