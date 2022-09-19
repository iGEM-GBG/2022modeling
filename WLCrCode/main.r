library('plotly')
library('pracma')

rm(list = ls())

# dCas9 diameter. 70 å = 7 nm
dCasLength <- 7# (nm) TODO: check actual length

N <- 20# Number of segments (less effect on smoothness than aT) #Doesn't seem to affect the result very much
aT <- 10^2 * 2# Number of iterations, 10^6 needed for calculating TEV contact probability
L <- 18# Protein linker chain length (1 amino acis: 3.4 - 4.0 Å usually for WLC)
bp <- 30# number of base pairs (optimal: 10,20,30)
DL <- 0.34*bp# DNA chain length (nm). 1 base pair: 0.34 nm (30bp: 10.2 nm)
rotPerbp <- 2*pi/10.5# Rotation per base pair (10.5 bp per turn - https://www.pnas.org/doi/10.1073/pnas.75.2.640)
lp <- 0.6# Protein persistence length

R <- integer(aT)# Radius vector
Rx <- matrix(0,aT,N)# x-direction
Ry <- matrix(0,aT,N)# y-direction
Rz <- matrix(0,aT,N)# z-direction

Position <- matrix(0,aT,3)# endpoint positions

DNAvector <- c(DL,0,0)# Vector representing DNA segment, length DL

sigma <- sqrt(L/(N*lp))# standard deviation

RxVis <- zeros(2,N)
RyVis <- zeros(2,N)
RzVis <- zeros(2,N)

RVis <- array(rep(0,2*N*3),dim=c(2,N,3))

nvis <- 0
start_time <- Sys.time()
for (j in 1:aT){
    # Create initial segment with random angles
    # theta_0= 2*pi*(rand(1)-1/2);
    # phi_0=2*pi*(rand(1)-1/2);

    # Create initial segment with given angles
    if (j <= aT/2){# Starting direction for first linker
        theta_0 <- pi/2-0.01
        phi_0 <- 0+0.01
    } else {# Starting direction for second linker
        theta_0 <- pi/2-0.01
        phi_0 <- pi-0.01# pi
    }

    # (Rx[i], Ry[i], Rz[i]) is a vector representing a segments direction and length.
    Rx[j,1] <- L/N *cos(phi_0)*sin(theta_0)
    Ry[j,1] <- L/N *sin(phi_0)*sin(theta_0)
    Rz[j,1] <- L/N *cos(theta_0)

    rho <- rnorm(N,0,sigma)
    #rho = random(ProbDist_rho,[1,N]);
    theta <- 2*pi*(rand(1,N)-1/2)

    # ALGORITM:
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

    # Save 2 values for visualization
    if (j == 1){
        RxVis[1,] <- Rx[j,]
        RyVis[1,] <- Ry[j,]
        RzVis[1,] <- Rz[j,]
    }
    start1 <- c(0,0,7)
    start2 <- c(10.2000000, -0.3567985, 7.0000000)
    pos1 = Position[1,]+start1
    pos2 = Position[j,]+start2
    if (Norm(pos1-pos2) < 1){
        RxVis[2,] <- Rx[j,]
        RyVis[2,] <- Ry[j,]
        RzVis[2,] <- Rz[j,]
    }
}
end_time <- Sys.time()
forloop_time <- end_time - start_time

# Add protein as balls on endpoints
r1 = 1        #radius (nm) for protein at first linker's endpoint 
r2 = 2        #radius (nm) for protein at second linker's endpoint
v = zeros(aT,3);    # Vector to middle of proteins, with same direction as last segment
for (j in 1:aT/2){
  v[j,] = r1*c(tail(Rx[j,],1),tail(Rx[j,],1),tail(Rx[j,],1))/Norm(c(tail(Rx[j,],1),tail(Rx[j,],1),tail(Rx[j,],1)))# vectors added to center of TEV halves
}
for (j in (aT/2+1):aT){
  v[j,] = r2*c(tail(Rx[j,],1),tail(Rx[j,],1),tail(Rx[j,],1))/Norm(c(tail(Rx[j,],1),tail(Rx[j,],1),tail(Rx[j,],1)))# vectors added to center of TEV halves
}

# Plot an example of generated linkers with endpoints in contact
start1 <- t(t(c(0,0,dCasLength)))# start position for first segment

# Calculate start point for second linker
dCas9Vector <- matrix(c(0, 0, dCasLength),nrow = 3)# Vector representing dCas9
n <- DNAvector/Norm(DNAvector)# Normalized vector for rotation
W <- array(c(0, -n[3], n[2],
          n[3], 0, -n[1],
         -n[2], n[1], 0), dim=c(3,3))
rotAngle <- rotPerbp*bp
RotMat <- eye(3) + sin(rotAngle)*W + (1-cos(rotAngle))*W^2
start2 <- DNAvector + RotMat %*% dCas9Vector # startposition second segment'



# cumsum is used for adding elements cumulatively for the coordinates of
# all segment chain joints
Rx_cumsum1 <- cumsum(RxVis[1,]) + start1[1]
Ry_cumsum1 <- cumsum(RyVis[1,]) + start1[2]
Rz_cumsum1 <- cumsum(RzVis[1,]) + start1[3]

Rx_cumsum2 <- cumsum(RxVis[2,]) + start2[1]
Ry_cumsum2 <- cumsum(RyVis[2,]) + start2[2]
Rz_cumsum2 <- cumsum(RzVis[2,]) + start2[3]


# df1 <- data.frame(c(start1[1],Rx_cumsum1), 
#                   c(start1[2],Ry_cumsum1), 
#                   c(start1[3],Rz_cumsum1))
# 
# df2 <- data.frame(c(start2[1],Rx_cumsum2), 
#                   c(start2[2],Ry_cumsum2),
#                   c(start2[3],Rz_cumsum2))

# Plot the total path
plot_ly() %>%
  add_trace(x = ~c(start1[1],Rx_cumsum1),
            y = ~c(start1[2],Ry_cumsum1),
            z = ~c(start1[3],Rz_cumsum1),
            line = list(color = "pink", width = 5), 
            type = 'scatter3d', mode = 'lines', name='Linker 1')%>%
  add_trace(x = ~c(start2[1],Rx_cumsum2), 
            y = ~c(start2[2],Ry_cumsum2), 
            z = ~c(start2[3],Rz_cumsum2),
            line = list(color = "lightblue", width = 5), 
            type = 'scatter3d', mode = 'lines', name='Linker 2')%>%
  add_markers(x = ~start1[1],
              y = ~start1[2],
              z = ~start1[3],
              name = 'Start point', marker = list(color = "blue", size=25))%>%
  add_markers(x = ~start2[1],
              y = ~start2[2],
              z = ~start2[3],
              name = 'Start point', marker = list(color = "blue", size=25))%>%
  add_markers(x = ~tail(Rx_cumsum1,1),
              y = ~tail(Ry_cumsum1,1),
              z = ~tail(Rz_cumsum1,1),
              name = 'End point linker 2', marker = list(color = "red", size=15))%>%
  add_markers(x = ~tail(Rx_cumsum2,1),
              y = ~tail(Ry_cumsum2,1),
              z = ~tail(Rz_cumsum2,1),
              name = 'End point linker 2', marker = list(color = "red", size=15))%>%
  layout(title = 'A Plotly Figure', paper_bgcolor='black',
         font = list(color = '#FFFFFF'),
         scene=list(xaxis=list(title="x-axis", gridcolor = '#FFFFFF'),
                    yaxis= list(title="y-axis", gridcolor = '#FFFFFF'),
                    zaxis= list(title="z-axis", gridcolor = '#FFFFFF')))

# plot_ly(df12) %>%
#   add_markers(x = ~c(0,start1[1],Rx_cumsum1), y = ~c(0,start1[2],Ry_cumsum1), z = ~c(0,start1[3],Rz_cumsum1), 
#               type = 'scatter3d', mode = 'lines')  %>% # Plot the first linker
#   add_markers(x = ~c(0,start1[1],Rx_cumsum1), y = ~c(0,start1[2],Ry_cumsum1), z = ~c(0,start1[3],Rz_cumsum1), 
#               type = 'scatter3d', mode = 'lines')  %>% # Plot the first linker

# Calculate distances between protein surfaces
distance <- integer(aT/2)
for (j in 1:(aT/2)){
    Position[j,] <- Position[j,] + start1 + v[j,]# Add start point coordinates
    Position[j+aT/2,] <- Position[j+aT/2,] + start2 + v[j+aT/2,]# Add start point coordinates
    distance[j] <- Norm(Position[j,]-Position[j+aT/2,]) - (r1 + r2)# Distance between protein midpoints - their radiuses
}

# Make histogram
#figure(3)
ti <- "Estimated probability density"
lb <- -(r1+r2)# Lower boundary for the histogram. >0 to avoid noise
ub <- max(distance)# Upper boundary for the histogram
step <- (ub-lb)/(sqrt(aT/2))# Step size for the points in the histogram
y <- histc(distance,seq(lb, ub, by = step))$cnt# Creating points for histogram
x <- seq(lb, ub, by = step)# Calculate the x-value for each data point in the histogram
k <- 1/trapz(x,y)# Riemann sum
dens <- y*k  # Density estimation
plot(x,dens,type='l', main="Histogram for endpoint distances", # Plot the data points in the histogram
     xlab="Endpoint distance (nm)", ylab="Approximated probability density")

# Model contact probability: consider points in contact if the
# distance is less than 0.5 nm
ContactProb <- 0
for (i in 1:50){
    if (lb+step*i < 0.5){# a distance of 3-5 å between proteins is considered as contact
        ContactProb <- ContactProb + dens[i]
    }
}
ContactProb <- ContactProb*step

disp("Probability of contact for linker length ", L, " nm:")
disp(ContactProb)

# Calculate number of iterations resulting in contact
num <- ContactProb * aT/2
disp("Number of iterations resulting in contact: ", floor(num))

