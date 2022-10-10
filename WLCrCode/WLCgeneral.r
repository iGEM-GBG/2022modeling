library('plotly')
library('pracma')

# Certain plots
#library('akima')
#library('rgl')
#library('mapplots')

# for 3D contour
library('misc3d')


Compound <- "prot"# "prot" or "DNA" to model protein or DNA respectively

N <- 50# Number of segments (less effect on smoothness than aT) #Doesn't seem to affect the result very much
aT <- 10^4# Number of chains    N=20  & 10^5 : 8.7 min  /  10^4: 
radiuses <- c(1,1) # Protein radiuses in nm

if (Compound == "prot"){
    L <- 18 # Protein chain length (nm)
    lp <- 0.45 # Protein linker persistence length (nm)
} else if (Compound == "DNA"){
    bp <- 30# Number of base pairs
    L <- 0.34*bp# DNA chain length (1 base pair: 0.34 nm, ~30 base pairs)
    lp <- 39# DNA persistence length (nm). Source: https://www.nature.com/articles/nphys2002
} else {
    L <- 10# Chain length (nm)
    lp <- 0.45# Persistence length (0.45 nm for some protein linkers, 40 nm for DNA)
}

sigma <- sqrt(L/(N*lp))# standard deviation

R <- integer(aT)# Radius vector
Rx <- zeros(aT,N)# x-direction
Ry <- zeros(aT,N)# y-direction
Rz <- zeros(aT,N)# z-direction

Position <- matrix(0,aT,3)# endpoint positions

start_time <- Sys.time()
for (j in 1:aT){

    # Create initial segment with direction in positive x direction
    theta_0 <- pi/2-0.001
    phi_0 <- 0+0.001
   
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

# Plot the chain generated first
# cumsum is used for adding elements cumulatively for the coordinates of
# all segment chain joints
Rx_cumsum <- cumsum(Rx[1,])
Ry_cumsum <- cumsum(Ry[1,])
Rz_cumsum <- cumsum(Rz[1,])

# Calculate protein midpoints
position1 <- c(-radiuses[1], 0, 0)
segmentCoord <- c(Rx[1,length(Rx[1,])], Ry[1,length(Ry[1,])], Rz[1,length(Rz[1,])])
position2 <- Position[1,] + radiuses[2]/Norm(segmentCoord) * segmentCoord

f <- function(x, y, z){
  x^2 + y^2 + z^2
}
get_sphere_contour <- function(R){
  x <- y <- z <- seq(-R, R, length.out = 100) 
  g <- expand.grid(x = x, y = y, z = z)
  voxel <- array(with(g, f(x, y, z)), dim = c(100, 100, 100))

  return (computeContour3d(voxel, level = R^2, x = x, y = y, z = z))
}

cont1 <- get_sphere_contour(radiuses[1])
cont2 <- get_sphere_contour(radiuses[2])

idx <- matrix(0:(nrow(cont1)-1), ncol=3, byrow=TRUE)

p <- plot_ly() %>%
  add_trace(x = c(0,Rx_cumsum),
            y = c(0,Ry_cumsum),
            z = c(0,Rz_cumsum),
            line = list(color = "pink", width = 5), 
            type = 'scatter3d', mode = 'lines',showlegend = FALSE)%>%
  layout(title = 'A Plotly Figure', paper_bgcolor='black',
         font = list(color = '#FFFFFF'),
         scene=list(xaxis=list(title="x-axis", gridcolor = '#FFFFFF'),
                    yaxis= list(title="y-axis", gridcolor = '#FFFFFF'),
                    zaxis= list(title="z-axis", gridcolor = '#FFFFFF')))%>%
  add_mesh(x = cont1[, 1] + position1[1], y = cont1[, 2] + position1[2], z = cont1[, 3] + position1[3],
           i = idx[, 1], j = idx[, 2], k = idx[, 3])%>%
  add_mesh(x = cont2[, 1] + position2[1], y = cont2[, 2] + position2[2], z = cont2[, 3] + position2[3],
           i = idx[, 1], j = idx[, 2], k = idx[, 3])

p

# Calculate vector from chain end to protein center
v = zeros(aT,3);
for (j in 1:aT) {
  normSegment <- c(Rx[j,length(Rx[j,])],Ry[j,length(Ry[j,])],Rz[j,length(Rz[j,])])
  normSegment = normSegment / Norm(normSegment)
  v[j,] = radiuses[2]*normSegment  # vectors added to center of TEV halves
}

# calculate positions of second protein and surface distances
prot_position <- zeros(aT,3) # center positions of second protein
distance <- integer(aT); # surface distances
for (j in (1:aT)){
  prot_position[j,] = Position[j,] + v[j,];
  distance[j] = Norm(position1 - prot_position[j,])-sum(radiuses)
}

# Create histogram
lb <- -sum(radiuses)# Lower boundary for the histogram.
ub <- max(distance)# Upper boundary for the histogram
step <- (ub-lb)/(sqrt(aT))# Step size for the points in the histogram
y <- histc(distance,seq(lb, ub, by = step))$cnt # Creating points for histogram
x <- seq(lb, ub, by = step)# Calculate the x-value for each data point in the histogram
k <- 1/trapz(x,y)# trapz: Trapezoidal numerical integration
dens <- y*k # Density estimation

histogram <- plot_ly()%>%
     add_trace(x = x,y = dens,
               line = list(color = "white"),
               type='scatter', mode = 'lines')%>%
     layout(title = paste("Surface distance probability density for polymer length ", L, "nm"), paper_bgcolor='black',
            plot_bgcolor = "black",
            font = list(color = '#FFFFFF'),
            xaxis=list(title=list(text ='Endpoint distance (nm)', font = list(color = "white")),gridcolor = 'white'),
            yaxis= list(title=list(text ='Probability density', font = list(color = "white")), gridcolor = 'white'))

histogram

# Calculate contact probability
ContactProb <- 0
for (i in 1:floor((ub-lb)/(step*2))){
  if (lb+step*i < 0.5){ # 5 å
    ContactProb = ContactProb + dens[i];
  }
}
ContactProb = ContactProb*step

