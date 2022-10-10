
# Model polymer endpoint distances and endpoint connection
onePolymer <- function(range, N, length_count, radiuses, aT, lp) {
  onePolymerSupportResults <- list()
  
  lengths <- linspace(range[1], range[2], length_count)
  for (i in 1:length_count) {
    onePolymerSupportResults = append(onePolymerSupportResults, 
                                      list(onePolymerSupport(N = N, aT = aT, L = lengths[i], lp = lp, radiuses = radiuses)))
  }
  # From onePolymerSupport: example_chain, distance_histogram, 
  
  # aT = 10^3: 1.3 s, 10^4: 28 s, 
  # resultList <- list("text" = "hello", "example_chains" = example_chains
  
  return (onePolymerSupportResults)
}

# Model distance and connection for two polymers
twoPolymers <- function(some_text) {
  text <- some_text
  number <- 8
  resultList <- list("text" = text, "number" = number)
  return (resultList)
}

# Single polymer support function
onePolymerSupport <- function(N, aT, L, lp, radiuses) {
  # N:  Number of segments
  # at: Number of chains
  # lp: Persistence length (nm)
  # L: Chain contour length
  
  sigma <- sqrt(L/(N*lp))# standard deviation
  
  R <- integer(aT)# Radius vector
  Rx <- zeros(aT,N)# x-direction
  Ry <- zeros(aT,N)# y-direction
  Rz <- zeros(aT,N)# z-direction
  
  Position <- matrix(0,aT,3)# endpoint positions
  
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
  
  cont1 <- get_sphere_contour(radiuses[1])
  cont2 <- get_sphere_contour(radiuses[2])
  
  #print(cont1)
  idx <- matrix(0:(nrow(cont1)-1), ncol=3, byrow=TRUE)
  
  p <- plot_ly() %>%
    add_trace(x = c(0,Rx_cumsum),
              y = c(0,Ry_cumsum),
              z = c(0,Rz_cumsum),
              line = list(color = "pink", width = 5), 
              type = 'scatter3d', mode = 'lines', showlegend = FALSE)%>%
    layout(title = paste('Example simulation with polymer length',L, "nm"), paper_bgcolor='black',
           font = list(color = '#FFFFFF'),
           scene=list(xaxis=list(title="x-axis", gridcolor = '#FFFFFF'),
                      yaxis= list(title="y-axis", gridcolor = '#FFFFFF'),
                      zaxis= list(title="z-axis", gridcolor = '#FFFFFF')))%>%
    add_mesh(x = cont1[, 1] + position1[1], y = cont1[, 2] + position1[2], z = cont1[, 3] + position1[3],
             i = idx[, 1], j = idx[, 2], k = idx[, 3])%>%
    add_mesh(x = cont2[, 1] + position2[1], y = cont2[, 2] + position2[2], z = cont2[, 3] + position2[3],
             i = idx[, 1], j = idx[, 2], k = idx[, 3])
  
  
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
  histogram <- plot_ly() %>%
    add_trace(x = x,y = dens,
              line = list(color = "white"),
              type='scatter', mode = 'lines')%>%
    layout(title = paste("Surface distance probability density for polymer length ", L, "nm"), paper_bgcolor='black',
           plot_bgcolor = "black",
           font = list(color = '#FFFFFF'),
           xaxis=list(title=list(text ='Endpoint distance (nm)', font = list(color = "white")),gridcolor = 'white'),
           yaxis= list(title=list(text ='Probability density', font = list(color = "white")), gridcolor = 'white'))
  

  # Calculate contact probability
  ContactProb <- 0
  for (i in 1:floor((ub-lb)/(step*2))){
    if (lb+step*i < 0.5){ # 5 å
      ContactProb = ContactProb + dens[i];
    }
  }
  ContactProb = ContactProb*step
  
  
  resultList <- list("example_chain" = p, "distance_histogram" = histogram, "contact_prob" = ContactProb)
  
  return (resultList)
}

# refactored code
get_sphere_contour <- function(R){
  if (R == 0) { R = 0.01 } # There will be errors if R = 0
  x <- y <- z <- seq(-R, R, length.out = 100) 
  g <- expand.grid(x = x, y = y, z = z)
  voxel <- array(with(g, f(x, y, z)), dim = c(100, 100, 100))
  
  #print(R)
  #print(voxel)
  cont <- computeContour3d(voxel, level = R^2, x = x, y = y, z = z)

  return (cont)
}
f <- function(x, y, z){
  return (x^2 + y^2 + z^2)
}

twoPolymersSupport <- function() {
  
}