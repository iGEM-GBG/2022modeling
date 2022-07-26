% Main file for TEV completion
% Calculate TEV completion probability by generating many paths, 
% DNA is modeled as a straight line
clf, clc, clear, close all

% dCas9 diameter. 70 å = 7 nm
dCasLength = 7; % (nm) TODO: check actual length

N = 50;                 % Number of segments (less effect on smoothness than aT) %Doesn't seem to affect the result very much           
aT = 10^6 * 2;          % Number of iterations, 10^6 needed for calculating TEV contact probability
L = 11.5;                 % Protein linker chain length (1 amino acis: 3.4 - 4.0 Å usually for WLC)
bp = 30;                % number of base pairs (optimal: 10,20,30)
DL = 0.34*bp;           % DNA chain length (nm). 1 base pair: 0.34 nm (30bp: 10.2 nm)
rotPerbp = 0.626573;    % Rotation per base pair (35.9° = 0.626573 radians)
lp = 0.6;               % Protein persistance length

R = zeros(1,aT);            % Radius vector
Rx=zeros(1,N);              % x-direction
Ry=zeros(1,N);              % y-direction
Rz=zeros(1,N);              % z-direction
Rx_sum = zeros(1,aT);       % x-dim distances
Ry_sum = zeros(1,aT);       % y-dim distances
Rz_sum = zeros(1,aT);       % z-dim distances
Position = zeros(aT,3);     % endpoint positions

DNAvector = [DL,0,0];       % Vector representing DNA segment, length DL
    
sigma = sqrt(L/(N*lp)); % standard deviation

%Create probability distribution for theta (vikta efter omkrets för given vinkel i sfär):  
N_Dist=1000;
lb=0;
ub=pi;
Angle=linspace(lb,ub,N_Dist);
probdist_angles_core=sin(Angle)/2;
cdf_angles=cumsum(probdist_angles_core);
cdf_angles=cdf_angles/max(cdf_angles); %Normalization
cdf_angles(1)=0;
cdf_angles(end)=1;
ProbDist_theta = makedist('PiecewiseLinear','x',Angle,'Fx',cdf_angles);

RxVis = zeros(2,N);
RyVis = zeros(2,N);
RzVis = zeros(2,N);

RVis = zeros(2,N,3);

nvis = 0;
for j=1:aT    
    % Create initial segment with random angles
%     theta_0=random(ProbDist_theta,1); % 2*pi*(randn(1)-1/2);
%     phi_0=2*pi*(randn(1)-1/2);

    % Create initial segment with given angles
    if j <= aT/2    % Starting direction for first linker
        theta_0 = pi/2-0.01;
        phi_0 = 0+0.01;
    else            % Starting direction for second linker
        theta_0 = pi/2-0.01;
        phi_0 = pi-0.01; % pi
    end

    % (Rx, Ry, Rz) är en vektor med längd som pekar från ett segment till nästa segment. Rx representerar ej segmentens koordinat i x! 
    Rx(1) = L/N *cos(phi_0).*sin(theta_0);
    Ry(1) = L/N *sin(phi_0).*sin(theta_0);
    Rz(1) = L/N *cos(theta_0);
    
    rho = [normrnd(0,sigma,[1,N])];
    %rho = random(ProbDist_rho,[1,N]);
    theta=2*pi*(randn(1,N)-1/2);

    % ALGORITM:   
    for i=1:N-1
        % Copy previous vector u=(a,b,c)
        u=[Rx(i); Ry(i); Rz(i)] * N/L;
        
        % Rotation of vectors is done with Rodrigues' formula, using the matrix notation u_rot = R*u.
        % See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

        % Vector n, orthogonal to u, n = (1,-a/b,0)/norm((1,-a/b,0)) 
        n = [1,-u(1)/u(2),0]/norm([1,-u(1)/u(2),0]);
            
        % Rotation matrice for angle rho between new section and previous
        % section
        W=[ 0,    -n(3),   n(2)
            n(3),  0,     -n(1)
           -n(2),  n(1),   0   ];

        % RotMat_rho=eye(3) + sin(rho(i))*W + 2*sin(rho(i)/2).^2.*W^2; % version 1
        RotMat_rho=eye(3) + sin(rho(i))*W + (1-cos(rho(i)))*W^2;
        u_bent=RotMat_rho*u;

        % Rotera runt den böjda vektorn kring riktningen av den föregående vektorn v med theta  
        W=[ 0,    -u(3),   u(2)
            u(3),  0,     -u(1)
           -u(2),  u(1),   0  ];
        % RotMat_theta=eye(3)+sin(theta(i))*W+2*sin(theta(i)/2).^2.*W^2; % verion 1
        RotMat_theta=eye(3)+sin(theta(i))*W+(1-cos(theta(i))).*W^2;
        u_new=RotMat_theta*u_bent;

        %Store new vector
        Rx(i+1)=L/N * u_new(1);
        Ry(i+1)=L/N * u_new(2);
        Rz(i+1)=L/N * u_new(3);
    end    
    
    % Save 2 values for visualization
    if j == 1
        RxVis(1,:) = Rx;
        RyVis(1,:) = Ry;
        RzVis(1,:) = Rz;
    end
    if j == (aT/2+1)
        RxVis(2,:) = Rx;
        RyVis(2,:) = Ry;
        RzVis(2,:) = Rz;
    end
    Rx_sum(j) = sum(Rx); %https://en.wikipedia.org/wiki/Spherical_coordinate_system   "Spherical coordinates (r, θ, φ) as commonly used in physics"
    Ry_sum(j) = sum(Ry);
    Rz_sum(j) = sum(Rz);
    Position(j,:) = [Rx_sum(j), Ry_sum(j), Rz_sum(j)];
    
    R(j)=sqrt(Rx_sum(j)^2+Ry_sum(j)^2+Rz_sum(j)^2);
end

% Plot the chain generated last, with linkers and DNA
figure(2)

start1 = [0,0,dCasLength];  % startposition first segment

% Calculate start point for second linker
dCas9Vector = [0; 0; dCasLength]; % Vector representing dCas9
n = DNAvector/norm(DNAvector);    % Normed vector for rotation
W=[ 0,    -n(3),   n(2)           % Rotation matrix
    n(3),  0,     -n(1)
   -n(2),  n(1),   0   ];
rotAngle = mod(rotPerbp*bp, 2*pi);
RotMat=eye(3) + sin(rotAngle)*W + (1-cos(rotAngle))*W^2;
start2 = DNAvector + (RotMat*dCas9Vector).'; % startposition second segment


% cumsum is used for adding elements cumulatively for the coordinates of
% all segment chain joints
Rx_cumsum1 = cumsum(RxVis(1,:)) + start1(1);
Ry_cumsum1 = cumsum(RyVis(1,:)) + start1(2);
Rz_cumsum1 = cumsum(RzVis(1,:)) + start1(3);

Rx_cumsum2 = cumsum(RxVis(2,:)) + start2(1);
Ry_cumsum2 = cumsum(RyVis(2,:)) + start2(2);
Rz_cumsum2 = cumsum(RzVis(2,:)) + start2(3);

% Plot the total path
hold on;
plot3([0,start1(1),Rx_cumsum1], ...              % Plot the first linker
    [0,start1(2),Ry_cumsum1],[0,start1(3),Rz_cumsum1],'k') 
plot3([DNAvector(1),start2(1),Rx_cumsum2], ...   % Plot the second linker
    [DNAvector(2),start2(2),Ry_cumsum2],[DNAvector(3),start2(3),Rz_cumsum2],'k') 
plot3([0,DNAvector(1)],[0,DNAvector(2)],[0,DNAvector(3)],'k')
view([377.60 36.63])

axis image
xlabel("X")
ylabel("Y")
zlabel("Z")

% Calculate distances between endpoints
distance = zeros(1,aT/2);
for j=1:(aT/2)
    Position(j,:) = Position(j,:) + start1;             % Add start point coordinates
    Position(j+aT/2,:) = Position(j+aT/2,:) + start2;   % Add start point coordinates
    distance(j) = sqrt(sum((Position(j,:)-Position(j+aT/2,:)).^2));
end

% Make histogram
figure(3)
ti = "Estimated probability density";
lb = 0; % Lower boundary for the histogram. >0 to avoid noise
ub = max(distance); % Upper boundary for the histogram
step = (ub-lb)/(sqrt(aT/2)); % Step size for the points in the histogram
[y,edges] = histcounts(distance,lb:step:ub); % Creating points for histogram
x=lb:step:(ub-step); % Calculate the x-value for each data point in the histogram
k=1/trapz(x,y);  % Riemann sum
dist = y*k;
plot(x,dist,'k') % Plot the data points in the histogram
xlabel("Endpoint distance (nm)", 'interpreter','latex','FontSize',14)
ylabel("Probability density", 'interpreter','latex','FontSize',14)
eval(sprintf("title('%s %s','interpreter','latex','FontSize',14)",ti,num2str(N)))

% Model contact probability: consider points in contact if the
% distance is less than 0.5 nm
ContactProb = 0;
for i = 1:50
    if step*i < 0.5 % 5 å
        ContactProb = ContactProb + dist(i);
    end
end
ContactProb = ContactProb*step;

disp("Probability of contact for linker length " + L + " nm:") 
disp(ContactProb)

% Calculate number of iterations resulting in contact
num = ContactProb * aT/2;
disp("Number of iterations resulting in contact: " + num)