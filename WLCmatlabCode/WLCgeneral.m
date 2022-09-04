clc, clear, close all

Compound = "prot"; % "prot" or "DNA" to model protein or DNA respectively
disp("The chosen compound is " + Compound)

N = 50;                 % Number of segments (less effect on smoothness than aT) %Doesn't seem to affect the result very much           
aT = 10^4;              % Number of chains    N=20  & 10^5 : 17 sek  /  10^6:169 sek

if Compound == "prot"
    L = 10;               % Protein chain length (nm)
    lp=0.6;                 % Persistence length (nm)
elseif Compound == "DNA"
    L = 10.2;               % DNA chain length (10.2 nm) (1 base pair: 0.34 nm, ~30 base pairs)
    lp=40;                  % Persistence length (nm)
else
    L = 1.75;               % Chain length (nm)
    lp=0.01;                % Persistence length (0.6 nm for protein, 40 nm for DNA)
end

R = zeros(1,aT);        % Radius vector
Rx=zeros(1,N);          % x-direction
Ry=zeros(1,N);          % y-direction
Rz=zeros(1,N);          % z-direction
Rx_sum = zeros(1,aT);       % x-dim distances
Ry_sum = zeros(1,aT);       % y-dim distances
Rz_sum = zeros(1,aT);       % z-dim distances
    

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
%     rw=random(ProbDist_angles,[1000000,1]);
%     histogram(rw)

%Skapa sannolikhetsfördelning för rho (Normalfördelning, men inom -pi:pi):  
N_Dist=1000;
lb=-pi;
ub=pi;
Angle=linspace(lb,ub,N_Dist);
probdist_angles_core=exp(-Angle.^2/(2*sigma^2)); %Ingen faktor framför exp() därför det löses automatiskt genom skalning i slutet då arean ska bli =1. Endast form är relevant. 
cdf_angles=cumsum(probdist_angles_core);
cdf_angles=cdf_angles/max(cdf_angles); %Normalization
cdf_angles(1)=0;
cdf_angles(end)=1;
ProbDist_rho = makedist('PiecewiseLinear','x',Angle,'Fx',cdf_angles);
%     rw=random(ProbDist_angles,[10000000,1]);
%     histogram(rw)


for j=1:aT    
    % Create initial segment with random angles
%     theta_0=random(ProbDist_theta,1); % 2*pi*(randn(1)-1/2);
%     phi_0 = 2*pi*(randn(1)-1/2);
    
    % Create initial segment with given angle
    theta_0 = 1;
    phi_0 = 1;

    % (Rx, Ry, Rz) är en vektor med längd som pekar från ett segment till nästa segment. Rx representerar ej segmentens koordinat i x! 
    Rx(1) = L/N *cos(phi_0).*sin(theta_0); %
    Ry(1) = L/N *sin(phi_0).*sin(theta_0); %
    Rz(1) = L/N *cos(theta_0);             %
    
    rho = [normrnd(0,sigma,[1,N])];
    %rho = random(ProbDist_rho,[1,N]);
    theta=2*pi*(randn(1,N)-1/2);

    % ALGORITM:   
    for i=1:N-1
        % Copy previous vector and normalize
        u=[Rx(i); Ry(i); Rz(i)] * N/L;
        
        % Rotation of vectors is done with Rodrigues' formula, using the matrix notation u_rot = R*u.
        % See https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

        % Vector n, orthogonal to u, n = (1,-a/b,0)/norm((1,-a/b,0)) 
        n = [1,-u(1)/u(2),0]/norm([1,-u(1)/u(2),0]);
            
        % Rotation matrice for angle rho between new section and previous section
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


    % Plot endpoints of all generated chains
    % Let aT be 10^2 for short running time
%     figure(1)
%     scatter3(Rx,Ry,Rz,'k.')  % Plot endpoints   
%     axis equal; hold on
%     scatter3(0,0,0,'r') % Plot origo as a red ring
        

    Rx_sum(j) = sum(Rx); %https://en.wikipedia.org/wiki/Spherical_coordinate_system   "Spherical coordinates (r, θ, φ) as commonly used in physics"
    Ry_sum(j) = sum(Ry);
    Rz_sum(j) = sum(Rz);
 
    R(j)=sqrt(Rx_sum(j)^2+Ry_sum(j)^2+Rz_sum(j)^2);
end

% Plot the chain generated last
figure(2)
% cumsum is used for adding elements cumulatively for the coordinates of
% all segment chain joints
Rx_cumsum = cumsum(Rx);     %https://en.wikipedia.org/wiki/Spherical_coordinate_system   "Spherical coordinates (r, θ, φ) as commonly used in physics"
Ry_cumsum = cumsum(Ry);
Rz_cumsum = cumsum(Rz);

plot3(Rx_cumsum,Ry_cumsum,Rz_cumsum,'k') %Plot the chain!
axis image
xlabel("X")
ylabel("Y")
zlabel("Z")


% Create histogram
% ti=char(input("Title? (Title + N) "));
ti = "Estimated probability density";
figure(3)
lb = 0; % Lower boundary for the histogram. >0 to avoid noise
ub = L; % Upper boundary for the histogram
step = (ub-lb)*5/(sqrt(aT)); % Step size for the points in the histogram
[y,edges] = histcounts(R,lb:step:ub); % Creating points for histogram
x=lb:step:(ub-step); % Calculate the x-value for each data point in the histogram
k=1/trapz(x,y);  % Riemann sum
plot(x,y*k,'b.') % Plot the data points in the histogram
xlabel("Chain reach [$nm$]", 'interpreter','latex','FontSize',14)
ylabel("Probability density [1/$nm$]", 'interpreter','latex','FontSize',14)
eval(sprintf("title('%s %s','interpreter','latex','FontSize',14)",ti,num2str(N)))
