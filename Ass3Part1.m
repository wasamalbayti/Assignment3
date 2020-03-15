% set up constant
q=-1.602*10^-19;
k=1.38064*10^-23;
w=200*10^-9;
l=100*10^-9;
m=0.26*9.109e-31;
%electric field
Ex= -(0.1-0)/(w)
%Force on each electron
Fx=Ex*q
%Accleration on electron
accelerationX=Fx/m
%Nominal size of region is 200nmx100nm
vt=@(t) sqrt(k*t/m);
thermalV= vt(300);
tmn=0.2*10^-12;
%mean free path
freepath=tmn*vt(300);
%set time step
deltaT=7.5*10^-16;
numParticles=1000;
partplot=randi([1,1000],1,20);

%initialize each particle's coordinates
verticalArray=rand(numParticles,1)*w;
horrizArray=rand(numParticles,1)*l;
% initizile velocity also
xVelocity=randn(numParticles,1).*thermalV/sqrt(2);
yVelocity=randn(numParticles,1).*thermalV/sqrt(2);
VoltageRms=sqrt(xVelocity.^2+yVelocity.^2);

currentx=zeros(1,1000);
averagexvel=zeros(1,1000);

velocityMean=mean(VoltageRms);

for i=1:1000
    
    %3. new Xvelocity with acceleration 
    xVelocity=xVelocity+accelerationX*deltaT;
    
    %boundary conditions in vertical array
    IT=(verticalArray>=w);
    yVelocity(IT)=-yVelocity(IT);
    verticalArray(IT)=(verticalArray(IT)-2*(verticalArray(IT)-w));
    
    IT=(verticalArray<=0);
    yVelocity(IT)=-yVelocity(IT);
    verticalArray(IT)=(verticalArray(IT)+2*(0-verticalArray(IT)));
    
    %Boundary conditions in horriziontal array
    horrizArray(horrizArray>=l)=horrizArray(horrizArray>=l) - l;
    horrizArray(horrizArray<=0)=horrizArray(horrizArray<=0)+l;
    
    %Update position with velocity
    horrizArray=horrizArray+xVelocity.*deltaT;
    verticalArray=verticalArray+yVelocity.*deltaT;
    
    %scattering of the electrons:
    pscat=1-exp(-deltaT/(0.2*10^-12));
    a=rand(numParticles,1);
    
    xVelocity(a < pscat)=randn(sum(a < pscat),1).*thermalV/sqrt(2);
    yVelocity(a < pscat)=randn(sum(a < pscat),1).*thermalV/sqrt(2);
    VoltageRms=sqrt(xVelocity.^2+yVelocity.^2);
    
    %Temperature vs time
    subplot(2,1,2)
    currentx(i)=mean(xVelocity)*-q*w*1e15*(100^2);
    plot(linspace(1,1000,1000),currentx);
    title(['Current Density'])
    subplot(2,1,1)
    plot(horrizArray(partplot),verticalArray(partplot),'.')
    title(['Electron movement']);
    xlim([0 l])
    ylim([0 w])
    hold on
    pause(.01)
end

[Ex,Ey]=meshgrid(0:l/20:l,0:w/20:w);
countPart=zeros(20,20);
tempcheck=zeros(20,20);
counter=0;
totalV=0;
for i=1:20
    txmn=Ex(1,i);
    txmx=Ex(1,i+1);
    for y =1:20
        tymn=Ey(y,1);
        tymx=Ey(y+1,1);
        for j=1:numParticles
            if(horrizArray(j)>txmn && horrizArray(j)<txmx && verticalArray(j)<tymx && verticalArray(j)>tymn)
                counter=counter+1;
                countPart(y,i)=countPart(y,i)+1;
                totalV=totalV+sqrt(xVelocity(j)^2+yVelocity(j)^2);
                if(counter~=0)
                    tempcheck(y,i)=m*(totalV^2)/(counter*k);
                end
            end
        end
        totalV=0;
        counter=0;
    end
end
%electron density map
figure(2)
surf(flipud(countPart))
title('Electron Density map')
zlabel('number of particles')

%temperature map
figure(3)
surf(flipud(tempcheck))
title('Temperature Density Map')
zlabel('Temperature')

