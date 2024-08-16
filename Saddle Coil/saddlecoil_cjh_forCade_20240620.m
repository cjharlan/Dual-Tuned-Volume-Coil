clear;
clc;

%Note units in meters
cr = 0.050;     % radius of cylinder surface on which helmholtz coil exists
cl = 0.10;     % length of the helmholtz coil

rFOV = 1;       % fraction of the cylinder radius to calculate fields over
rnocalc=0.001;  % minimum distance away from conductor for calculations

Ro=zeros(3,128,128);

%Axial view:
% Ro(1,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128),128,1);
% Ro(2,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128)',1,128);

%Sagittal view:
Ro(2,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128)',1,128);
Ro(3,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128),128,1);

%Coronal view:
% Ro(1,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128)',1,128);
% Ro(3,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128),128,1);

b1a=bfieldca(pi/2,pi/2,cr,cl,Ro,rnocalc);   % top element
b1b=bfieldca(3*pi/2,pi/2,cr,cl,Ro,rnocalc); % bottom element
bt=squeeze((b1a(1,:,:)+1i*b1a(2,:,:)) - (b1b(1,:,:)+1i*b1b(2,:,:)));

%find value at center of bt
centervalue=bt(round(size(bt,1)/2),round(size(bt,2)/2));

%normalize to value at center 
bt = (bt/centervalue)*100;

%Note below: real(bt) are B1x, imag(bt) are B1y

%contour levels for +-10% deviation from isocenter of the coil
level = [90,110];
level2 = [95,105];
level3 = [97,103];
%level4 = [99,101];

%image
imagesc(squeeze(Ro(3,1,:)),squeeze(Ro(2,:,1)),abs(squeeze(bt)));
colorbar
hold on

%contour plot
contour(squeeze(Ro(3,1,:)),squeeze(Ro(2,:,1)),abs(squeeze(bt)),level,'LineColor','r')
hold on
contour(squeeze(Ro(3,1,:)),squeeze(Ro(2,:,1)),abs(squeeze(bt)),level2,'--','LineColor','r')
hold on
contour(squeeze(Ro(3,1,:)),squeeze(Ro(2,:,1)),abs(squeeze(bt)),level3,':','LineColor','r')
hold on
%contour(squeeze(Ro(3,1,:)),squeeze(Ro(2,:,1)),abs(squeeze(bt)),level4,'-.','LineColor','b')
grid on
set(gca, 'fontsize', 12, 'fontweight', 'bold')
hold on
plot(0,0,'x','markersize',15,'color','r')
colormap(flipud(brewermap([],'Spectral')))
zMin = 20;
zMax = 180;
caxis([zMin, zMax]);
title(['Sagittal View, r = ',num2str(cr),' m, L = ',num2str(cl),' m']);

xlim([-0.01 0.01])
ylim([-0.01 0.01])