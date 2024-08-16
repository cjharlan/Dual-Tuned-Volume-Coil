
%Note units in meters
cr = 0.010;     % radius of cylinder surface on which helmholtz coil exists
cl = 0.020;     % length of the helmholtz coil
rFOV = 1;       % fraction of the cylinder radius to calculate fields over
rnocalc=0.001;  % minimum distance away from conductor for calculations

Ro=zeros(3,128,128);

%Axial view:
Ro(1,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128),128,1);
Ro(2,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128)',1,128);

%Sagittal view:
%Ro(2,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128)',1,128);
%Ro(3,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128),128,1);

%Coronal view:
%Ro(1,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128)',1,128);
%Ro(3,:,:)=repmat(linspace(-cr*rFOV,cr*rFOV,128),128,1);

b1a=b1coil(0.1,0.1,2.094,Ro);
%b1a=bfieldca(pi/2,pi/2,cr,cl,Ro,rnocalc);   % top element
%b1b=bfieldca(3*pi/2,pi/2,cr,cl,Ro,rnocalc); % bottom element
bt=squeeze((b1a(1,:,:)+1i*b1a(2,:,:))); %- (b1b(1,:,:)+1i*b1b(2,:,:)))

%Note below: real(bt) are B1x, imag(bt) are B1y
imagesc(squeeze(Ro(3,1,:)),squeeze(Ro(2,:,1)),abs(squeeze(bt)))
