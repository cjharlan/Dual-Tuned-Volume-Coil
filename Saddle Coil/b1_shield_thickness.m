clear;
clc;
close all;

%blengths are all 2.5x the value of the bdiameter
blen1 = 0.1375; %mm
bdiam1 = 0.055; %mm

blen2 = 0.1125; %mm
bdiam2 = 0.045; %mm

blen3 = 0.100; %mm
bdiam3 = 0.040; %mm

blen4 = 0.0875; %mm
bdiam4 = 0.035; %mm

blen5 = 0.095; %mm
bdiam5 = 0.038; %mm

nleg = 12; %number of rungs 
npts = 64;
fovz = 100;
fovx = 100;

one = birdfield_cjh(blen1, bdiam1, nleg, npts, fovz, fovx);
two = birdfield_cjh(blen2, bdiam2, nleg, npts, 110, 110);
three = birdfield_cjh(blen3, bdiam3, nleg, npts, fovz, fovx);
four = birdfield_cjh(blen4, bdiam4, nleg, npts, fovz, fovx);
five = birdfield_cjh(blen5, bdiam5, nleg, npts, fovz, fovx);
