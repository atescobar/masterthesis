k = 400; % thermal conductivity of copper, W/(m-K)
rho = 8960; % density of copper, kg/m^3
specificHeat = 386; % specific heat of copper, J/(kg-K)
thick = .01; % plate thickness in meters
stefanBoltz = 5.670373e-8; % Stefan-Boltzmann constant, W/(m^2-K^4)
hCoeff = 1; % Convection coefficient, W/(m^2-K)
% The ambient temperature is assumed to be 300 degrees-Kelvin.
ta = 300;
emiss = .5; % emissivity of the plate surface


numberOfPDE = 1;
model = createpde(numberOfPDE);

width = 1;
height = 1;


gdm = [3 4 0 width width 0 0 0 height height]';
g = decsg(gdm, 'S1', ('S1')');


geometryFromEdges(model,g);


figure;
pdegplot(model,'EdgeLabels','on');
axis([-.1 1.1 -.1 1.1]);
title 'Geometry With Edge Labels Displayed';

c = thick*k;

a = @(~,state) 2*hCoeff + 2*emiss*stefanBoltz*state.u.^3;
f = 2*hCoeff*ta + 2*emiss*stefanBoltz*ta^4;
d = thick*rho*specificHeat;
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);

applyBoundaryCondition(model,'dirichlet','Edge',1,'u',1000);

setInitialConditions(model,0);

hmax = .1; % element size
msh = generateMesh(model,'Hmax',hmax);
figure;
pdeplot(model);
axis equal
title 'Plate With Triangular Element Mesh'
xlabel 'X-coordinate, meters'
ylabel 'Y-coordinate, meters'

R = solvepde(model);
u = R.NodalSolution;
figure;
pdeplot(model,'XYData',u,'Contour','on','ColorMap','jet');
title 'Temperature In The Plate, Steady State Solution'
xlabel 'X-coordinate, meters'
ylabel 'Y-coordinate, meters'
axis equal
p = msh.Nodes;
plotAlongY(p,u,0);
title 'Temperature As a Function of the Y-Coordinate'
xlabel 'Y-coordinate, meters'
ylabel 'Temperature, degrees-Kelvin'
fprintf('Temperature at the top edge of the plate = %5.1f degrees-K\n', ...
  u(4));

specifyCoefficients(model,'m',0,'d',d,'c',c,'a',a,'f',f);
endTime = 5000;
tlist = 0:50:endTime;
numNodes = size(p,2);


u0(1:numNodes) = 300;

setInitialConditions(model,1000,'edge',1);

model.SolverOptions.RelativeTolerance = 1.0e-3;
model.SolverOptions.AbsoluteTolerance = 1.0e-4;


R = solvepde(model,tlist);
u = R.NodalSolution;
figure;
plot(tlist,u(3, :));
grid on
title 'Temperature Along the Top Edge of the Plate as a Function of Time'
xlabel 'Time, seconds'
ylabel 'Temperature, degrees-Kelvin'

figure;
pdeplot(model,'XYData',u(:,end),'Contour','on','ColorMap','jet');
title(sprintf('Temperature In The Plate, Transient Solution( %d seconds)\n', ...
  tlist(1,end)));
xlabel 'X-coordinate, meters'
ylabel 'Y-coordinate, meters'
axis equal;

fprintf('\nTemperature at the top edge(t=%5.1f secs)=%5.1f degrees-K\n', ...
  tlist(1,end), u(4,end));

