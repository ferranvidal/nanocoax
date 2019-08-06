d0=fileparts([pwd,filesep]);
addpath([d0,'/meshes']);
addpath([d0,'/supplemental_functions']);

clear all; close all; clc;

%% Run one example
G = 1; % in nanometers
b = 1.5; % factor than multiplies beta2, set to 0 to run LRA
lambda = 5.932; % in microns

[mesh,setup,Tp,FE,EDG,HDG,JDG,RDG] = run_example(G,b,lambda);

%% Visualize |Ex| on a 2D slice of the mesh, similarly for other fields
plotmesh = 1; % plotting the mesh (choose 0 or 1)

ztop = mesh.zfilm(2); % upper interface gold air
figure()
faceplot(mesh,abs(EDG(:,1,:)),ztop,plotmesh)

zmid = mesh.zgoldmid; % middle of gold layer
figure()
faceplot(mesh,abs(EDG(:,1,:)),zmid,plotmesh)

save(['dis_',num2str(G),num2str(b),'t.mat'],'mesh','setup','Tp','FE','EDG','HDG')