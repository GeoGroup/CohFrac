clear;clc;close all;
% 
% Young_modulus=17e9;
% Poisson_ratio=0.2;
% Fracture_energy=50;
% Tensile_strength=9e5;

Young_modulus=15e9;
Poisson_ratio=0.25;
Fracture_energy=3000;
Tensile_strength=2e6;

cz_length=((9*pi*Young_modulus)/(32*(1-Poisson_ratio^2)))*Fracture_energy/Tensile_strength^2
