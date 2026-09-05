clear;
clc;
close all;

addpath(genpath('src'));

run('scripts/01_load_transform.m');
run('scripts/02_nominal_vs_returns.m');
run('scripts/03_basic_statistics.m');
