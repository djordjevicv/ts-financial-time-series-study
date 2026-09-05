clear;
clc;
close all;

addpath(genpath('src'));

run('scripts/load_transform.m');
run('scripts/nominal_vs_returns.m');
run('scripts/basic_statistics.m');
