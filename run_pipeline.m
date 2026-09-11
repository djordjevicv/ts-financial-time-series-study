clear;
clc;
close all;

addpath(genpath('src'));

run('scripts/load_transform.m');
run('scripts/nominal_vs_returns.m');
run('scripts/basic_statistics.m');
run('scripts/ar_model.m');
run('scripts/arch_effect.m');

if hasARCH
    run('scripts/garch_model.m');
end

run('scripts/forecasting.m');
