# fullSKT
Related to the paper "On the influence of cross-diffusion in pattern formation" by M. Breden, C. Kuehn, C. Soresina

We selected three figures of the paper, namely Figures 1, 9, and 15. The provided code in each folder generates these figures. In particular:
- Figure 1: triangular case, weak competition.
- Figure 9: non-triangular case, strong competition.
- Figure 15: cross-diffusion and self-diffusion, weak competition.

The four parameter sets used in the paper are:
1) r1=5; r2=2; a1=3; a2=3; b1=1; b2=1;
2) r1=2; r2=5; a1=1; a2=1; b1=0.5; b2=3;
3) r1=15/2; r2=16/7; a1=4; a2=2; b1=6; b2=1;
4) r1=5; r2=5; a1=2; a2=3; b1=5; b2=4;
If you change the parameter set, you should also check and adapt the continuation parameter in the -init.m file.

*Please note*: the main file 'SKT1DX.m' sets pde2path, if you specify the folder path (lines 10,11).
