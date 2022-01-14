The code for "Channel Estimation and Projection for RIS-assisted MIMO Using Zadoff-Chu Sequences"

The relationship between main function and figures in paper:
Fig3: main_H_CRB.m
Fig4: main_H.m
Fig5: main_Hd.m
Fig6: main_H_K

Introduction for functions:
multi_signal_sameBeta：generate channel from UE to RIS；
multi_signal_Hd：generate channel from UE to BS；
Hd_est_fft_sage：Coarse Solution to (88) Using FFTs；
Hd_est_newton_sage：Refined Channel Estimation for Hd Using Newton's Method;
allD_est_fft_sage：Coarse Solution to (39) Using FFTs；
allD_est_newton_sage：Refined Channel Estimation for H Using Newton’s Method
