# EnhancedPHAT
Generalized cross-correlation phase transform (GCC-PHAT) is a cornerstone for time difference of arrival (TDOA) estimation, direction-of-arrival (DOA) estimation, and source localization. It is simple, robust, and generally performs well. However, GCC assumes spatially uncorrelated noises, which is unlikely to be true in reality. For example, diffuse noises and reverberations are omnipresent, and their spatial coherence is known to be sinc(angular_frequency*mic_distance/sound_speed). Such spatial correlation introduces bias to the cross phase spectra used in GCC-PHAT. Fortunately, there exists a simple closed-form solution for compensating this bias just using the noise coherence function (please check the enclosed pdf file details).   

We have compared the standard GCC-PHAT, enhanced GCC-PHAT with bias removed cross phase spectra and wide band generalized eigenvalue decomposition multiple signal classification (GEVD-MUSIC) on a speech TDOA task. A set of typical comparison results are shown as below (please check enclosed data and code for details. The STFT window is designed by [this tool](https://sites.google.com/site/lixilinx/home/psmfb)).

![alt text](https://github.com/lixilinx/EnhancedPHAT/blob/master/phat_vs_music.png)

Enhanced GCC-PHAT is virtually as simple as the original GCC-PHAT, and yet, performs as well as GEVD-MUSIC, which is significantly more complicated and requires the estimation of noise power spectra density matrices.
