# Enhanced PHAT

These supplementary materials are prepared for paper:

Xi-Lin Li, On correcting the phase bias of GCC in spatially correlated noise fields, Signal Processing, 2020.

Generalized cross-correlation (GCC) is widely used for time difference of arrival (TDOA) estimation, direction-of-arrival (DOA) estimation, and source localization. It is simple, robust, and generally performs well. However, GCC assumes spatially uncorrelated noises, which is unlikely to be true in reality. This paper proposes a simple method to remove the phase bias in GCC by assuming only the knowledge of noise coherence function.  

### On the microphone recordings

The acoustic data have been recorded in a typical conference room with size about (4, 5, 2.5) in meters. The microphone array is put on a desk located in the middle of room. The speaker is close to a wall. The distance between microphone array and speaker is about 2.5 meters. The distance between the two microphones is 0.11 meter. The noise coherence function can be predicted as   

sinc(angular_frequency*mic_distance/sound_speed)

where mic_distance = 0.11 meter. Note that if the microphones are put inside certain chambers or tubes, this mic_distance is supposed to be slightly larger than the actual physical distance. 

We have done two sets of recordings, i.e., broadside and endfire. For the broadside case, speech arrives at the array from front direction; for the endfire case, speech arrives from the side direction.

The STFT window is designed by [this method](https://ieeexplore.ieee.org/document/8304771). Below is an estimated noise coherence function sample. We can see that the predicted noise coherence function is quite accurate, except at frequencies round 1600 Hz. The imaginary part of noise coherence function apparently deviates from 0 there. This possibly is caused by the air conditioner noise, which is somewhat directional. Anyway, we use the theoretically predicted noise coherence function for enhanced GCC-PHAT in the paper, and it works well.     

![alt text](https://github.com/lixilinx/EnhancedPHAT/blob/master/coherence.png)

### Further details on the derivation of ML solutions

Those are given in file ML_solution_derivation.pdf.
