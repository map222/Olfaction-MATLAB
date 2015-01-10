# Olfaction-MATLAB
This repository contains sample code from my work at the University of Geneva on how the brain processes odors. I have included the paper I published from my time there (Patterson et al (2013)), as well as five scripts I used in writing the paper. Here is a brief description of each script.

MP_CalcPercentResponsiveCells10: This was a large script I used that called other scripts to analyze all my data. Of statistical relevance is the ANOVA that occurs in the middle of it, which compares the average firing rate of neurons before and after a breath. This ANOVA was used to throughout the paper, but for an example see the middle panels of Figures 1A-D, and Fig. 1E.

MP_KSTwoOdorBreaths: Neurons can change the timing of their firing in response to an odor. To look at changes in firing, I performed a Kolmogorov-Smirnov test on the distribution of spikes, contained in this script. To see changes in firing phase, see the bottom panels of Fig. 1A-D.

BreathSpikeDistanceBoot3: The two scripts above look at how individual cells respond to odor. To see how the population responded to odor, I created a population vector, and then calculated the distance between odors. The results of this script are shown in Fig. 2.

MP_PredictOdorBreath8: Another way to see what information was in the population is using a template matching algorithm. The results of this script are shown in Fig. 2C.

ReadBehavCSV: A large part of data science is simply parsing and cleaning data. This is an example script I used to parse .csv files.
