# Sex differences in functional modulation of microglia by early-life physical stress in a rat model of chronic primary low back pain

Early-life physical stress reprograms microglia and shapes adult pain vulnerability in a sex-specific manner. 
This repo contains the analysis pipeline (R) to reproduce the paper’s key figures: microglia state classification, dimensionality reduction, and group comparisons.

# What happened in this study?

Chronic primary low back pain is a leading cause of disability.
Early-life stress (ELS) is a strong risk factor that can leave long-lasting immune-brain imprints.
Microglia are central to neuroimmune modulation of pain, and sex differences are increasingly recognized as critical to understanding mechanisms and tailoring therapies.
Hence, in this study the funcional modulation of microglia was studied according to its phenotype in the dorsal horn spinal cord region in rats. The microglia were classified as: Surveillant, Primed, Activated. 
The proportion of microglia each state was identified to find out the response of ELS on the predisposition of chronic primary low back pain in both the sexes.

### Description
The microglia were segmented using the 3DMorph MATLAB script.

<img width="473" height="226" alt="image" src="https://github.com/user-attachments/assets/b4a110e1-0c54-416e-9463-ad82aa7c5d0d" />

An unsupervised PAM clustering method was used for the clustering of the cells accorsing to teh extracted features into 3 main clusters:
<img width="473" height="294" alt="image" src="https://github.com/user-attachments/assets/6a66fa71-bcc5-4c40-9248-fb3e685d2c6f" />
(Details in the paper)

As a results males were found to posses higher proportion of primed and activated microglia and females had higher surveillant state of microglia.
