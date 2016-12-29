MODIFICATIONS FROM THE PREVIOUS VERSION:
————————————————————————————————————————
- New formula for the borders Shannon-Nyquist Exclusion zone
=> changes in freq_analysis_prelims and timefreq_analysis_prelims
- Add new options for the tapering windows (frequency analysis)
- in ‘carma_params’ (in Wavepal.py), when signif_level_type=« a », there also a first round to estimate « mylength », and a second round generating nmcmc*mylength samples. It is now exactly like when « n » is in signif_level_type. 