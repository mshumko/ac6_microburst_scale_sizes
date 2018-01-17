This repo will house files and intermediate data used to:
1) Catalogue microbursts seen with AC6 (> 35 keV dosimiters)
    using the burst parameter or the wavelet detection algorithms.
    - Catalogue will contain time of largest amplitude,
        amplitude, L, MLT, lat, lon, alt, (pitch angle?),
        in-track lag and distance, total distance, and maybe Kp and DST.
2) Merge the catalogues from the two units to determine if the
    microburst observation was spatial (curtain) or temporal 
    (flash).
    - If the microburst deemed to be both a spatial and temporal 
        event, flag it for the user and attemp to rectify this 
        ambiguity via:
        - Cross corralating the timeseries from both units to find
          the events with the highest correlation coefficient. 
    - Create separate catalogues for flashes and curtains containing
        the following information: times, amplitude, L, MLT, 
        spacecraft separation, spacecraft lag, Kp, and maybe DST.
3) Do statistical studies on the curtains and flashes to determine
    if they have different L-MLT distributions, or spatial
    scale size. 
    - Take care to correctly normalize the distributions.
4) Test the hypothesis if curtains are remnant microbursts by 
    checking if they, on average get gobbled up by the SAA 
    (There will be less seen east of the SAA, before the chorus 
    region in the dawn MLT sector).
    - Take care to correctly normalize the distributions.

Data tree structure is:
data
    microburst_catalogues - Results from step 1 (catalogues of bursts)
    flash_catalogues - Results from step 2 (catalogues of flashes)
    curtain_catalogues - Results from step 2 (catalogues of curtains)