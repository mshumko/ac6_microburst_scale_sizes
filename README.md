# OVERVIEW

This repo contains the code and the processed data used to
estimate the microburst scale size distribution using the
pair of AC-6 CubeSats. The basic idea is to build up a 
statistical distribution of microbursts observed 
coincidently by both AC-6 units and save their separation
for each coincident event.

Then we remove noisy detections (radio noise, and other 
dosimiter artifacts), and calculate the spatial
normalization distribution. This normalization distribution
is calculated by finding the number of 10 Hz data time stamps
at which A) both units were taking 10 Hz data, and B) the data
quality flag is 0 (good data). 

Then the LEO spatial scale size distribution is the histogram
of the number of coincident microbursts as a function of 
spacecraft separation, divided by the observation time as a 
function of separation (the normalization distributions). In 
theory, the spacecraft separation at which no more coincident
microbursts are observed is the largest microburst scale size.

# REPO STRUCTURE
├── data/ **Contains the processed data files** <br />
│   ├── catalogue_version_notes.txt **User info on different catalog versions.** <br />
│   ├── curtain_catalogues/ **Folder with spatial event catalogs** <br />
│   ├── coincident_microbursts_catalogues/ **Folder with temporal event catalogs** <br />
│   ├── microburst_catalogues/ **Folder with microburst event catalogs (uncombined from both units)** <br />
│   └── norm/ **Normalization csv files** <br />
|   └── plots/ **Contains validation plots (not pushed)** <br />
├── docs/ **In-depth repo documentation. Manuscript will be added here** <br />
├── flash-curtain-determination/ <br />
│   ├── append_goemag_indicies.py <br />
│   ├── plotMatchedBursts.py
│   ├── replace_error_sep_lags.py
│   ├── sort_curtains_flashes.py
│   └── test_sort_curtains_flashes.py
├── logs
│   └── microburst_detection.log
├── microburst-detection
│   ├── merge_daily_data.py
│   ├── microburst_detection.py
│   └── microburst_detection_wrapper.py
├── README.md **This file**
└── stats **Folder for statistical analysis of the catalogs**
    ├── ac6_mission_separation.py
    ├── equatorial_scale_sizes.py
    ├── leo_scale_sizes.py
    ├── L_MLT_distribution.py
    ├── L_scale_size_wrapper.py
    ├── microburst_browser.py
    ├── norm.py
    ├── old
    │   ├── bar_test.py
    │   ├── dist_norm.py
    │   ├── gaus_fit.py
    │   ├── global_distribution.py
    │   ├── L_MLT_distribution.py
    │   ├── l_mlt_norm.py
    │   ├── MLT_lon_distribution.py
    │   ├── normalization
    │   │   ├── dist_norm.npy
    │   │   ├── dist_norm.txt
    │   │   ├── dist_vals.npy
    │   │   ├── L_edges_norm.npy
    │   │   ├── L_MLT_sec_norm.npy
    │   │   └── MLT_edges_norm.npy
    │   ├── scale_size_norm.py
    │   └── scale_sizes.py

# HOW TO RUN THE ENTIRE PIPELINE


MORE 
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