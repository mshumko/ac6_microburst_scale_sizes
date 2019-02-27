# Versioning information for coincident_microburst_catalog_v#.csv
v0 = ('First run with all of the AC6 data. Many NaNs'
     ' appear when 0 coincident microbursts and 0 '
     'random microbursts were detected (0/0 division). '
     'In addition, infinites appear when one or the other '
     'spacecraft random occurance rate is 0. Example if '
     'AC6A observed microbursts during the pass, but not '
     'AC6B.') 
v1 = ('Second run with all of the AC6 data. Many NaNs'
     ' appear when 0 coincident microbursts and 0 '
     'random microbursts were detected (0/0 division).'
     ' Infinities have been replaced by square of the '
     'individual occurance rate.')   
v2 = ('CC_threshold = 0.8 (I think), I also changed the data '
      'structure. Now it is "pass_start", "total_sep",'
      '"uburst_rate", "chance_rate", "rate_ratio"') 
v3 = ('Same as v2, but with CC_threshold = 0.9.') 
v4 = ('Wavelet detector + peak_std variable added.') 
v5 = ('Same as v4 except used burst parameter.') 
v6 = ('Used burst parameter + first itertation to find peak width. '
      'A non-nan peak value is saved if there was a peak near the'
      ' center index +/- 1. Then the peak width taken at half of '
      'the prominence is saved.') 
    
