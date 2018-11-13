# Versioning information for microburst_catalog_v#.csv
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
v2 = ('Same as v1, except that I changed the data '
      'structure. Now it is "pass_start", "total_sep",'
      '"uburst_rate", "chance_rate", "rate_ratio"') 
    