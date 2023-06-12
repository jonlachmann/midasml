# Model w intercept works.

    Code
      res
    Output
      $lambda
      [1] 210.13057000   2.10130570   0.02101306
      
      $gamma
      [1] 0.8
      
      $cvm
      [1] 436.902128   1.942887   1.854322
      
      $cvsd
      [1] 65.6609218  0.5180706  0.4032474
      
      $cvupper
      [1] 502.563049   2.460958   2.257569
      
      $cvlower
      [1] 371.241206   1.424816   1.451074
      
      $nzero
      list()
      
      $name
      [1] "Single outcome sg-LASSO"
      
      $lamin
      $lamin$lambda.min
      [1] 0.02101306
      
      $lamin$lambda.1se
      [1] 2.101306
      
      
      $sgl.fit
      $b0
             s0        s1        s2 
       25.48788 -49.93209 -50.68629 
      
      $a0
      NULL
      
      $beta
      5 x 3 sparse Matrix of class "dgCMatrix"
         s0      s1           s2
      V1  0 0.33298 3.363098e-01
      V2  0 0.33298 3.363098e-01
      V3  0 0.33298 3.363098e-01
      V4  0 0.00000 3.697493e-17
      V5  0 0.00000 3.697493e-17
      
      $df
      [1] 0 3 5
      
      $dim
      [1] 5 3
      
      $lambda
      [1] 210.13057000   2.10130570   0.02101306
      
      $nf
      [1] 1
      
      $npasses
      [1] 6
      
      $jerr
      [1] 0
      
      $dimx
      [1] 50  5
      
      $call
      sglfit(x = x, y = y, gamma = gamma, nlambda = 3, method = "single", 
          lambda = lambda, gindex = gindex, intercept = ..1)
      
      attr(,"class")
      [1] "sglfit"  "sglpath"
      
      $cv.fit
      $cv.fit$lam.min
      $cv.fit$lam.min$b0
             s2 
      -50.68629 
      
      $cv.fit$lam.min$beta
                V1           V2           V3           V4           V5 
      3.363098e-01 3.363098e-01 3.363098e-01 3.697493e-17 3.697493e-17 
      
      
      $cv.fit$lam.1se
      $cv.fit$lam.1se$b0
             s1 
      -49.93209 
      
      $cv.fit$lam.1se$beta
           V1      V2      V3      V4      V5 
      0.33298 0.33298 0.33298 0.00000 0.00000 
      
      
      
      attr(,"class")
      [1] "reg.sgl"     "tscv.sglfit"

# Model w/o intercept works.

    Code
      res
    Output
      $lambda
      [1] 5522.7740255   55.2277403    0.5522774
      
      $gamma
      [1] 0.8
      
      $cvm
      [1] 939.871753  82.985915   1.865999
      
      $cvsd
      [1] 211.2955792  15.3543185   0.4381448
      
      $cvupper
      [1] 1151.167332   98.340234    2.304144
      
      $cvlower
      [1] 728.576173  67.631597   1.427854
      
      $nzero
      list()
      
      $name
      [1] "Single outcome sg-LASSO"
      
      $lamin
      $lamin$lambda.min
      [1] 0.5522774
      
      $lamin$lambda.1se
      [1] 0.5522774
      
      
      $sgl.fit
      $b0
      s0 s1 s2 
       0  0  0 
      
      $a0
      NULL
      
      $beta
      5 x 3 sparse Matrix of class "dgCMatrix"
         s0        s1            s2
      V1  0 0.5116853  1.0059188675
      V2  0 0.1680012  0.0000000000
      V3  0 0.0000000  0.0000000000
      V4  . .          .           
      V5  0 0.0000000 -0.0006957098
      
      $df
      [1] 0 2 2
      
      $dim
      [1] 5 3
      
      $lambda
      [1] 5522.7740255   55.2277403    0.5522774
      
      $nf
      [1] 0
      
      $npasses
      [1] 56
      
      $jerr
      [1] 0
      
      $dimx
      [1] 50  5
      
      $call
      sglfit(x = x, y = y, gamma = gamma, nlambda = 3, method = "single", 
          lambda = lambda, gindex = gindex, intercept = ..1)
      
      attr(,"class")
      [1] "sglfit"  "sglpath"
      
      $cv.fit
      $cv.fit$lam.min
      $cv.fit$lam.min$b0
      s2 
       0 
      
      $cv.fit$lam.min$beta
                 V1            V2            V3            V4            V5 
       1.0059188675  0.0000000000  0.0000000000  0.0000000000 -0.0006957098 
      
      
      $cv.fit$lam.1se
      $cv.fit$lam.1se$b0
      s2 
       0 
      
      $cv.fit$lam.1se$beta
                 V1            V2            V3            V4            V5 
       1.0059188675  0.0000000000  0.0000000000  0.0000000000 -0.0006957098 
      
      
      
      attr(,"class")
      [1] "reg.sgl"     "tscv.sglfit"

