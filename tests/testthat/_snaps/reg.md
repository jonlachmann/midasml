# Model w intercept works.

    Code
      res
    Output
      $lambda
      [1] 207.96733826   2.07967338   0.02079673
      
      $gamma
      [1] 0.8
      
      $cvm
      [1] 421.8696869   0.7125367   0.7995840
      
      $cvsd
      [1] 58.4168818  0.1824065  0.2239527
      
      $cvupper
      [1] 480.2865687   0.8949433   1.0235367
      
      $cvlower
      [1] 363.4528051   0.5301302   0.5756313
      
      $nzero
      list()
      
      $name
      [1] "Single outcome sg-LASSO"
      
      $lamin
      $lamin$lambda.min
      [1] 2.079673
      
      $lamin$lambda.1se
      [1] 2.079673
      
      
      $sgl.fit
      $b0
             s0        s1        s2 
       25.53440 -49.10914 -49.85558 
      
      $a0
      NULL
      
      $beta
      5 x 3 sparse Matrix of class "dgCMatrix"
         s0           s1           s2
      V1  0 3.295521e-01 3.328476e-01
      V2  0 3.295521e-01 3.328476e-01
      V3  0 3.295521e-01 3.328476e-01
      V4  0 5.676626e-17 4.378542e-17
      V5  0 5.676626e-17 4.378542e-17
      
      $df
      [1] 0 5 5
      
      $dim
      [1] 5 3
      
      $lambda
      [1] 207.96733826   2.07967338   0.02079673
      
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
             s1 
      -49.10914 
      
      $cv.fit$lam.min$beta
                V1           V2           V3           V4           V5 
      3.295521e-01 3.295521e-01 3.295521e-01 5.676626e-17 5.676626e-17 
      
      
      $cv.fit$lam.1se
      $cv.fit$lam.1se$b0
             s1 
      -49.10914 
      
      $cv.fit$lam.1se$beta
                V1           V2           V3           V4           V5 
      3.295521e-01 3.295521e-01 3.295521e-01 5.676626e-17 5.676626e-17 
      
      
      
      attr(,"class")
      [1] "reg.sgl"     "tscv.sglfit"

# Model w/o intercept works.

    Code
      res
    Output
      $lambda
      [1] 5530.4174529   55.3041745    0.5530417
      
      $gamma
      [1] 0.8
      
      $cvm
      [1] 997.9771069  77.0119420   0.7486296
      
      $cvsd
      [1] 205.849974  12.874938   0.198751
      
      $cvupper
      [1] 1203.8270807   89.8868804    0.9473807
      
      $cvlower
      [1] 792.1271331  64.1370035   0.5498786
      
      $nzero
      list()
      
      $name
      [1] "Single outcome sg-LASSO"
      
      $lamin
      $lamin$lambda.min
      [1] 0.5530417
      
      $lamin$lambda.1se
      [1] 0.5530417
      
      
      $sgl.fit
      $b0
      s0 s1 s2 
       0  0  0 
      
      $a0
      NULL
      
      $beta
      5 x 3 sparse Matrix of class "dgCMatrix"
         s0        s1          s2
      V1  0 0.4970415 0.991879761
      V2  0 0.1734487 0.003248937
      V3  . .         .          
      V4  . .         .          
      V5  . .         .          
      
      $df
      [1] 0 2 2
      
      $dim
      [1] 5 3
      
      $lambda
      [1] 5530.4174529   55.3041745    0.5530417
      
      $nf
      [1] 0
      
      $npasses
      [1] 5
      
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
               V1          V2          V3          V4          V5 
      0.991879761 0.003248937 0.000000000 0.000000000 0.000000000 
      
      
      $cv.fit$lam.1se
      $cv.fit$lam.1se$b0
      s2 
       0 
      
      $cv.fit$lam.1se$beta
               V1          V2          V3          V4          V5 
      0.991879761 0.003248937 0.000000000 0.000000000 0.000000000 
      
      
      
      attr(,"class")
      [1] "reg.sgl"     "tscv.sglfit"

