# metafrontier print output is stable

    Code
      print(fit)
    Output
      
      Metafrontier Model
      ------------------
      Method:           sfa 
      Metafrontier:     deterministic 
      Groups:           G1, G2 
      Total obs:        100 
        G1: 50 obs
        G2: 50 obs
      
      Group log-likelihoods:
        G1: 3.8358
        G2: -4.2613
      
      Mean TGR by group:
        G1: 1
        G2: 0.7922

# metafrontier summary output is stable

    Code
      summary(fit)
    Output
      
      Metafrontier Model Summary
      ==========================
      
      Call:
      metafrontier(formula = log_y ~ log_x1 + log_x2, data = sim$data, 
          group = "group", method = "sfa")
      
      Method:        sfa 
      Metafrontier:  deterministic 
      
      --- Group: G1 (n = 50) ---
                  Estimate Std. Error z value Pr(>|z|)    
      (Intercept)  0.76557    0.44512   1.720   0.0854 .  
      log_x1       0.52537    0.02328  22.564  < 2e-16 ***
      log_x2       0.16220    0.02329   6.964 3.32e-12 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      Log-likelihood: 3.8358 
      
      --- Group: G2 (n = 50) ---
                  Estimate Std. Error z value Pr(>|z|)    
      (Intercept)  0.44596    0.10769   4.141 3.45e-05 ***
      log_x1       0.51329    0.02282  22.490  < 2e-16 ***
      log_x2       0.20886    0.02453   8.513  < 2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      Log-likelihood: -4.2613 
      
      --- Metafrontier ---
                  Estimate
      (Intercept)   0.7656
      log_x1        0.5254
      log_x2        0.1622
      
      --- Efficiency Decomposition ---
       Group Mean_TE Mean_TGR Mean_TE_star
          G1  0.9959   1.0000       0.9959
          G2  0.7505   0.7922       0.5947
      
      --- Technology Gap Ratio Summary ---
       Group  N   Mean     SD    Min     Q1 Median     Q3    Max
          G1 50 1.0000 0.0000 1.0000 1.0000  1.000 1.0000 1.0000
          G2 50 0.7922 0.0589 0.6973 0.7501  0.784 0.8323 0.9162
      

# malmquist_meta print output is stable

    Code
      print(malm)
    Output
      
      Metafrontier Malmquist TFP Index
      ================================
      Orientation:  output 
      RTS:          crs 
      Groups:       G1, G2 
      Periods:      1 -> 2 
      Observations: 40 
      
      Mean decomposition (M* = TEC x TGC x TC*):
        MPI  = 1.295 
        TEC  = 3.315 
        TGC  = 2.129 
        TC*  = 0.3311 

# malmquist_meta summary output is stable

    Code
      summary(malm)
    Output
      
      Metafrontier Malmquist TFP Index Summary
      =========================================
      
      Call:
      malmquist_meta(formula = log_y ~ log_x1 + log_x2, data = pd, 
          group = "group", time = "time")
      
      Orientation: output 
      RTS:         crs 
      Groups:      G1, G2 
      Periods:     1 -> 2 
      
      --- Three-Way Decomposition by Group ---
      M* = TEC x TGC x TC*
      
      Group: G1 (n = 20 )
         MPI    TEC    TGC     TC 
      1.2811 1.1009 3.6293 0.3271 
      
      Group: G2 (n = 20 )
         MPI    TEC    TGC     TC 
      1.3121 5.5288 0.6295 0.3362 
      
      --- Technology Gap Ratios ---
      
      Group: G1 
        Mean TGR (from): 0.3377 
        Mean TGR (to):   0.9936 
        Mean TGC:        3.629 
      
      Group: G2 
        Mean TGR (from): 1 
        Mean TGR (to):   0.6295 
        Mean TGC:        0.6295 
      

# DEA metafrontier print output is stable

    Code
      print(fit)
    Output
      
      Metafrontier Model
      ------------------
      Method:           dea 
      Metafrontier:     deterministic 
      Groups:           G1, G2 
      Total obs:        60 
        G1: 30 obs
        G2: 30 obs
      
      Mean TGR by group:
        G1: 1
        G2: 0.8398

