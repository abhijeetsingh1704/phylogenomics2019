      seqfile = seqfile.txt         * sequence data filename
      outfile = results_0.001.txt   * main result file name

        noisy = 9      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = -2     * -2:pairwise

      seqtype = 1      * 1:codons
    CodonFreq = 3      * 0:equal, 1:F1X4, 2:F3X4, 3:F61
        model = 0      *
      NSsites = 0      * 
        icode = 0      * 0:universal code

    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * initial or fixed kappa

    fix_omega = 1      * 1:omega fixed, 0:omega to be estimated 
        omega = 0.001  * 1st fixed omega value [change this]
       
       * EXCERCISE 1
       *alternate fixed omega values
       *omega = 0.005  * 2nd fixed value 
       *omega = 0.01   * 3rd fixed value
       *omega = 0.05   * 4th fixed value
       *omega = 0.10   * 5th fixed value
       *omega = 0.20   * 6th fixed value
       *omega = 0.40   * 7th fixed value
       *omega = 0.80   * 8th fixed value
       *omega = 1.60   * 9th fixed value
       *omega = 2.00   * 10th fixed value
