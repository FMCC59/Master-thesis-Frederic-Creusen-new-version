C Hardening Curves for VUMAT
C ==================================================================
C
      PARAMETER (NCC=33,NCT=19,NCS=20)
      DIMENSION YCC(NCC,2),YCT(NCT,2),YCS(NCS,2)
C
C     Compression hardening curve
      data (YCC(i,1), i=1,NCC) /
     *  0.000000000,
     *  0.002239247,
     *  0.005701262,
     *  0.009938009,
     *  0.014639938,
     *  0.020521984,
     *  0.026802194,
     *  0.033066911,
     *  0.040204922,
     *  0.047022298,
     *  0.054138719,
     *  0.061194225,
     *  0.067668820,
     *  0.074486195,
     *  0.080994020,
     *  0.087196699,
     *  0.093997463,
     *  0.100377351,
     *  0.106143105,
     *  0.112528532,
     *  0.118792128,
     *  0.125338153,
     *  0.131236249,
     *  0.137665982,
     *  0.143741287,
     *  0.149308245,
     *  0.155947846,
     *  0.161713038,
     *  0.167645488,
     *  0.173697521,
     *  0.179706937,
     *  0.185461052,
     *  10.000000000 /
      data (YCC(i,2), i=1,NCC) /
     *  130.253256600,
     *  148.302697800,
     *  159.776678700,
     *  169.066195700,
     *  176.502209500,
     *  178.333636600,
     *  179.481028700,
     *  178.885275400,
     *  177.517220300,
     *  176.524285600,
     *  174.339814100,
     *  172.398059200,
     *  171.868489000,
     *  170.875561800,
     *  170.213587800,
     *  171.669902100,
     *  170.743158000,
     *  171.493383700,
     *  173.788175400,
     *  174.516332500,
     *  175.729922700,
     *  175.818181900,
     *  178.488094000,
     *  179.039717700,
     *  181.003548700,
     *  183.187997500,
     *  183.805842000,
     *  187.005309200,
     *  187.733466300,
     *  191.594900000,
     *  192.918825500,
     *  196.162429900,
     *  200.000000000 /
C
C     Tension hardening curve
      data (YCT(i,1), i=1,NCT) /
     *  0.000000000,
     *  0.000018161,
     *  0.000034868,
     *  0.000086233,
     *  0.000163815,
     *  0.000240308,
     *  0.000317509,
     *  0.000395132,
     *  0.000493514,
     *  0.000636065,
     *  0.000790903,
     *  0.000952196,
     *  0.001127912,
     *  0.001319586,
     *  0.001545202,
     *  0.001774464,
     *  0.002049382,
     *  0.002339382,
     *  10.000000000 /
      data (YCT(i,2), i=1,NCT) /
     *  16.446230840,
     *  19.085320600,
     *  21.866976990,
     *  24.695045290,
     *  27.430281850,
     *  30.208618660,
     *  33.063210800,
     *  36.083573820,
     *  39.047578640,
     *  41.878964260,
     *  44.829704330,
     *  47.813598340,
     *  50.843909930,
     *  53.778073590,
     *  56.718865090,
     *  59.878479690,
     *  62.862378230,
     *  65.876108880,
     *  66.000000000 /
C
C     Shear hardening curve
      data (YCS(i,1), i=1,NCS) /
     *  0.000000000,
     *  0.000216694,
     *  0.000466873,
     *  0.000519947,
     *  0.000722971,
     *  0.000919799,
     *  0.001199673,
     *  0.001394640,
     *  0.001703889,
     *  0.001956767,
     *  0.002368569,
     *  0.002919021,
     *  0.003338433,
     *  0.003838918,
     *  0.004393124,
     *  0.005186539,
     *  0.005717747,
     *  0.006303682,
     *  0.007170760,
     *  10.000000000 /
      data (YCS(i,2), i=1,NCS) /
     *  19.769928170,
     *  22.226491080,
     *  24.587649400,
     *  27.401952660,
     *  29.825137620,
     *  32.410487220,
     *  34.805040040,
     *  37.385619000,
     *  39.794475480,
     *  42.375058580,
     *  44.564494410,
     *  46.515443570,
     *  48.709645900,
     *  50.755987240,
     *  52.745089140,
     *  54.242868520,
     *  56.274902070,
     *  58.197219090,
     *  59.766553990,
     *  60.000000000 /
C