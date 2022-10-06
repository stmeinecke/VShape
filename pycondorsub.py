import numpy as np
import os

DIR="/net/granat/users2/meinecke/Vshape/"
PROG="Vshape"
RUN="5_JG_tau_long_int_delta_10000_100x200"

if(os.path.exists(PROG)):
  print "program " + PROG + " found"
  
  PATH = DIR+"run"+RUN
  
  if( os.path.exists(PATH)  ):
    print "Error: directory " + PATH + "already exists"
  
  else:
     
    os.mkdir(PATH)
    os.mkdir(PATH+"/data")
    os.mkdir(PATH+"/data/TS")
    os.mkdir(PATH+"/data/netgain")
    os.mkdir(PATH+"/data/bin")
    os.mkdir(PATH+"/data/powerSpec")
    os.mkdir(PATH+"/data/IntAC")
    os.mkdir(PATH+"/data/LTTJ")
    os.mkdir(PATH+"/data/LTTJ/powerspecs")
    
    print "Directory " + PATH + " has been created"
    
    argfile = open('pyargfile','w')
    
    

      
    ##JG JQ
    #MEM=2
    ##for J_Q in np.linspace(0,-100,201, endpoint=True):
    #for J_Q in [-65.0]:
      #argfile.write("-sweep_J_G_with_J_Q -J_Q " + "{:1.3f}".format(J_Q) +" -sIntTime 5000 -sOutTime 500 -sStart 0.0 -sEnd 4.0 -sSteps 400")
      #argfile.write(" -noise -nStr 1E-1")
      ##argfile.write(" -T01 0.26666667")
      ##argfile.write(" -wpowerSpec")
      #argfile.write(" -wIntAC")
      #argfile.write(" -wExtrema")
      #argfile.write(" -wExtrema_all")
      ##argfile.write(" -wDeltaT")
      #argfile.write(" -wNetgain")
      ##argfile.write(" -wNoSweep")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_FML_J1.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML2LEI_J2.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML4LEI_J3.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML6LEI_J5.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC2_J2.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC3_J4.4.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC4_J6.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC5_J5.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC6_J6.0.bin")
      #argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML5_J4.5.bin")
      #argfile.write("\n")
      
    ##JG Q_0
    #MEM=1
    #for J_Q in np.linspace(0,-0.05,101, endpoint=True):
      #argfile.write("-sweep_J_G_with_Q_0 -Q_0 " + "{:1.6f}".format(J_Q) +" -sIntTime 2000 -sOutTime 100 -sStart 0.0 -sEnd 0.3 -sSteps 100")
      #argfile.write(" -noise -nStr 1E-1")
      ##argfile.write(" -wpowerSpec")
      ##argfile.write(" -wIntAC")
      #argfile.write(" -wExtrema")
      ##argfile.write(" -wExtrema_all")
      ##argfile.write(" -wDeltaT")
      ##argfile.write(" -wNetgain")
      ##argfile.write(" -wNoSweep")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_FML_J1.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML2LEI_J2.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML4LEI_J3.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML6LEI_J5.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC2_J2.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC3_J4.4.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC4_J6.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC5_J5.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC6_J6.0.bin")
      #argfile.write("\n")
    
    #JG T01
    MEM=4
    for T01 in np.linspace(0,0.5,101, endpoint=True):
      argfile.write("-sweep_J_G_with_T01 -T01 " + "{:1.6f}".format(T01) +" -sIntTime 80000 -sOutTime 1000 -sStart 0.015 -sEnd 0.13 -sSteps 200")
      #argfile.write(" -noise -nStr 1E-1")
      argfile.write(" -dt 4e-5 -nout 2")
      argfile.write(" -Q_0 -0.03 -gamma_Q 200 -delta 10000")
      #argfile.write(" -alpha_G 2.0 -alpha_Q 2.0")
      #argfile.write(" -wpowerSpec")
      #argfile.write(" -wIntAC")
      #argfile.write(" -wExtrema")
      #argfile.write(" -wExtrema_all")
      #argfile.write(" -wDeltaT")
      #argfile.write(" -wNetgain")
      #argfile.write(" -wNoSweep")
      #argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_FML.bin")
      argfile.write("\n")
      
    ##T01 JG
    #MEM=2
    ##for JG in np.linspace(0,4.0,201, endpoint=True):
    #for JG in [1.48,1.49,1.5,1.51,1.52]:
      #argfile.write(" -sweep_TO1_with_J_G -J_G " + "{:1.6f}".format(JG) +" -sIntTime 1000 -sOutTime 100 -sStart 0.0 -sEnd 0.5 -sSteps 200")
      #argfile.write(" -noise -nStr 1E-1")
      ##argfile.write(" -wpowerSpec")
      ##argfile.write(" -wIntAC")
      #argfile.write(" -wExtrema")
      ##argfile.write(" -wExtrema_all")
      ##argfile.write(" -wDeltaT")
      #argfile.write(" -wNetgain")
      ##argfile.write(" -wNoSweep")
      #argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_FML_J1.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML2LEI_J2.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML4LEI_J3.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML6LEI_J5.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC2_J2.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC3_J4.4.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC4_J6.0.bin")
      #argfile.write("\n")
    
      
    ##alpah_G alpha_Q
    #MEM=2
    #for aQ in np.linspace(0,6.0,201, endpoint=True):
      #argfile.write("-sweep_alpha_G_with_alpha_Q -alpha_Q " + "{:1.6f}".format(aQ) +" -sIntTime 2000 -sOutTime 500 -sStart 0.0 -sEnd 6.0 -sSteps 200")
      #argfile.write(" -noise -nStr 1E-1")
      #argfile.write(" -J_G 1.0")
      ##argfile.write(" -alpha_G 5.0 -alpha_Q 5.0")
      ##argfile.write(" -wpowerSpec")
      ##argfile.write(" -wIntAC")
      ##argfile.write(" -wExtrema")
      ##argfile.write(" -wExtrema_all")
      ##argfile.write(" -wDeltaT")
      ##argfile.write(" -wNetgain")
      ##argfile.write(" -wNoSweep")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_FML_J1.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML2LEI_J2.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML4LEI_J3.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_HML6LEI_J5.0.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC2_J2.5.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC3_J4.4.bin")
      ##argfile.write(" -wNoSweep -loadHist -histFile " + DIR + "Xhist_PC4_J6.0.bin")
      #argfile.write("\n")
      
      
    ##JG nStr
    #MEM=1
    #for nStr in np.geomspace(1e-50,1e3,201, endpoint=True):
      #argfile.write("-sweep_J_G_with_nStr -nStr " + "{:1.6e}".format(nStr) +" -sIntTime 5000 -sOutTime 200 -sStart 0.0 -sEnd 6.0 -sSteps 300")
      ##argfile.write(" -wpowerSpec -wIntAC -wExtrema -wDeltaT")
      #argfile.write(" -noise")
      #argfile.write(" -strSuffix "+"{:1.6e}".format(nStr))
      ##argfile.write(" -wNoSweep")
      ##argfile.write(" -loadHist -histFile " + DIR + "Xhist_stdHML3.bin")
      #argfile.write("\n")
      
      
    ##JG delta
    #MEM=8
    #for nStr in np.geomspace(50,10000,201, endpoint=True):
      #argfile.write("-sweep_J_G_with_delta -delta " + "{:1.6e}".format(nStr) +" -sIntTime 5000 -sOutTime 400 -sStart 0.0 -sEnd 4.0 -sSteps 200")
      #argfile.write(" -calc_dt")
      ##argfile.write(" -wpowerSpec")
      ##argfile.write(" -wIntAC")
      ##argfile.write(" -wExtrema")
      ##argfile.write(" wDeltaT")
      ##argfile.write(" -wExtrema_all")
      ##argfile.write(" -noise -nStr 1E-1")
      ##argfile.write(" -wNoSweep")
      #argfile.write("\n")
      
    
    ##JG delta - low gammma - low noise
    #MEM=4
    #for nStr in np.geomspace(50,10000,201, endpoint=True)[100:]:
      #argfile.write("-sweep_J_G_with_delta -delta " + "{:1.6e}".format(nStr) +" -sIntTime 10000 -sOutTime 200 -sStart 0.0 -sEnd 4.0 -sSteps 200")
      #argfile.write(" -calc_dt")
      ##argfile.write(" -wpowerSpec")
      ##argfile.write(" -wIntAC")
      ##argfile.write(" -wExtrema")
      ##argfile.write(" wDeltaT")
      ##argfile.write(" -wExtrema_all")
      #argfile.write(" -noise -nStr 1E-5")
      #argfile.write(" -wNoSweep")
      #argfile.write("\n")
    
    #alphas = np.linspace(0.0,6.0,201)
    #MEM = 2
    #for k in range(alphas.size):
      #argfile.write("-TS -netGain -outTime 100 -nout 1 -noTSout -noise -intTime 5000") 
      #argfile.write(" -alpha_G " + "{:1.6f}".format(alphas[k]))
      #argfile.write(" -alpha_Q " + "{:1.6f}".format(alphas[k]))
      #argfile.write(" -strSuffix _aG_" + "{:1.6f}".format(alphas[k]) + "_aQ_" + "{:1.6f}".format(alphas[k]))       
      #argfile.write(" -loadHist -histFile " + DIR + "Xhist_FML_J1.0.bin")
      #argfile.write("\n")
      
 
  
    argfile.close()
    print "pyargfile created"
    
    substr = "qsub -mem " + str(MEM) + " -m n -speed 4 -w " + PATH + " -argfile pyargfile " + PROG
    print "Submit string: "
    print substr
    
    os.system( substr )

else:
  print "Error: program not found!"

