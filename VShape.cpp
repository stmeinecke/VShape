// #define REALTIME
#define SHIFTEDTIME

#include "incl/global.cpp"
#include "incl/get_params.cc"
#include "incl/get_specs_sm_v6.cc"
#include "incl/get_correlations_sm.cc"
#include "VShape_parameters.cpp"
#include "VShape_vars.cpp"
#include "incl/vars_vec_v9.cpp"
#include "incl/DDEintegrator_v5.cpp"
#include "incl/evalMLTS_v14.cpp"
#include "VShape_eqs.cpp"



using namespace std;

int main(int argc, char* argv[]){
  
//   //initialize exp lookup table
//   sm::init_expf();
  
    std::string str_suffix = katana::getCmdOption(argv, argv+argc, "-strSuffix" , "");
  
  //noise
  katana::Seed s;
  s.set_seed();
  
  //solver
  DDEintegrator DDEsolver;
  //solver parameters
  unsigned long long int tn = 0;
  double intTime = katana::getCmdOption(argv, argv+argc, "-intTime" , 300.0);
  double dt = katana::getCmdOption(argv, argv+argc, "-dt" , 4e-5); //should be roughly 1/(5*delta) 

  
  //thresholds and triggers for time series analysis, i.e. for pulse detection
  double T_maxPW = katana::getCmdOption(argv, argv+argc, "-T_maxPW" , 0.3); //max pulse width search radius in units of T0
  double T_dCountTol = katana::getCmdOption(argv, argv+argc, "-T_dCountTol" , 0.05); // tolerance for detecting unique maxima
  double T_relMax = katana::getCmdOption(argv, argv+argc, "-T_relMax" , 0.01); //threshold for maxima = pulses
  double T_TVPD_trig = katana::getCmdOption(argv, argv+argc, "-T_TVPD_trig" , 0.05); //trigger value for pulse detection
  double T_TVPD_bound = katana::getCmdOption(argv, argv+argc, "-T_TVPD_bound" , 0.001); //pulse boundaries for weighted pulse position/pulse area
  
  //system parameters
  parameters p; 
  //set default parameters:
  

  //Time rescaled with T
  p.TR = katana::getCmdOption(argv, argv+argc, "-TR" , 1.0);
  p.T = katana::getCmdOption(argv, argv+argc, "-T" , 1.0);
  p.T01 = katana::getCmdOption(argv, argv+argc, "-T01" , 0.25);
  p.T12 = katana::getCmdOption(argv, argv+argc, "-T12" , 0.25);
  p.alpha_G = katana::getCmdOption(argv, argv+argc, "-alpha_G" , 0.0);
  p.alpha_Q = katana::getCmdOption(argv, argv+argc, "-alpha_Q" , 0.0);
  p.kappa = katana::getCmdOption(argv, argv+argc, "-kappa" , 0.99);
  p.delta = katana::getCmdOption(argv, argv+argc, "-delta" , 10000)/p.TR;
  p.gamma_Q = katana::getCmdOption(argv, argv+argc, "-gamma_Q" , 200)/p.TR; 
  p.gamma_G = katana::getCmdOption(argv, argv+argc, "-gamma_G" , 0.625)/p.TR;
  p.r_s = katana::getCmdOption(argv, argv+argc, "-r_s" , 2.0);
  p.J_Q = katana::getCmdOption(argv, argv+argc, "-J_Q" , -6)/p.TR;
  p.J_G = katana::getCmdOption(argv, argv+argc, "-J_G" , 0.04)/p.TR;
  p.omega = katana::getCmdOption(argv, argv+argc, "-omega" , 0.0)/p.TR;
  
  p.Q_0 = katana::getCmdOption(argv, argv+argc, "-Q_0" , -0.03)/p.TR;
  bool set_Q_0 = katana::getCmdOption_bool(argv, argv+argc, "-Q_0" , false);
  
  p.nStr = katana::getCmdOption(argv, argv+argc, "-nStr" , 1E-1);
  
  
  bool bool_calc_dt = katana::getCmdOption_bool(argv, argv+argc, "-calc_dt" , false);
  if(bool_calc_dt == true){
    dt = 1.0/(5.0*p.delta);
    cout << "Time step calculated: dt = " << dt << endl;
  }
  

  auto computeParas = [&](){
    p.T12 = 0.5*p.T - p.T01;
    p.SqrtNoiseStr = sqrt( p.nStr / dt );
    if(set_Q_0) p.J_Q = p.Q_0 * p.gamma_Q;
  };
  computeParas();
  
  //dynamical variables with delay
  double max_delay = katana::getCmdOption(argv, argv+argc, "-tau_max" , p.T);
  vars_vec_wdX Xhist(max_delay,dt);
  cout << "History initilized with max_delay / dt: " << (max_delay)/dt << endl;

  //   double intpart;
//   if(modf(max_delay/dt, &intpart) > 1E-3){
//     cout << "Caution: max_tau is not an integer multiple of dt" << endl;
//     Xhist.dt_divides_tau = false;
//   }
//   else Xhist.dt_divides_tau = true;
  
  //////////////////////////////////////////
  //noise
  //////////////////////////////////////////
  
  p.SqrtNoiseStr = sqrt( p.nStr * dt );
  
  std::function<void (vars*, vars_vec*, parameters*)> noise;
  
  std::function<void (vars*, vars_vec*, parameters*)> GWnoise = [&](vars* X, vars_vec* Xhist, parameters* p){ 
      X->E += p->SqrtNoiseStr*katana::gwNoise(); 
  };
    
  if(katana::getCmdOption_bool(argv, argv+argc, "-noNoise" , false)){
    noise = after_step_empty;
  }
  else{
    noise = GWnoise;
  }
 
  
  
  //////////////////////////////////////////
  //outputfunctions
  //////////////////////////////////////////
  
  double outTime = katana::getCmdOption(argv, argv+argc, "-outTime" , 100.0);
  if(outTime < 0.0) outTime = intTime;
  unsigned long long int outTime_ntn = (unsigned long long int)(outTime/dt);
  
  
  //for timeseries output
  outputfile.precision(10); //outputfile must be declared in global.hpp
  unsigned int out_k = 0;
  unsigned int out_noutk = (int)katana::getCmdOption(argv, argv+argc, "-nout" , 1);
  double out_dt = dt * (double)out_noutk;
  
  //output timeseries of time and all dynamical variables
  //slows down programm quite a lot due to constant writing to the constant file ->>> bad!
  auto outputAllToFile = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outputfile << tn*dt << "\t" << norm(x->E) << "\t" << x->G << "\t" << x->Q << "\t" << x->E.real() << "\t" << x->E.imag() << std::endl; 
    }
    out_k = (out_k+1)%out_noutk;
  };
  
  
#ifdef SHIFTEDTIME
  auto NetGain = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outputfile << tn*dt << "\t" << norm(x->E) << "\t";
      outputfile << Xhist.at_Rnd_rPtr(p.T01)->G << "\t" << Xhist.at_Rnd_rPtr(p.T01+p.T12)->Q << "\t" ;
      outputfile << 0.5*(Xhist.at_Rnd_rPtr(2.0*p.T01)->G + Xhist.at_Rnd_rPtr(p.T)->G) + Xhist.at_Rnd_rPtr(p.T)->Q + log(p.kappa) << "\t";
//       outputfile << real(x->E) << "\t" << imag(x->E) << "\t";
//       const std::complex<double> EE = Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*( Xhist.at_Rnd_rPtr(p.T)->G + Xhist.at_Rnd_rPtr(2.0*p.T01)->G ) + (1.0-sm::img*p.alpha_Q)*Xhist.at_Rnd_rPtr(p.T)->Q );
//       outputfile << real(EE) << "\t" << imag(EE) << "\t";
      outputfile << std::endl;
    }
    out_k = (out_k+1)%out_noutk;
  };
  
  auto EfieldEvo = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outputfile << tn*dt << "\t" << norm(Xhist.at_Rnd_rPtr(p.T)->E) << "\t";
      outputfile << norm( Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*(Xhist.at_Rnd_rPtr(p.T)->G) ) ) << "\t";
      outputfile << norm( Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*(Xhist.at_Rnd_rPtr(p.T)->G) + (1.0-sm::img*p.alpha_Q)*Xhist.at_Rnd_rPtr(p.T)->Q ) ) << "\t";
      outputfile << norm( Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*( Xhist.at_Rnd_rPtr(p.T)->G + Xhist.at_Rnd_rPtr(2.0*p.T01)->G ) + (1.0-sm::img*p.alpha_Q)*Xhist.at_Rnd_rPtr(p.T)->Q ) ) << "\t";
      
      outputfile << "\t" << norm(x->E) << "\t";
      outputfile << std::endl;
    }
    out_k = (out_k+1)%out_noutk;
  };
#endif
  
  
#ifdef REALTIME
  auto NetGain = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outputfile << tn*dt << "\t" << norm(x->E) << "\t" << x->G << "\t" << x->Q << "\t" << 0.5*(Xhist.at_Rnd_rPtr(p.T01)->G + Xhist.at_Rnd_rPtr(p.T01+2.0*p.T12)->G) + Xhist.at_Rnd_rPtr(p.T01+p.T12)->Q << std::endl;
    }
    out_k = (out_k+1)%out_noutk;
  };
  
  
  auto EfieldEvo = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outputfile << tn*dt << "\t" << norm(Xhist.at_Rnd_rPtr(p.T)->E) << "\t";
      outputfile << norm( Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*(Xhist.at_Rnd_rPtr(p.T01+2.0*p.T12)->G) ) ) << "\t";
      outputfile << norm( Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*(Xhist.at_Rnd_rPtr(p.T01+2.0*p.T12)->G) + (1.0-sm::img*p.alpha_Q)*Xhist.at_Rnd_rPtr(p.T01+p.T12)->Q ) ) << "\t";
      outputfile << norm( Xhist.at_Rnd_rPtr(p.T)->E * sm::ExpC( 0.5*(1.0-sm::img*p.alpha_G)*(Xhist.at_Rnd_rPtr(p.T01)->G+Xhist.at_Rnd_rPtr(p.T01+2.0*p.T12)->G) + (1.0-sm::img*p.alpha_Q)*Xhist.at_Rnd_rPtr(p.T01+p.T12)->Q ) ) << "\t";
      outputfile << norm(x->E) << "\t";
      outputfile << std::endl;
    }
    out_k = (out_k+1)%out_noutk;
  };
#endif
  
  //output to vectors for TS analysis
  std::vector<double> outDataVec;
  std::vector<std::complex<double>> outDataVecComplex;
  outDataVec.reserve((int)(outTime/out_dt));
  outDataVecComplex.reserve((int)(outTime/out_dt));
  auto outputToVector = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
    if(tn >= tn_final - outTime_ntn && out_k == 0){
      outDataVec.push_back(norm(x->E));
      outDataVecComplex.push_back(x->E);
    }
    out_k = (out_k+1)%out_noutk;
  };
  
  
  //////////////////////////////////////////
  //for sweeps
  //////////////////////////////////////////
  
  std::ofstream out_sweep_ML_State, out_sweep_Ext, out_sweep_Ext_all, out_sweep_DeltaT, out_sweep_TS;
  out_sweep_ML_State.precision(10), out_sweep_Ext_all.precision(10), out_sweep_TS.precision(10), out_sweep_Ext.precision(10), out_sweep_DeltaT.precision(10);
  
  evalMLTS outEv_sweep;
  outEv_sweep.PW_rel_searchRadius = T_maxPW;
  outEv_sweep.doubleCountTol = T_dCountTol;
  outEv_sweep.pulse_rel_max_th = T_relMax;
  outEv_sweep.TVPD_rel_pulse_trig = T_TVPD_trig;
  outEv_sweep.TVPD_rel_pulse_bounds = T_TVPD_bound;
  outEv_sweep.TSSmoothing = 1;
  
  double* PtrToSweepPar;
  std::string str_sweep_parameter;
  double* PtrToSecPar;
  std::string str_sweep_sec_parameter;
  
  double sweep_IntTime = katana::getCmdOption(argv, argv+argc, "-sIntTime" , 15000.0);
  double sweep_OutTime = katana::getCmdOption(argv, argv+argc, "-sOutTime" , 5000.0);
  int nsweep_steps = katana::getCmdOption(argv, argv+argc, "-sSteps" , 100);
  double sweep_start = katana::getCmdOption(argv, argv+argc, "-sStart" , 0);
  double sweep_end = katana::getCmdOption(argv, argv+argc, "-sEnd" , 1);
  
  int sweep_nAllExt_out = 2000;
  int sweep_nDeltaT_out = 100;
  
  std::string str_MLL_state;
  std::string str_extrema;
  std::string str_extrema_all;
  std::string str_DeltaT;
  std::string str_powerSpec;
  std::string str_opticalSpec;
  std::string str_AC;
  std::string str_IntAC;
  std::string str_TS;
  std::string str_bin;
  std::string str_netgain;
  
  std::string str_sweep_data = "data/";
  std::string str_sweep_updown = "up";
  
  bool bool_powerSpec = katana::getCmdOption_bool(argv, argv+argc, "-wpowerSpec" , false);
  bool bool_opticalSpec = katana::getCmdOption_bool(argv, argv+argc, "-wopticalSpec" , false);
  bool bool_AC = katana::getCmdOption_bool(argv, argv+argc, "-wAC" , false);
  bool bool_IntAC = katana::getCmdOption_bool(argv, argv+argc, "-wIntAC" , false);
  bool bool_TS = katana::getCmdOption_bool(argv, argv+argc, "-wTS" , false);
  bool bool_bin = katana::getCmdOption_bool(argv, argv+argc, "-wBin" , false);
  bool bool_Extrema = katana::getCmdOption_bool(argv, argv+argc, "-wExtrema" , false);
  bool bool_Extrema_all = katana::getCmdOption_bool(argv, argv+argc, "-wExtrema_all" , false);
  bool bool_DeltaT = katana::getCmdOption_bool(argv, argv+argc, "-wDeltaT" , false);
  bool bool_NoSweep = katana::getCmdOption_bool(argv, argv+argc, "-wNoSweep" , false);
  bool bool_LogSweep = katana::getCmdOption_bool(argv, argv+argc, "-wLogSweep" , false);
  bool bool_netgain = katana::getCmdOption_bool(argv, argv+argc, "-wNetgain" , false);
  
  auto sweepMLL = [&](){
    
    
    std::vector<double> SweepVals(nsweep_steps+1);
    SweepVals[0] = sweep_start;

    if(bool_LogSweep){
      double SweepLogMult = pow( (sweep_end/sweep_start), 1.0/((double)(nsweep_steps)) );
      for(unsigned int k = 1; k < SweepVals.size(); k++){
        SweepVals[k] = SweepVals[k-1]*SweepLogMult;
//         cout << k << " " <<  SweepVals[k] << endl;
      }
    }
    else{
      double SweepLinAdd = (sweep_end-sweep_start)/(double)nsweep_steps;
      for(unsigned int k = 1; k < SweepVals.size(); k++){
        SweepVals[k] = SweepVals[k-1]+SweepLinAdd;
//         cout << k << " " <<  SweepVals[k] << endl;
      }
    }
    
    
    std::ostringstream ostringstream_DoubleToStr;
    ostringstream_DoubleToStr << std::fixed << std::setprecision(8);
//     ostringstream_DoubleToStr << std::scientific << std::setprecision(8);
    ostringstream_DoubleToStr << *PtrToSecPar;
//     std::string StrSecParVal = "";
    std::string StrSecParVal = ostringstream_DoubleToStr.str();
    ostringstream_DoubleToStr.str("");
    ostringstream_DoubleToStr.clear();

    
    str_MLL_state = str_sweep_data + str_sweep_updown + "_sweep_MLL_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_extrema = str_sweep_data + str_sweep_updown + "_sweep_extrema_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_powerSpec = str_sweep_data + "powerSpec/powerSpec_"+str_sweep_updown+"_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_opticalSpec = str_sweep_data + "opticalSpec/opticalSpec_"+str_sweep_updown+"_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_AC = str_sweep_data + "AC/AC_" + str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_IntAC = str_sweep_data + "IntAC/IntAC_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_TS = str_sweep_data + "TS/TS_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_bin = str_sweep_data + "bin/bin_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_extrema_all = str_sweep_data + str_sweep_updown + "_sweep_all_extrema_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_DeltaT = str_sweep_data + str_sweep_updown + "_sweep_DeltaT_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_netgain = str_sweep_data + str_sweep_updown + "_sweep_netgain_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    str_netgain = str_sweep_data + "netgain/netgain_"+str_sweep_updown + "_"+str_sweep_sec_parameter+"_" + StrSecParVal + str_suffix;
    
    outTime = sweep_OutTime;
    outTime_ntn = (unsigned long long int)(outTime/dt);
    
    outDataVec.reserve(outTime_ntn);
    outDataVecComplex.reserve(outTime_ntn);
    
    out_sweep_ML_State.open(str_MLL_state);
    if(bool_Extrema) out_sweep_Ext.open(str_extrema);
    if(bool_Extrema_all) out_sweep_Ext_all.open(str_extrema_all);
    if(bool_DeltaT) out_sweep_DeltaT.open(str_DeltaT);

    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    
    for(unsigned int k = 0; k < SweepVals.size(); k++){
      *PtrToSweepPar = SweepVals[k];  
      
      ostringstream_DoubleToStr << std::fixed << std::setprecision(8);
      ostringstream_DoubleToStr << *PtrToSweepPar;
      std::string StrSweepParVal = ostringstream_DoubleToStr.str();
      ostringstream_DoubleToStr.str("");
      ostringstream_DoubleToStr.clear();
      
        
      sm::clck = clock();
      if(bool_NoSweep) {
        //read history file
        if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
          std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
          Xhist.load(histFile);
          tn = Xhist.tn;
        }
        else{
          tn=0;
          Xhist.setToCnst(1E-6);
        }
      }
        
      tn=0;
      outDataVec.resize(0), outDataVecComplex.resize(0);
      computeParas();
      DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, sweep_IntTime, dt, outputToVector);
      
      outEv_sweep.loadNewTS(&outDataVec, out_dt, p.T);
      
      
      cout << str_sweep_parameter << ": " << *PtrToSweepPar;
      out_sweep_ML_State << *PtrToSweepPar << "\t" << *PtrToSecPar << "\t"; // 2
      //basics
      out_sweep_ML_State << outEv_sweep.TSAverage << "\t" << outEv_sweep.TSSmallestMin << "\t" << outEv_sweep.TSGreatestMax << "\t"; // 5
      //from max detection:
      out_sweep_ML_State << outEv_sweep.MLnMax() << "\t" << outEv_sweep.meanMax <<  "\t" << outEv_sweep.meanMaxTSep << "\t"; // 8
      out_sweep_ML_State << outEv_sweep.PW_mean_PFWHM << "\t" << outEv_sweep.PW_mean_PSD << "\t"; // 10
      out_sweep_ML_State << outEv_sweep.PtPJ_MaxPos() << "\t" << outEv_sweep.RelAmpJitter() << "\t"; // 12
      //from IntAC:
      out_sweep_ML_State << outEv_sweep.ML_Fund_RR_from_IntAC() << "\t" << outEv_sweep.ML_FWHM_IntAC() << "\t" << outEv_sweep.ML_StdDev_IntAC() << "\t"; //15
      //from TriggerValPulseDetection_Statistics
      out_sweep_ML_State << outEv_sweep.TVPD_meanPulseMax << "\t" << outEv_sweep.TVPD_meanPulseArea << "\t" << outEv_sweep.TVPD_meanPulseSD << "\t" << outEv_sweep.TVPD_meanPulseSep << "\t"; // 19
      out_sweep_ML_State << sqrt(outEv_sweep.Variance(outEv_sweep.TVPD_PM_vec))/outEv_sweep.TVPD_meanPulseMax << "\t" << sqrt(outEv_sweep.Variance(outEv_sweep.TVPD_PA_vec))/outEv_sweep.TVPD_meanPulseArea << "\t" << outEv_sweep.PtPJitter(outEv_sweep.TVPD_PS_vec) << "\t"; //22
      // changing triggers and thresholds for evaluation
      double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv_sweep.TriggerValPulseDetection_Statistics_forEval(0.05, 0.002);
      out_sweep_ML_State << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //28
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv_sweep.TriggerValPulseDetection_Statistics_forEval(0.05, 0.01);
      out_sweep_ML_State << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //34
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv_sweep.TriggerValPulseDetection_Statistics_forEval(0.10, 0.02);
      out_sweep_ML_State << ev_TSep << "\t" << ev_TSsd << "\t" << ev_PM << "\t" << ev_PMsd << "\t" << ev_PA << "\t" << ev_PAsd << "\t"; //40
      //IntAC Eval
      double ev_fMa, ev_fMap, ev_fMi, ev_fMip;
      tie(ev_fMa, ev_fMap, ev_fMi, ev_fMip) = outEv_sweep.EvalIntAC();
      out_sweep_ML_State << ev_fMa << "\t" << ev_fMap << "\t" << ev_fMi << "\t" << ev_fMip << "\t"; //44
      //IntAC max pos
      std::vector<double> IACmp;
      const unsigned int IACnmoout = 5;
      outEv_sweep.IntAC_MaxAboveTh(IACmp, 0.2);
      for(unsigned int k = 1; k < IACnmoout+1; k++){
        if(k < IACmp.size()){ out_sweep_ML_State << IACmp[k] << "\t"; } //49
        else{ out_sweep_ML_State << -1.0 << "\t"; }
      }
      outEv_sweep.IntAC_MaxAboveTh(IACmp, 0.1);
      for(unsigned int k = 1; k < IACnmoout+1; k++){
        if(k < IACmp.size()){ out_sweep_ML_State << IACmp[k] << "\t"; } //54
        else{ out_sweep_ML_State << -1.0 << "\t"; }
      }
      outEv_sweep.IntAC_MaxAboveTh(IACmp, 0.01);
      for(unsigned int k = 1; k < IACnmoout+1; k++){
        if(k < IACmp.size()){ out_sweep_ML_State << IACmp[k] << "\t"; } //59
        else{ out_sweep_ML_State << -1.0 << "\t"; }
      }
      //new line
      out_sweep_ML_State << std::endl;
      
      
      if(bool_Extrema){
        out_sweep_Ext << *PtrToSweepPar << "\t" << outEv_sweep.TSGreatestMax << std::endl;
        for(std::size_t k = 0; k < outEv_sweep.uniqueMax.size(); k++) out_sweep_Ext << *PtrToSweepPar << "\t" << outEv_sweep.uniqueMax[k] << std::endl;      
        out_sweep_Ext << *PtrToSweepPar << "\t" <<  outEv_sweep.TSSmallestMin << std::endl;
        for(std::size_t k = 0; k < outEv_sweep.uniqueMin.size(); k++) out_sweep_Ext << *PtrToSweepPar << "\t" << outEv_sweep.uniqueMin[k] << std::endl;
      }
      
      if(bool_Extrema_all){
        out_sweep_Ext_all << *PtrToSweepPar << "\t" << outEv_sweep.TSGreatestMax << std::endl;
        for(int k = 0; k < std::min((int)outEv_sweep.maxima.size(),sweep_nAllExt_out); k++) out_sweep_Ext_all << *PtrToSweepPar << "\t" << outEv_sweep.maxima[k] << std::endl;      
        out_sweep_Ext_all << *PtrToSweepPar << "\t" <<  outEv_sweep.TSSmallestMin << std::endl;
        for(int k = 0; k < std::min((int)outEv_sweep.uniqueMin.size(),sweep_nAllExt_out); k++) out_sweep_Ext_all << *PtrToSweepPar << "\t" << outEv_sweep.uniqueMin[k] << std::endl;
      }
    
      if(bool_DeltaT){
        for(int k = 0; k < std::min((int)((int)outEv_sweep.maxima_tvec.size()-1),sweep_nDeltaT_out); k++) out_sweep_DeltaT << *PtrToSweepPar << "\t" << (outEv_sweep.maxima_tvec[k+1] - outEv_sweep.maxima_tvec[k]) << std::endl;
      }
      
      if(bool_TS){
        out_sweep_TS.open(str_TS+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
        for(std::size_t k = 0; k < outDataVec.size()/10; k++) out_sweep_TS << k*out_dt << "\t" << outDataVec[k] << std::endl;
        out_sweep_TS.close();
      }
    
      if(bool_IntAC){
        katana::sm_print_correlation(outEv_sweep.IntAC, out_dt, str_IntAC+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix, 1.1);
      }
      
      if(bool_AC){
        std::vector<double> AC;
        katana::get_autocorr_norm_out(outDataVecComplex, AC);
        katana::sm_print_correlation(AC, out_dt, str_AC+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix, 10.0);
      }
      
      if(bool_powerSpec){
        std::vector<double> powerSpec;
        std::vector<double> WindowedData;
        sm::window_Hann(outDataVec, WindowedData);
        sm::get_power_spec(WindowedData, powerSpec, out_dt);
        sm::dump_power_spec(powerSpec, out_dt, 5.0, str_powerSpec+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
      }
      
      if(bool_opticalSpec){
        std::vector<double> opticalSpec;
        sm::get_optical_spec(outDataVecComplex, opticalSpec);
        sm::dump_optical_spec(opticalSpec, out_dt, 5000.0, str_opticalSpec+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
      }
      
      if(bool_bin){
        Xhist.tn = tn;
        Xhist.save(str_bin+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
      }
      
      
      if(bool_netgain){
        outputfile.open(str_netgain+"_"+str_sweep_parameter+"_"+StrSweepParVal + str_suffix);
        DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, 4.0, dt, NetGain);
        outputfile.close();
      }

      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_p = std::chrono::system_clock::to_time_t(now);
      cout << ",\tintegration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << ",\tcurrent time and date: " << std::put_time(std::localtime(&now_p), "%F %T") << endl;
      
    }
    out_sweep_Ext.close();
    out_sweep_ML_State.close();
    if(bool_Extrema_all) out_sweep_Ext_all.close();
    if(bool_DeltaT) out_sweep_DeltaT.close();
    
  };
  
  
  //////////////////////////////////////////
  //simulations
  //////////////////////////////////////////
  
  
  //timeseries to file
  if(katana::getCmdOption_bool(argv, argv+argc, "-simpleTS" , false)){
    outputfile.open("data/out_simpleTS");
    tn=0;
    
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string binhist = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(binhist);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    
    sm::clck = clock();
//     DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, EfieldEvo);
    DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, NetGain);
//     DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, outputAllToFile);
//     DDEsolver.DDE_euler(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, outputAllToFile);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    outputfile.close();
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      Xhist.save("data/Xhist.bin");
    }
  }
  
  //compute and analize time series
  if(katana::getCmdOption_bool(argv, argv+argc, "-TS" , false)){
    evalMLTS outEv;
    outEv.PW_rel_searchRadius = T_maxPW;
    outEv.doubleCountTol = T_dCountTol;
    outEv.pulse_rel_max_th = T_relMax;
    outEv.TVPD_rel_pulse_trig = T_TVPD_trig;
    outEv.TVPD_rel_pulse_bounds = T_TVPD_bound;
    outEv.TSSmoothing = 1;
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      tn=0;
      Xhist.setToCnst(1E-6);
    }
    
    outDataVec.reserve(outTime_ntn);
    outDataVecComplex.reserve(outTime_ntn);

    sm::clck = clock();
    DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
//     DDEsolver.DDE_euler(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
    cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
    
    outEv.loadNewTS(&outDataVec, out_dt, p.T);
    

    if(!katana::getCmdOption_bool(argv, argv+argc, "-noTSout" , false)){
//       double TSoutTime = katana::getCmdOption(argv, argv+argc, "-TSoutTime" , outTime);
      unsigned int TSoutn = (int)(katana::getCmdOption(argv, argv+argc, "-TSoutTime" , outTime)/out_dt);
      outputfile.open("data/out_TS");
      for(std::size_t k = outDataVec.size() - TSoutn; k < outDataVec.size(); k++){
        outputfile << k*out_dt << "\t" << outDataVec[k] << "\t";
  //       outputfile << real(outDataVecComplex[k]) << "\t" << imag(outDataVecComplex[k]) << "\t" << arg(outDataVecComplex[k]) << "\t";
        outputfile << endl;
      }
      outputfile.close();
    }

    
    
//     std::vector<double> Phase,Freq;
//     double avOptFreq = outEv.opticalFreq(outDataVecComplex, Phase, Freq, out_dt);
//     cout << "avOptFreq: " << avOptFreq << endl;
//     
//     outputfile.open("data/out_TS_phase");
//     for(std::size_t k = 0; k < Freq.size(); k++) outputfile << k*out_dt << "\t" << Phase[k] << "\t" << Freq[k] << endl;
//     outputfile.close();
   
//     
//     cout << "----input----" << endl;
    cout << "----results----" << endl;
    cout << "average output power: " << outEv.TSAverage << " greatest peak power: " << outEv.TSGreatestMax << " TSSmallestMin: " << outEv.TSSmallestMin << endl;
    cout << "----from maxima detection:" << endl;
    cout << "number of unique ML max: " << outEv.MLnMax() << " mean peak power: " << outEv.meanMax <<  " mean max t seperation: " << outEv.meanMaxTSep << endl;
    cout << "PtP jitter from max pos: " << outEv.PtPJ_MaxPos() << " amplitude jitter from max pos: " << outEv.RelAmpJitter() << endl;
    cout << "PW from pulse FWHM: " << outEv.PW_mean_PFWHM << " PW from pulse StdDev: " << outEv.PW_mean_PSD << endl;
    cout << "----from the intensity AC:" << endl;
    cout << "ML fundamental RR from IntAC: " << outEv.ML_Fund_RR_from_IntAC() << endl;
    cout << "IntAC FWHM: " << outEv.ML_FWHM_IntAC() << " IntAC StdDev: " << outEv.ML_StdDev_IntAC() << endl;
    cout << "----from weighted pulse positions:" << endl;
    cout << "mean TVPD peak power: " << outEv.TVPD_meanPulseMax << " mean TVPD pulse area: " << outEv.TVPD_meanPulseArea << " TVPD PW from PStdDev: " << outEv.TVPD_meanPulseSD << " TVPD mean pulse sep: " << outEv.TVPD_meanPulseSep << endl;
    cout << "TVPD peak power rel jitter: " << sqrt(outEv.Variance(outEv.TVPD_PM_vec))/outEv.TVPD_meanPulseMax << " TVPD pulse area rel jitter: " << sqrt(outEv.Variance(outEv.TVPD_PA_vec))/outEv.TVPD_meanPulseArea << endl;
    cout << "TVPD pulse position jitter: " << outEv.PtPJitter(outEv.TVPD_PS_vec) << endl;
    double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
    tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.005, 0.001);
    cout << "ev_TSep: " <<  ev_TSep << " " << ev_TSsd << " ev_TSsd: " << ev_PM << " ev_PM: " << ev_PMsd << " ev_PMsd: " << ev_PA << " ev_PA: " << ev_PAsd << endl;
    
    
//     outEv.dump_IACmaxpos(0.2);
//     outEv.dump_IACmaxpos(0.1);
//     outEv.dump_IACmaxpos(0.01);
    
//     std::vector<double> IACmp;
//     outEv.IntAC_MaxAboveTh(IACmp, 0.2);
//     cout << "IACmp size: " << IACmp.size() << endl;
//     for(unsigned int k = 0; k < IACmp.size(); k++){
//       cout << k << " " << IACmp[k] << endl;
//     }
    
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-powerSpec" , false)){
      std::vector<double> powerSpec;
      std::vector<double> WindowedData;
      
      sm::get_power_spec(outDataVec, powerSpec, out_dt);
      sm::dump_power_spec(powerSpec, out_dt, 10, "data/out_TS_powerSpec");
      
      sm::window_Hann(outDataVec, WindowedData);
      sm::get_power_spec(WindowedData, powerSpec, out_dt);
      sm::dump_power_spec(powerSpec, out_dt, 10, "data/out_TS_powerSpec_Hann");
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-opticalSpec" , false)){
      std::vector<double> opticalSpec;
      std::vector<std::complex<double>> WindowedData;
      
      sm::get_optical_spec(outDataVecComplex, opticalSpec);
      sm::dump_optical_spec(opticalSpec, out_dt, 5000, "data/out_TS_opticalSpec"+str_suffix);
      
      sm::window_Hann(outDataVecComplex, WindowedData);
      sm::get_optical_spec(WindowedData, opticalSpec);
      sm::dump_optical_spec(opticalSpec, out_dt, 5000, "data/out_TS_opticalSpec_Hann"+str_suffix);
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-AC" , false)){
      std::vector<double> AC;
      katana::get_autocorr_norm_out(outDataVecComplex, AC);
      katana::sm_print_correlation(AC, out_dt, "data/out_TS_AC", 1.5);
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-IntAC" , false)){
      katana::sm_print_correlation(outEv.IntAC, out_dt, "data/out_TS_IntAC", 1.5);
    }
    
    //save history to binary file
    if(katana::getCmdOption_bool(argv, argv+argc, "-saveHist" , false)){
      Xhist.tn = tn;
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-shistFile" , "data/Xhist.bin");
      Xhist.save(histFile);
    }
    
    if(katana::getCmdOption_bool(argv, argv+argc, "-netGain" , false)){
      outputfile.open("data/out_netGain"+str_suffix);
      DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, 5.0, dt, NetGain);
      outputfile.close();
    }
    
  }
  

  //compute LTTJ
  if(katana::getCmdOption_bool(argv, argv+argc, "-LTTJ" , false)){
    noise = GWnoise;
    
    evalMLTS outEv;
    outEv.PW_rel_searchRadius = T_maxPW;
    outEv.doubleCountTol = T_dCountTol;
    outEv.pulse_rel_max_th = T_relMax;
    outEv.TVPD_rel_pulse_trig = T_TVPD_trig;
    outEv.TVPD_rel_pulse_bounds = T_TVPD_bound;
    outEv.TSSmoothing = 1;
    
    
    bool LTTJpowerSpec = katana::getCmdOption_bool(argv, argv+argc, "-LTTJ_wPS" , false);
    double LTTJpowerSpecOut = katana::getCmdOption(argv, argv+argc, "-LTTJ_PSOut" , 51);
    
    //read history file
    if(katana::getCmdOption_bool(argv, argv+argc, "-loadHist" , false)){
      std::string histFile = katana::getCmdOption(argv, argv+argc, "-histFile" , "data/Xhist.bin");
      Xhist.load(histFile);
      tn = Xhist.tn;
    }
    else{
      cout << "calculate new TS..." << endl;
        
      tn=0;
      Xhist.setToCnst(1E-6);
      
      outDataVec.reserve(outTime_ntn);
      outDataVecComplex.reserve(outTime_ntn);
      
      outDataVec.resize(0), outDataVecComplex.resize(0);
      computeParas();
      
      sm::clck = clock();
      DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, intTime, dt, outputToVector);
      cout << "integration time: " << clock() - sm::clck << "(" << (float)(clock() - sm::clck)/CLOCKS_PER_SEC << " seconds)" << endl;
      
      outEv.loadNewTS(&outDataVec, out_dt, p.T);
        
      outputfile.open("data/LTTJ_TS");
      for(std::size_t k = 0; k < outDataVec.size(); k++){
        outputfile << k*out_dt << "\t" << outDataVec[k] << "\t";
        outputfile << endl;
      }
      outputfile.close();
      
      cout << "----results----" << endl;
      cout << "average output power: " << outEv.TSAverage << " greatest peak power: " << outEv.TSGreatestMax << " TSSmallestMin: " << outEv.TSSmallestMin << endl;
      cout << "----from maxima detection:" << endl;
      cout << "number of unique ML max: " << outEv.MLnMax() << " mean peak power: " << outEv.meanMax <<  " mean max t seperation: " << outEv.meanMaxTSep << endl;
      cout << "PtP jitter from max pos: " << outEv.PtPJ_MaxPos() << " amplitude jitter from max pos: " << outEv.RelAmpJitter() << endl;
      cout << "PW from pulse FWHM: " << outEv.PW_mean_PFWHM << " PW from pulse StdDev: " << outEv.PW_mean_PSD << endl;
      cout << "----from the intensity AC:" << endl;
      cout << "ML fundamental RR from IntAC: " << outEv.ML_Fund_RR_from_IntAC() << endl;
      cout << "IntAC FWHM: " << outEv.ML_FWHM_IntAC() << " IntAC StdDev: " << outEv.ML_StdDev_IntAC() << endl;
      cout << "----from weighted pulse positions:" << endl;
      cout << "mean TVPD peak power: " << outEv.TVPD_meanPulseMax << " mean TVPD pulse area: " << outEv.TVPD_meanPulseArea << " TVPD PW from PStdDev: " << outEv.TVPD_meanPulseSD << " TVPD mean pulse sep: " << outEv.TVPD_meanPulseSep << endl;
      cout << "TVPD peak power rel jitter: " << sqrt(outEv.Variance(outEv.TVPD_PM_vec))/outEv.TVPD_meanPulseMax << " TVPD pulse area rel jitter: " << sqrt(outEv.Variance(outEv.TVPD_PA_vec))/outEv.TVPD_meanPulseArea << endl;
      cout << "TVPD pulse position jitter: " << outEv.PtPJitter(outEv.TVPD_PS_vec) << endl;
      double ev_TSep, ev_TSsd, ev_PM, ev_PMsd, ev_PA, ev_PAsd;
      tie(ev_TSep,ev_TSsd,ev_PM,ev_PMsd,ev_PA,ev_PAsd) = outEv.TriggerValPulseDetection_Statistics_forEval(0.005, 0.001);
      cout << ev_TSep << " " << ev_TSsd << " " << ev_PM << " " << ev_PMsd << " " << ev_PA << " " << ev_PAsd << endl;
      
    }

    
    if(katana::getCmdOption_bool(argv, argv+argc, "-LTTJ_buffer" , false)){
      double LTTJ_bufferIntTime = katana::getCmdOption(argv, argv+argc, "-LTTJ_BIntTime" , 500.0*p.T);
      cout << "integrate buffer time " << LTTJ_bufferIntTime << endl;
      DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, LTTJ_bufferIntTime, dt, outputfunc_empty);
    }
    
    
    //get initial pulse statistics
    cout << "test mode-locking state..." << endl;
    double testTime = katana::getCmdOption(argv, argv+argc, "-LTTJ_testTime" , 100.0*p.T);
    outTime_ntn = (unsigned long long int)(testTime/dt);
    
    outDataVec.reserve(outTime_ntn);
    outDataVecComplex.reserve(outTime_ntn);
    
    outDataVec.resize(0), outDataVecComplex.resize(0);
    computeParas();
    
    DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, testTime, dt, outputToVector);
    
    outEv.loadNewTS(&outDataVec, out_dt, p.T);
    
    if(outEv.TVPD_meanPulseSep < 0.9*p.T || outEv.TVPD_meanPulseSep > 1.1*p.T || outEv.PtPJitter(outEv.TVPD_PS_vec) > 0.1*p.T ){
      cout << "no evidence for FML found" << endl;
      return 0;
    }
    else{
      cout << "evidence for FML found" << endl;
    }

    
      
    cout << "starting long term timing jitter calculation!" << endl;
    
    tn = 0;
    double jitterCalcT = katana::getCmdOption(argv, argv+argc, "-LTTJ_intTime" , 1000.0*p.T);
    outTime = jitterCalcT;
    outTime_ntn = (unsigned long long int)(outTime/dt);
    unsigned int nNoiseRel = katana::getCmdOption(argv, argv+argc, "-LTTJ_nNoiseRel" , 50);
    std::vector<double> T_CN (nNoiseRel);
    std::vector<std::vector<double>> PulseTimes (nNoiseRel);
    
    std::string str_out_LTTJ_suffix = katana::getCmdOption(argv, argv+argc, "-LTTJ_str_suf" , "");
    
//     //output only intensity to vector for analysis
//     std::vector<double> outDataVecLTTJ;
//     outDataVecLTTJ.reserve(outTime_ntn+5);
//     auto outputToVectorLTTJ = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
//       if(tn >= tn_final - outTime_ntn){
//         outDataVecLTTJ.push_back(norm(x->E));
//       }
//     };
//     
    
    const double nsigmas = 10;
    double SampleDt = katana::getCmdOption(argv, argv+argc, "-SampleDt" , 10*dt);
    unsigned int SampleDt_dtsteps = (int)(SampleDt/dt);
    unsigned int ConvKernelSize = (int)(nsigmas*SampleDt/dt);
    
    std::vector<double> ConvKernel(ConvKernelSize);
    std::vector<double> LTTJPS_ShortHist(ConvKernelSize);
    unsigned int LTTJPS_SHind = 0;
    std::vector<double> ConvDataVec;
    
    if(LTTJpowerSpec){
      jitterCalcT += nsigmas*SampleDt;
      outTime_ntn += ConvKernelSize;
      
      cout << "SampleDt: " << SampleDt << ", SampleDt / dt: " << SampleDt_dtsteps << ", max frequency 1000/(2*SampleDt): " << 1.0/(2.0*SampleDt) << endl;
      
      double ConvNorm = 0.0;
      for(unsigned int k = 0; k<ConvKernelSize; k++){
          double v = exp(-0.5*(k-ConvKernelSize/2.0)*(k-ConvKernelSize/2.0)/((double)(0.5*SampleDt_dtsteps*0.5*SampleDt_dtsteps)));
          ConvKernel[k] = v;
          ConvNorm += v;
      }
      for(unsigned int k = 0; k<ConvKernelSize; k++){
        ConvKernel[k] = ConvKernel[k]/ConvNorm;
      }
      
      ConvDataVec.reserve((int)(outTime/SampleDt));
      cout << "reseverd ConvDataVec size = 8byte * " << (int)(outTime/SampleDt) << " = " << 8e-6*outTime/SampleDt << " Mbyte" << endl;
    }
    
    
    double LTTJ_TSGreatestMax = outEv.TSGreatestMax;
    double LTTJ_PulseTrigTh = 0.1;
    double LTTJ_PulseBound = 1e-3;
    bool LTTJ_pulse = false;
    
    double LTTJ_PulseNorm = 0.0;
    double LTTJ_FirstMoment = 0.0;
    
    unsigned int LTTJ_ShortHist_Size = (int)(0.5*p.T/dt);
    std::vector<double> LTTJ_ShortHist(LTTJ_ShortHist_Size);
    std::vector<double> PulsePositions;
    PulsePositions.resize(0);
    
    auto LTTJ_MaxPosOut = [&](vars* x, unsigned long long int tn, unsigned long long int tn_final){
      if(tn >= tn_final - outTime_ntn){
        
        const double LTTJ_TSVal = norm(x->E);
        LTTJ_ShortHist[tn%LTTJ_ShortHist_Size] = LTTJ_TSVal;
        
        if(LTTJ_pulse == false 
          && LTTJ_TSVal > LTTJ_TSGreatestMax*LTTJ_PulseTrigTh 
          && (PulsePositions.size() == 0 || (tn*dt) > PulsePositions.back() + 0.5*p.T )
        ){
          LTTJ_pulse = true;
          unsigned int k = 1;
          unsigned int SHind = ((tn-k)%LTTJ_ShortHist_Size);
          while(LTTJ_ShortHist[SHind] > LTTJ_TSGreatestMax*LTTJ_PulseBound){
            LTTJ_FirstMoment += LTTJ_ShortHist[SHind] * (tn-k) * dt;
            LTTJ_PulseNorm += LTTJ_ShortHist[SHind];
            k++;
            SHind = ((tn-k)%LTTJ_ShortHist_Size);
          }
        }
        if(LTTJ_pulse == true && LTTJ_TSVal > LTTJ_TSGreatestMax*LTTJ_PulseBound){
          LTTJ_FirstMoment += LTTJ_TSVal * tn * dt;
          LTTJ_PulseNorm += LTTJ_TSVal;
        }
        if(LTTJ_pulse == true
          && LTTJ_TSVal < LTTJ_TSGreatestMax*LTTJ_PulseBound 
          && LTTJ_ShortHist[((tn-1)%LTTJ_ShortHist_Size)] < LTTJ_TSGreatestMax*LTTJ_PulseBound 
          && LTTJ_ShortHist[((tn-2)%LTTJ_ShortHist_Size)] < LTTJ_TSGreatestMax*LTTJ_PulseBound){
          LTTJ_pulse = false;
//           cout << LTTJ_FirstMoment / LTTJ_PulseNorm << endl;
          PulsePositions.push_back(LTTJ_FirstMoment / LTTJ_PulseNorm);
          LTTJ_FirstMoment = 0.0;
          LTTJ_PulseNorm = 0.0;
        }
      }
      
      if(LTTJpowerSpec && tn >= tn_final - outTime_ntn){
        LTTJPS_ShortHist[LTTJPS_SHind] = norm(x->E);
        LTTJPS_SHind = (LTTJPS_SHind+1)%ConvKernelSize;
        
        if(tn >= tn_final - outTime_ntn + ConvKernelSize){
          if( (LTTJPS_SHind%SampleDt_dtsteps) == 0 ){
            double ConvVal = 0.0;
            for(unsigned int k = 0; k<ConvKernelSize; k++){
              ConvVal += LTTJPS_ShortHist[(k+LTTJPS_SHind)%ConvKernelSize]*ConvKernel[k];
            }
            ConvDataVec.push_back(ConvVal);
          }
        }
      }
    };
    
    
    for(std::size_t k = 0; k < nNoiseRel; k++){
      
      PulsePositions.resize(0);
      ConvDataVec.resize(0);
      LTTJPS_SHind = 0;
      
      cout << "run: " << k+1 << " out of " << nNoiseRel << " noise realizations; T_C["<<k<<"]: ";

      DDEsolver.DDE_RK4(MLL_derivs, noise, &Xhist, &p, &tn, jitterCalcT, dt, LTTJ_MaxPosOut);
      
      PulseTimes[k] = PulsePositions;

      T_CN[k] = (PulsePositions[PulsePositions.size()-1]-PulsePositions[0])/(PulsePositions.size()-1);
      
      cout << T_CN[k];
      std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
      std::time_t now_p = std::chrono::system_clock::to_time_t(now);
      cout << ",\tcurrent time and date: " << std::put_time(std::localtime(&now_p), "%F %T") << endl;
      
      if(LTTJpowerSpec){
      
        cout << "ConvDataVec capacity after computation: " << 8e-6*ConvDataVec.capacity() << " Mbyte" << endl;
        cout << "ConvKernel capacity after computation: " << 8e-6*ConvKernel.capacity() << " Mbyte" << endl;
        cout << "LTTJPS_ShortHist capacity after computation: " << 8e-6*LTTJPS_ShortHist.capacity() << " Mbyte" << endl;
        
        if(katana::getCmdOption_bool(argv, argv+argc, "-LTJJ_PowerSpec_Hann" , false)){
          sm::window_Hann(ConvDataVec);
        }

        sm::get_power_spec(ConvDataVec, SampleDt);
        sm::dump_power_spec(ConvDataVec, SampleDt*1E-3, LTTJpowerSpecOut, "data/LTTJ/powerspecs/powerspec"+str_out_LTTJ_suffix + "_noiseRel_" + to_string(k));
      }
    
    }
    
    //calc and print clock time from all noise realizations
    double T_C = 0;
    for (std::size_t k = 0; k < nNoiseRel; k++) T_C += T_CN[k]/nNoiseRel;
    cout << "T_C: " << T_C << endl;
    
    //determine minimal number of pulses in the pulseTimes vectors (all noise realizations)
    std::size_t PulseTimes_maxsize = PulseTimes[0].size();
    for (std::size_t k = 1; k < nNoiseRel; k++) if(PulseTimes[k].size() < PulseTimes_maxsize) PulseTimes_maxsize = PulseTimes[k].size();
    
    
    //write pulse position matrix
    if(katana::getCmdOption_bool(argv, argv+argc, "-LTTJ_wDeltaT" , false)){
      std::string str_out_LTTJ_PulsePositions = katana::getCmdOption(argv, argv+argc, "-str_out_LTTJ_PulsePositions" , "data/LTTJ/LTTJ_PulsePositions");
      std::ofstream LTTJ_PulsePositions (str_out_LTTJ_PulsePositions + str_out_LTTJ_suffix, std::ios::out | std::ios::trunc);
      LTTJ_PulsePositions.precision(15);
      
      for (std::size_t m = 0; m < PulseTimes_maxsize; m++){
        LTTJ_PulsePositions << m << "\t";
        for (std::size_t k = 0; k < nNoiseRel; k++) LTTJ_PulsePositions << PulseTimes[k][m] << "\t";
        LTTJ_PulsePositions << endl;
      }
    }
    
    //write pulse interval matrix
    std::string str_out_LTTJ_PulseIntervals = katana::getCmdOption(argv, argv+argc, "-str_out_LTTJ_PulseIntervals" , "data/LTTJ/LTTJ_PulseIntervals");
    std::ofstream LTTJ_PulseIntervals (str_out_LTTJ_PulseIntervals + str_out_LTTJ_suffix, std::ios::out | std::ios::trunc);
    LTTJ_PulseIntervals.precision(12);
    for (std::size_t m = 1; m < PulseTimes_maxsize; m++){
      LTTJ_PulseIntervals << m << "\t";
      for (std::size_t k = 0; k < nNoiseRel; k++) LTTJ_PulseIntervals << PulseTimes[k][m] - PulseTimes[k][m-1]<< "\t";
      LTTJ_PulseIntervals << endl;
    }
    
    //create delta_t matrix
    std::vector<std::vector<double>> delta_t (nNoiseRel , std::vector<double>(PulseTimes_maxsize,0.0));
    for (std::size_t m = 1; m < PulseTimes_maxsize; m++){
      for (std::size_t k = 0; k < nNoiseRel; k++) delta_t[k][m] = PulseTimes[k][m] - (PulseTimes[k][0] + m*T_C);
    }
    
    //vectors for delta t statistics
    std::vector<double> mean_delta_t (PulseTimes_maxsize, 0.0);
    std::vector<double> var_delta_t (PulseTimes_maxsize, 0.0);
    
    //output delta_t matrix
    if(katana::getCmdOption_bool(argv, argv+argc, "-LTTJ_wPulsePositions" , false)){
      std::string str_out_LTTJ_deltaT = katana::getCmdOption(argv, argv+argc, "-str_out_LTTJ_deltaT" , "data/LTTJ/LTTJ_deltaT");
      std::ofstream LTTJ_deltaT (str_out_LTTJ_deltaT + str_out_LTTJ_suffix, std::ios::out | std::ios::trunc);
      LTTJ_deltaT.precision(12);
      for (std::size_t m = 1; m < PulseTimes_maxsize; m++){
        LTTJ_deltaT << m << "\t";
        for (std::size_t k = 0; k < nNoiseRel; k++){
          LTTJ_deltaT << delta_t[k][m] << "\t";
        }
        LTTJ_deltaT << endl;
      }
    }
    
    //calc pulse number dependent statistics over noise realizations
    for (std::size_t m = 1; m < PulseTimes_maxsize; m++){
      for (std::size_t k = 0; k < nNoiseRel; k++) mean_delta_t[m] += delta_t[k][m] / nNoiseRel; 
      for (std::size_t k = 0; k < nNoiseRel; k++) var_delta_t[m] += (mean_delta_t[m] - delta_t[k][m]) * (mean_delta_t[m] - delta_t[k][m]) / nNoiseRel;
    }
    
//     for (std::size_t k = 0; k < nNoiseRel; k++) cout << delta_t[k][PulseTimes_maxsize-1] << endl;
    cout << "mean_delta_t: " << mean_delta_t[PulseTimes_maxsize-1] << endl;
    cout << "var_delta_t: " << var_delta_t[PulseTimes_maxsize-1] << endl;
    cout << "std dev delta_t/m: " << sqrt(var_delta_t[PulseTimes_maxsize-1]/(PulseTimes_maxsize-1)) << endl;
    
    //output pulse number dependent mean, var and LTTJ
    std::string str_out_LTTJ_jitter_evo = katana::getCmdOption(argv, argv+argc, "-str_out_LTTJ_jitter_evo" , "data/LTTJ/LTTJ_jitter_evo");
    std::ofstream LTTJ_jitter_evo (str_out_LTTJ_jitter_evo + str_out_LTTJ_suffix, std::ios::out | std::ios::trunc);
    LTTJ_jitter_evo.precision(12);
    for (std::size_t m = 1; m < PulseTimes_maxsize;m++) LTTJ_jitter_evo << m << "\t" << mean_delta_t[m] << "\t" << var_delta_t[m] << "\t" << sqrt(var_delta_t[m]/m) << endl;
  
  }
  
  
  
  //J_G sweep with T01
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_G_with_T01" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.J_G;
    str_sweep_parameter = "J_G";
    PtrToSecPar = &p.T01;
    str_sweep_sec_parameter = "T01";
    
    sweepMLL();
  }
  
  //T01 sweep with JG
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_TO1_with_J_G" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.T01;
    str_sweep_parameter = "T01";
    PtrToSecPar = &p.J_G;
    str_sweep_sec_parameter = "J_G";
    
    sweepMLL();
  }
  
  //J_G sweep with J_Q
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_G_with_J_Q" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.J_G;
    str_sweep_parameter = "J_G";
    PtrToSecPar = &p.J_Q;
    str_sweep_sec_parameter = "J_Q";
    
    sweepMLL();
  }
  
  //J_G sweep with nStr
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_G_with_nStr" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.J_G;
    str_sweep_parameter = "J_G";
    PtrToSecPar = &p.nStr;
    str_sweep_sec_parameter = "nStr";
    
    sweepMLL();
  }
  
  //J_G sweep with delta
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_G_with_delta" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.J_G;
    str_sweep_parameter = "J_G";
    PtrToSecPar = &p.delta;
    str_sweep_sec_parameter = "delta";
    
    sweepMLL();
  }
  
  //J_G sweep with Q_0
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_J_G_with_Q_0" , false)){
  
    tn = 0;
    Xhist.setToCnst(1E-6);
    set_Q_0 = true;
    
    PtrToSweepPar = &p.J_G;
    str_sweep_parameter = "J_G";
    PtrToSecPar = &p.Q_0;
    str_sweep_sec_parameter = "Q_0";
    
    sweepMLL();
  }
  
  //alpha_G sweep with alpha_Q
  if(katana::getCmdOption_bool(argv, argv+argc, "-sweep_alpha_G_with_alpha_Q" , false)){

    tn = 0;
    Xhist.setToCnst(1E-6);
    
    PtrToSweepPar = &p.alpha_G;
    str_sweep_parameter = "alpha_G";
    PtrToSecPar = &p.alpha_Q;
    str_sweep_sec_parameter = "alpha_Q";
    
    sweepMLL();
  }

  

  
  
  return 0;
} 
