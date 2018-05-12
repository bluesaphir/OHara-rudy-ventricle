//ventr
#include <stdio.h>

#define NTABLE 1000
#define ETMIN -100.
#define ETMAX 100.

typedef double REAL;

struct Table { int n1; float frac; } T1;
struct Tlocal { REAL to; REAL uo; REAL tn; REAL un; };

struct State {
    double v;
    double nai;
    double nass;
    double ki;
    double kss;
    double cai;
    double cass;
    double cansr;
    double cajsr;
    double m;
    double hf;
    double hs;
    double j;
    double hsp;
    double jp;
    double mL;
    double hL;
    double hLp;
    double a;
    double iF;
    double iS;
    double ap;
    double iFp;
    double iSp;
    double d;
    double ff;
    double fs;
    double fcaf;
    double fcas;
    double jca;
    double nca;
    double ffp;
    double fcafp;
    double xrf;
    double xrs;
    double xs1;
    double xs2;
    double xk1;
    double Jrelnp;
    double Jrelp;
    double CaMKt;
} *CurrentState, *NextState;

int iter;
REAL ee;

// Coefficients  for ventricle 10 23 12+/ 4 17 16 +/5 23 16+/ 11 2 16+/2 19 17+/3 11 17/ 9 30 16
//'SCN5A', 'KCNA4', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'SLC8A1', 'ATP1A1' ,'KCNJ11' ,'ATP2B4' ,'ATP2A2', 'CALM1', 'RYR2', 'CAMK2D'
	
/*ref

double const c_gnal = 1;
double const c_gto =1; 
double const c_pca = 1;
double const c_gkr = 1;
double const c_gks =1;
double const c_gk1 =1;
double const c_gncx =1;
double const c_pnak =1;
double const c_gkb = 1;
double const c_gkatp = 1;
double const c_gpca = 1;
double const c_jrelinf =1;
double const c_jup = 1;
double const c_CMDN = 1;
double const c_arel = 1;
double const c_CaMKo =1; */
	
						
// 10 23 12						
/*double const c_gnal = 0.8;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =1.102;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 1.064;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 1.386; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 0.727; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =0.852; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx =0.793;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 3.790; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 0.835; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 1.984; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1.;//0.881; 
double const c_jup = 0.906;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 0.814;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 0.230; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo =0.897;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */

/*4 17 16 						
double const c_gnal = 1.019;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =0.894;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 0.705;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 0.847; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 0.883; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =0.686; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx =0.792;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 1.090; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 1.184; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 0.999; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1.;//0.881; 
double const c_jup = 0.827;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 0.978;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 1.322; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo =1.012;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */
						
/*5 23 16 
double const c_gnal = 0.762;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =0.793;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 0.979;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 0.999; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 0.696; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =1.120; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx =1.113;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 0.613; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 0.691; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 0.458; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1.;//0.881; 
double const c_jup = 0.942;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 0.993;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 1.136; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo =0.922;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */
						
//11_2_16
						
double const c_gnal = 1.191;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =0.853;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 0.982;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 0.977; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 1.125; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =1.184; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx = 1.150;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 0.727; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 1.175; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 0.668; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1.;//0.881; 
double const c_jup = 1.25;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 1.118;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 1.159; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo =0.882;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */																
									
/*2_19_17
						
double const c_gnal = 1.167;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =1.283;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 1.317;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 2.137; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 1.496; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =1.098; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx = 0.990;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 1.473; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 2.658; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 2.632; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1;//0.881; 
double const c_jup = 0.75;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 0.60;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 0.98; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo = 2.25;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */
						
						
/*3_11_17
						
double const c_gnal = 1.057;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =0.877;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 1.068;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 0.698; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 1.249; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =1.474; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx = 1.234;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 0.709; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 0.999; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 0.994; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1.;//0.881; 
double const c_jup = 1.41;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 1.06;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 1.51; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo = 1.14;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */	
						 	
/*9_30_16

double const c_gnal = 1.003;//0.776; SCN5A 0.80-1.019-0.761-1.191-1.167-1.056-1.002
double const c_gto =1.196;//0.824; for fast epicard KCND3 1.119-1.229-0.905-1.013-0.902-0.759-1.070
						//        ventr SLOW endocard for KCNA4 1.102-0.894-0.792-0.853-1.283-0.877-1.196
double const c_pca = 0.905;//0.941; CACNA1C 1.064-0.705-0.978-0.982-1.317-1.068-0.905
double const c_gkr = 1.091; // galina KCNH2 1.386-0.847-0.998-0.977-2.137-0.697-1.091
double const c_gks = 0.822; // galina KCNQ1 0.727-0.883-0.695-1.125-1.496-1.249-0.822
double const c_gk1 =1.059; // galina KCNJ2 0.852-0.686-1.119-1.183-1.098-1.474-1.059
double const c_gncx = 0.927;//1.101; SLC8A1 0.793-0.792-1.113-1.149-0.990-1.234-0.927
double const c_pnak = 1.388; //0.711; ATP1A1 3.790-1.090-0.612-0.726-1.473-0.708-1.388
double const c_gkb = 1; //0.473;
double const c_gkatp = 1.115; //galina KCNJ11 0.835-1.184-0.691-1.174-2.658-0.999-1.115
double const c_gpca = 0.896; //galina ATP2B4 1.984-0.998-0.458-0.668-2.632-0.993-0.896
double const c_jrelinf =1.;//0.881; 
double const c_jup = 0.917;//1.175; ATP2A2 0.91-0.83-0.94-1.25-0.75-1.41-0.92
double const c_CMDN = 1.04;//1.184; CALM1 0.81-0.98-0.99-1.12-0.60-1.06-1.04
double const c_arel = 0.67; //RYR2 0.23-1.32-1.35-1.16-0.98-1.51-0.67
double const c_CaMKo = 1.15;// (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
						// CAMK2D 0.9-1.01-0.92-0.88-2.25-1.14-1.15 */	
						 	




//iso coefficients
double const c_INa_dV_act=1.046;
double const c_INa_dV_inact=0.543;
double const c_ICaL_dV_act=0.959;
double const c_ICaL_dV_inact=1.228;
double const c_Knai=1.026;
double const c_GKb_iso=0.617;
double const c_SERCA=0.995;
double const c_GNa_ISO=1.012;
double const c_Tnl=1.155;
double const c_IKs_ISO=0.958;//1.652;
double const c_txs_1_iso=0.459;//1.987;
double const c_a_rel=1.051;
double const c_tau_rel=0.993;

//constants
double const nao=140.0;//extracellular sodium in mM
double const cao=1.8;//extracellular calcium in mM
double const ko=5.4;//extracellular potassium in mM

//buffer paramaters
double const BSRmax=0.047;
double const KmBSR=0.00087;
double const BSLmax=1.124;
double const KmBSL=0.0087;
double const cmdnmax=0.05;
double const kmcmdn=0.00238;
double const trpnmax=0.07;
double const kmtrpn=0.0005;
double const csqnmax=10.0;
double const kmcsqn=0.8;

//CaMK paramaters
double const aCaMK=0.05;
double const bCaMK=0.00068;
double const CaMKo=c_CaMKo*0.05;
double const KmCaM=0.0015;
double const KmCaMK=0.15;

//physical constants
double const R=8314.0;
double const T=310.0;
double const F=96485.0;

//cell geometry
double const L=0.01;
double const rad=0.0011;
double const vcell=1000*3.14*rad*rad*L;
double const Ageo=2*3.14*rad*rad+2*3.14*rad*L;
double const Acap=2*Ageo;
double const vmyo=0.68*vcell;
double const vmito=0.26*vcell;
double const vsr=0.06*vcell;
double const vnsr=0.0552*vcell;
double const vjsr=0.0048*vcell;
double const vss=0.02*vcell;

//introduce varaibles for reversal potentials, currents, fluxes, and CaMK
double ENa,EK,EKs;
double dss,fss;
double INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist,IKatp;
double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
double CaMKa,CaMKb;

//introduce APD, timing, and counting parameters
int APD_flag=0;
double APD;
double t_vdot_max;
double vrest;
double vo = -87.5;
double dt=0.005;
double t0=0.;
double t=0;
double dto;
double vdot_old;
double vdot=0;
double vdot_max;
int p=1;
int n_stim=0;
int Count=0;


const double amp = - 80;//-250; //-250;//250;//stimulus amplitude in uA/uF
const double start = 0;//start time of the stimulus, relative to each beat
const double duration = 0.5;//duration of teh stimulus in ms

double g_gap_junc=5.0;
