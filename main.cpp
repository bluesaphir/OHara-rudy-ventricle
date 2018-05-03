#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "sa.h"
#include <time.h>
#include <assert.h>

#define TABLE 1

#define READIS 1 // 9.03.2017
#define WRITEIS 0 // 9.03.2017
#define CHAIN 0 // 28.06.2017
#define ADAPTIVE_STEP 1 // 06.07.2017
#define IKR_GRANDI 0 //17.07.2017 IKr formulation by Grandi


//Article: Simulation of Action Potentials From Metabolically Impaired Cardiac Myocytes Role of ATP-Sensitive K+ Current
//Article: COMPUTER MODEL OF THE EFFECTS OF PINACIDIL ON ATP-SENSITIVE POTASSIUM CURRENT
//Authors: José M. Ferrero, Javier Sáiz, José M. Ferrero, Nitish V. Thakor

#define IKATP 0 // ATP-sensitive potassium current


/*void revpots(struct State *CurrentState);//compute reversal potentials
void RGC(struct State *CurrentState);//compute rates, gates, and currents
void stimulus(struct State *CurrentState);//determine the value for the periodic stimulus
void voltage(struct State *CurrentState);//calculate teh new membrane voltage
void dVdt_APD(struct State *CurrentState);//caluculate voltage derivative and APD90
void FBC(struct State *CurrentState);//calculate fluxes, buffers, and concentrations
void open(); // 3.03.2017
void close(); // 3.03.2017
void stateInitialization(struct State *CurrentState);// 9.03.2017*/

const double CL = 1000;//pacing cycle length
const double np = 100.;
const double ft = np * CL;//1000*CL;//final time
const int skip = 1./(dt);//number of timesetps to skip in sampling of data in output file
const double safetime = 25.0;//time from the beginning of each beat during which dt is fixed to small values
const double beatssave = 100;//number of beats to save in the output

const int celltype = 0;  //endo = 0, epi = 1, M = 2
const int iso = 0;
FILE*t_,*v_,*nai_,*nass_,*ki_,*kss_,*cai_,*cass_,*cansr_,*cajsr_,*Jrel_,*CaMKt_,*Jup_,*Jtr_,*Jdiff_,*JdiffNa_,*JdiffK_,*Jleak_,*INa_,*INaL_,*Ito_,*ICaL_,*ICaNa_,*ICaK_,*IKr_,*IKs_,*IK1_,*INaCa_i_,*INaCa_ss_,*INaCa_,*INaK_,*IKb_,*INab_,*IpCa_,*ICab_,*Ist_,*dt_,*APD_, *ENa_, *EK_, *EKs_,*f_,*d_,*nca_,*fca_,*jca_,*fp_,*fcap_,*tff_,*tfs_,*tffp_,*dss_, *fss_,*Ical_V;


void tss_xs1(REAL ee, double *txs1, double *xs1ss)
{
     *xs1ss=1.0/(1.0+exp((-(ee + 11.60))/8.932));
     *txs1=817.3+1.0/(2.326e-4*exp((ee + 48.28)/17.80)+0.001292*exp((-(ee + 210.0))/230.0));
    if (iso==1)
    {
        *txs1*=5.475*c_txs_1_iso;//0.6;
    }
    
}

void tss_xs2(REAL ee, double *txs2)
{
     *txs2=1.0/(0.01*exp((ee - 50.0)/20.0)+0.0193 * exp((-(ee + 66.54))/31.0));
}

void tss_m(REAL ee, double *tm , double  *mss){
    if (iso==1)
    {
        *mss=1.0/(1.0+exp((-(ee + 39.57+3.7*c_INa_dV_act))/9.871));
    }
    else  *mss=1.0/(1.0+exp((-(ee + 39.57))/9.871));
    
    *tm=1.0/(6.765 * exp((ee + 11.64)/34.77) + 8.552 * exp(-(ee + 77.42)/5.955));
}

void tss_hf(REAL ee, double  *thf, double *hss){
    if (iso==1)
    {
        *hss=1.0/(1+exp((ee+82.90+4.9*c_INa_dV_inact)/6.086));
        //hss=1.0/(1+exp((CurrentState[z].v+82.90+5.0)/6.086));
    }
    else  *hss=1.0/(1+exp((ee+82.90)/6.086));
    
    *thf=1.0/(1.432e-5 * exp(-(ee + 1.196)/6.285)+6.149 * exp((ee + 0.5096)/20.27));
}

void tss_hs(REAL ee, double  *ths){
    *ths=1.0/(0.009794*exp(-(ee + 17.95)/28.05)+0.3343*exp((ee + 5.730)/56.66));
}

void tss_j(REAL ee, double *tj){
    *tj=2.038+1.0/(0.02136*exp(-(ee +100.6)/8.281)+0.3052*exp((ee + 0.9941)/38.45));
    
}

void tss_hsp(REAL ee, double  *hssp){
    if (iso==1)
    {
        //hssp=1.0/(1+exp((CurrentState[z].v + 89.1+5.0)/6.086));
        *hssp=1.0/(1+exp((ee + 89.1+4.9*c_INa_dV_inact)/6.086));
    }
    else *hssp=1.0/(1+exp((ee + 89.1)/6.086));
}

void tss_mL(REAL ee, double *mLss){
    *mLss=1.0/(1.0 + exp((-(ee + 42.85))/5.264));
}

void tss_hL(REAL ee, double *hLss){
    *hLss=1.0/(1.0+exp((ee + 87.61)/7.488));
}

void tss_hLp(REAL ee, double *hLssp){
    *hLssp=1.0/(1.0+exp((ee +93.81)/7.488));
}

void tss_a(REAL ee, double *ta, double *ass){
    *ass=1.0/(1.0+exp((-(ee - 14.34))/14.82));
    *ta=1.0515/(1.0/(1.2089*(1.0+exp(-(ee -18.4099)/29.3814)))+3.5/(1.0+exp((ee + 100.0)/29.3814)));
}

void tss_iF(REAL ee, double *tiF,double *iss, double *tiFp){
    *iss=1.0/(1.0+exp((ee + 43.94)/5.711));
    double delta_epi;
    if (celltype==1)
    {
        delta_epi=1.0-(0.95/(1.0+exp((ee + 70.0)/5.0)));
    }
    else
    {
        delta_epi=1.0;
    }
    *tiF=4.562+1/(0.3933*exp((-(ee +100.0))/100.0)+0.08004*exp((ee + 50.0)/16.59));
    *tiF *= delta_epi;
    
    double dti_develop=1.354+1.0e-4/(exp((ee - 167.4)/15.89)+exp(-(ee - 12.23)/0.2154));
    double dti_recover=1.0-0.5/(1.0+exp((ee + 70.0)/20.0));
    *tiFp = dti_develop*dti_recover * (*tiF);
}

void tss_iS(REAL ee, double *tiS, double *tiSp){
    double delta_epi;
    if (celltype==1)
    {
        delta_epi=1.0-(0.95/(1.0+exp((ee + 70.0)/5.0)));
    }
    else
    {
        delta_epi=1.0;
    }
    *tiS=23.62+1/(0.001416*exp((-(ee + 96.52))/59.05)+1.780e-8*exp((ee + 114.1)/8.079));
    *tiS *= delta_epi;
    
    double dti_develop = 1.354+1.0e-4/(exp((ee - 167.4)/15.89)+exp(-(ee - 12.23)/0.2154));
    double dti_recover = 1.0-0.5/(1.0+exp((ee + 70.0)/20.0));
    *tiSp = dti_develop * dti_recover * (*tiS);
}

void tss_ap(REAL ee, double *assp){
    *assp=1.0/(1.0+exp((-(ee -24.34))/14.82));
}

void tss_d(REAL ee, double *td, double *dss){
    if (iso==1)
    {
        // dss=1.0/(1.0+exp((-(CurrentState[z].v+3.940+16.0))/4.230));
        // dss=1.0/(1.0+exp((-(CurrentState[z].v+3.940+9.9))/4.230));
        *dss=1.0/(1.0+exp((-(ee+3.940+9.9*c_ICaL_dV_act))/4.230));
    }
    else *dss=1.0/(1.0+exp((-(ee+3.940))/4.230));
    
    *td=0.6+1.0/(exp(-0.05*(ee +6.0))+exp(0.09*(ee +14.0)));
}

void tss_ff(REAL ee, double *tff, double *fss){
    if (iso==1)
    {
        //fss=1.0/(1.0+exp((CurrentState[z].v + 19.58+2.8)/3.696));
        *fss=1.0/(1.0+exp((ee + 19.58+8.0*c_ICaL_dV_inact)/3.696));
    }
    else *fss=1.0/(1.0+exp((ee + 19.58)/3.696));
    
    *tff=7.0+1.0/(0.0045*exp(-(ee + 20.0)/10.0)+0.0045*exp((ee + 20.0)/10.0));
}

void tss_fs(REAL ee, double *tfs){
    *tfs=1000.0+1.0/(0.000035*exp(-(ee+5.0)/4.0)+0.000035*exp((ee + 5.0)/6.0));
}

void tss_fcaf(REAL ee, double *tfcaf){
    *tfcaf=7.0+1.0/(0.04*exp(-(ee - 4.0)/7.0)+0.04 * exp((ee - 4.0)/7.0));
}

void tss_fcas(REAL ee, double *tfcas, double *Afcaf){
    *tfcas=100.0+1.0/(0.00012 * exp(-ee /3.0) + 0.00012 * exp(ee /7.0));
    *Afcaf=0.3+0.6/(1.0+exp((ee - 10.0)/10.0));
}

void tss_xrs(REAL ee, double *txrs, double *xrss, double *rkr){
    #if IKR_GRANDI
        *xrss =   1.0/(1.0 + exp(-(ee  + 10.0)/5.0));
        *txrs  = 550.0/(1.0 + exp(-(ee + 22.0)/9.0)) * 6.0/(1.0 + exp((ee + 11.0)/9.0)) + 230.0/(1.0 + exp((ee + 40.0)/20.0));
        *rkr  =   1.0/(1.0 + exp((ee + 74.0)/24.0));
    #else
        *xrss= 1.0/(1.0+exp((-(ee + 8.337))/6.789));
        *txrs= 1.865+1.0/(0.06629*exp((ee -34.70)/7.355)+1.128e-5 * exp((-(ee - 29.74))/25.94));
        *rkr = 1.0/(1.0+exp((ee + 55.0)/75.0))*1.0/(1.0+exp((ee - 10.0)/30.0));
    #endif
}

void tss_xrf(REAL ee, double *txrf, double *Axrf){
    *txrf=12.98+1.0/(0.3652*exp((ee - 31.66)/3.869)+4.123e-5 * exp((-(ee -47.78))/20.38));
    *Axrf=1.0/(1.0+exp((ee +54.81)/38.21));
}

void tss_xk1(REAL ee, double *txk1, double *xk1ss,double *rk1){
    *xk1ss=1.0/(1.0+exp(-(ee + 2.5538 * ko + 144.59)/(1.5692 * ko + 3.8115)));
    *txk1=122.2/(exp((-(ee + 127.2))/20.36)+exp((ee + 236.8)/69.33));
    *rk1=1.0/(1.0+exp((ee +105.8-2.6*ko)/9.493));
}

clock_t Ina_start, Ina_stop, Ina_late_start, Ina_late_stop, Ito_start, Ito_stop, Ical_start, Ical_stop, Ikr_start, Ikr_stop, Iks_start, Iks_stop, Ik1_start, Ik1_stop, Inaca_start, Inaca_stop, Inak_stop, Inak_start, Inab_stop, Inab_start, ical_test_start, ical_test_stop, ina_t_start, ina_t_stop, ina_wt_start, ina_wt_stop;

void RGC(struct State *CurrentState,  struct State *NextState, int z,  int chain_length)
{
    
    
    for (z=0; z<chain_length; z++)
    {
        FILE *foutput;
        foutput = fopen("foutput.txt","w");
        ical_test_start = clock();
	    CurrentState[z]=NextState[z];
	    //---------old: revpots()---------//
	    ENa=(R*T/F)*log(nao/CurrentState[z].nai);
	    EK=(R*T/F)*log(ko/CurrentState[z].ki);
	    EKs=(R*T/F)*log((ko + 0.01833 * nao)/(CurrentState[z].ki + 0.01833 * CurrentState[z].nai));

	    //-------------RGC()--------------//
        
        /* 1. Calculate INa_fast*/
        Ina_start = clock();
        #if TABLE
            CaMKb = CaMKo * (1.0 - CurrentState[z].CaMKt)/(1.0+KmCaM/CurrentState[z].cass);
            CaMKa = CaMKb + CurrentState[z].CaMKt;
            double vffrt = CurrentState[z].v * F * F/(R*T);
            double vfrt = CurrentState[z].v * F/(R*T);
            double kmtrpn=0.0005;
        
            if (iso==1)
            {
                kmtrpn*=1.6*c_Tnl;
            }
        
            static double mss_t[NTABLE], tm_t[NTABLE], hss_t[NTABLE], thf_t[NTABLE], ths_t[NTABLE], tj_t[NTABLE], hssp_t[NTABLE];
        
            iter = 0;
        
            static int first_inaf = 1;
            for ( ; first_inaf; first_inaf = 0){
                for (iter = -1; ++iter<NTABLE; )
                {
                    ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                    tss_m(ee, tm_t+iter, mss_t+iter);
                    tss_hf(ee, thf_t+iter, hss_t+iter);
                    tss_hs(ee, ths_t+iter);
                    tss_j(ee, tj_t+iter);
                    tss_hsp(ee, hssp_t+iter);
                }
            }
        
            register int iii = T1.n1;
        
            ina_t_start = clock();
            double mss = mss_t[iii]+T1.frac*(mss_t[iii+1] - mss_t[iii]);
            double tm = tm_t[iii]+T1.frac*(tm_t[iii+1] - tm_t[iii]);
            double thf = thf_t[iii]+T1.frac*(thf_t[iii+1] - thf_t[iii]);
            double hss = hss_t[iii]+T1.frac*(hss_t[iii+1] - hss_t[iii]);
            double ths = ths_t[iii]+T1.frac*(ths_t[iii+1] - ths_t[iii]);
            double tj = tj_t[iii]+T1.frac*(tj_t[iii+1] - tj_t[iii]);
            double hssp = hssp_t[iii]+T1.frac*(hssp_t[iii+1] - hssp_t[iii]);
            ina_t_stop = clock();
        
        #else
            CaMKb = CaMKo * (1.0 - CurrentState[z].CaMKt)/(1.0+KmCaM/CurrentState[z].cass);
            CaMKa = CaMKb + CurrentState[z].CaMKt;
            double vffrt = CurrentState[z].v * F * F/(R*T);
            double vfrt = CurrentState[z].v * F/(R*T);
            double kmtrpn=0.0005;
        
            if (iso==1)
            {
                kmtrpn*=1.6*c_Tnl;
            }
        
            double hss, hssp;
    
            double mss;
            if (iso==1)
            {
                mss=1.0/(1.0+exp((-(CurrentState[z].v + 39.57+3.7*c_INa_dV_act))/9.871));
            }
        
            else  mss=1.0/(1.0+exp((-(CurrentState[z].v + 39.57))/9.871));
        
            double tm=1.0/(6.765 * exp((CurrentState[z].v + 11.64)/34.77) + 8.552 * exp(-(CurrentState[z].v + 77.42)/5.955));
        
            if (iso==1)
            {
                hss=1.0/(1+exp((CurrentState[z].v+82.90+4.9*c_INa_dV_inact)/6.086));
                //hss=1.0/(1+exp((CurrentState[z].v+82.90+5.0)/6.086));
            }
            else  hss=1.0/(1+exp((CurrentState[z].v+82.90)/6.086));
	   
            double thf=1.0/(1.432e-5 * exp(-(CurrentState[z].v + 1.196)/6.285)+6.149 * exp((CurrentState[z].v + 0.5096)/20.27));
            double ths=1.0/(0.009794*exp(-(CurrentState[z].v + 17.95)/28.05)+0.3343*exp((CurrentState[z].v + 5.730)/56.66));
            double tj=2.038+1.0/(0.02136*exp(-(CurrentState[z].v +100.6)/8.281)+0.3052*exp((CurrentState[z].v + 0.9941)/38.45));
        
            if (iso==1)
            {
              //hssp=1.0/(1+exp((CurrentState[z].v + 89.1+5.0)/6.086));
                hssp=1.0/(1+exp((CurrentState[z].v + 89.1+4.9*c_INa_dV_inact)/6.086));
            }
            else hssp=1.0/(1+exp((CurrentState[z].v + 89.1)/6.086));
        #endif
        
        double Ahf=0.99;
        double Ahs=1.0-Ahf;
        double jss=hss;
        double thsp=3.0*ths;
        double tjp=1.46*tj;
        
        
        NextState[z].m = mss - (mss - CurrentState[z].m)*exp(-dt/tm);
        NextState[z].hf = hss - (hss - CurrentState[z].hf) * exp(-dt/thf);
        NextState[z].hs = hss - (hss - CurrentState[z].hs) * exp(-dt/ths);
        NextState[z].j = jss - ( jss - CurrentState[z].j) * exp(-dt/tj);
        NextState[z].hsp = hssp -(hssp - CurrentState[z].hsp) * exp(-dt/thsp);
        NextState[z].jp = jss - (jss - CurrentState[z].jp) * exp(-dt/tjp);
        
        double h = Ahf * CurrentState[z].hf + Ahs * CurrentState[z].hs;
        double hp = Ahf * CurrentState[z].hf + Ahs * CurrentState[z].hsp;
        
        double GNa=75.;
        if (iso==1)
        {
            GNa*=2.35*c_GNa_ISO;//2.5;//2.2;//2.7;
        }
        else GNa*=c_gnal;
        double fINap=(1.0/(1.0+KmCaMK/CaMKa));
        
        INa = GNa * (CurrentState[z].v-ENa) * CurrentState[z].m * CurrentState[z].m * CurrentState[z].m * ((1.0-fINap) * h * CurrentState[z].j + fINap * hp * CurrentState[z].jp);
        
        Ina_stop = clock();
        fprintf(foutput,"Ina all functions: %e\n",(double)(ina_wt_stop - ina_wt_start) / (CLOCKS_PER_SEC));
        fprintf(foutput,"Ina: %e\n",(double)(Ina_stop - Ina_start) / (CLOCKS_PER_SEC));
        
        /* 2. Calculate INa_late*/
        Ina_late_start = clock();
        #if TABLE
            static double mLss_t[NTABLE], hLss_t[NTABLE], hLssp_t[NTABLE];
        
            iter = 0;
        
            static int first_inal = 1;
            for ( ; first_inal; first_inal = 0){
                for (iter = -1; ++iter<NTABLE; )
                {
                    ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                    tss_mL(ee, mLss_t+iter);
                    tss_hL(ee, hLss_t+iter);
                    tss_hLp(ee, hLssp_t+iter);
                }
            }
        
            double mLss = mLss_t[iii]+T1.frac*(mLss_t[iii+1] - mLss_t[iii]);
            double hLss = hLss_t[iii]+T1.frac*(hLss_t[iii+1] - hLss_t[iii]);
            double hLssp = hLssp_t[iii]+T1.frac*(hLssp_t[iii+1] - hLssp_t[iii]);
        
        #else
            double mLss=1.0/(1.0+exp((-(CurrentState[z].v + 42.85))/5.264));
            double hLss=1.0/(1.0+exp((CurrentState[z].v + 87.61)/7.488));
            double hLssp=1.0/(1.0+exp((CurrentState[z].v +93.81)/7.488));
	    #endif
        
        double tmL=tm;
        double thL=200.0;
        double thLp=3.0*thL;
        
        NextState[z].mL = mLss - (mLss - CurrentState[z].mL) * exp(-dt/tmL);
        NextState[z].hL = hLss - (hLss - CurrentState[z].hL) * exp(-dt/thL);
        NextState[z].hLp = hLssp - (hLssp - CurrentState[z].hLp) * exp(-dt/thLp);
        
        double GNaL=0.0075*c_gnal;
        if (celltype==1)
        {
            GNaL*=0.6;
        }
        double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
        INaL = GNaL * (CurrentState[z].v - ENa) * CurrentState[z].mL * ((1.0-fINaLp) * CurrentState[z].hL + fINaLp * CurrentState[z].hLp);
        
        Ina_late_stop = clock();
        fprintf(foutput,"Ina_late: %e\n",(double)(Ina_late_stop - Ina_late_start) / (CLOCKS_PER_SEC));
        //printf("INaL: %lf\n", INaL);
        
        
        /* 3. Calculate Ito*/
        Ito_start = clock();
        #if TABLE
            static double ta_t[NTABLE], ass_t[NTABLE], tiF_t[NTABLE], tiS_t[NTABLE], iss_t[NTABLE], assp_t[NTABLE], tiFp_t[NTABLE], tiSp_t[NTABLE];
        
            iter = 0;
        
            static int first_ito = 1;
            for ( ; first_ito; first_ito = 0){
                for (iter = -1; ++iter<NTABLE; )
                {
                    ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                    tss_a(ee, ta_t+iter, ass_t+iter);
                    tss_iF(ee, tiF_t+iter, iss_t+iter, tiFp_t+iter);
                    tss_iS(ee, tiS_t+iter, tiSp_t+iter);
                    tss_ap(ee, assp_t+iter);
                }
            }
            double ta = ta_t[iii]+T1.frac*(ta_t[iii+1] - ta_t[iii]);
            double ass = ass_t[iii]+T1.frac*(ass_t[iii+1] - ass_t[iii]);
            double tiF = tiF_t[iii]+T1.frac*(tiF_t[iii+1] - tiF_t[iii]);
            double iss = iss_t[iii]+T1.frac*(iss_t[iii+1] - iss_t[iii]);
            double tiS = tiS_t[iii]+T1.frac*(tiS_t[iii+1] - tiS_t[iii]);
            double assp = assp_t[iii]+T1.frac*(assp_t[iii+1] - assp_t[iii]);
            double tiFp = tiFp_t[iii]+T1.frac*(tiFp_t[iii+1] - tiFp_t[iii]);
            double tiSp = tiSp_t[iii]+T1.frac*(tiSp_t[iii+1] - tiSp_t[iii]);

        
        #else
            double ass=1.0/(1.0+exp((-(CurrentState[z].v - 14.34))/14.82));
            double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(CurrentState[z].v -18.4099)/29.3814)))+3.5/(1.0+exp((CurrentState[z].v + 100.0)/29.3814)));
        
            double iss=1.0/(1.0+exp((CurrentState[z].v + 43.94)/5.711));
            double delta_epi;
            if (celltype==1)
            {
                delta_epi=1.0-(0.95/(1.0+exp((CurrentState[z].v + 70.0)/5.0)));
            }
            else
            {
                delta_epi=1.0;
            }
            double tiF=4.562+1/(0.3933*exp((-(CurrentState[z].v +100.0))/100.0)+0.08004*exp((CurrentState[z].v + 50.0)/16.59));
            double tiS=23.62+1/(0.001416*exp((-(CurrentState[z].v + 96.52))/59.05)+1.780e-8*exp((CurrentState[z].v + 114.1)/8.079));
            tiF *= delta_epi;
            tiS *= delta_epi;
        
            double assp=1.0/(1.0+exp((-(CurrentState[z].v -24.34))/14.82));
            double dti_develop=1.354+1.0e-4/(exp((CurrentState[z].v - 167.4)/15.89)+exp(-(CurrentState[z].v - 12.23)/0.2154));
            double dti_recover=1.0-0.5/(1.0+exp((CurrentState[z].v + 70.0)/20.0));
            double tiFp=dti_develop*dti_recover*tiF;
            double tiSp=dti_develop*dti_recover*tiS;
        #endif
        
        double AiF=1.0/(1.0+exp((CurrentState[z].v - 213.6)/151.2));
        double AiS=1.0-AiF;
        double i = AiF * CurrentState[z].iF + AiS * CurrentState[z].iS;
        double ip = AiF * CurrentState[z].iFp + AiS * CurrentState[z].iSp;
        
        NextState[z].a = ass - (ass - CurrentState[z].a) * exp(-dt/ta);
        NextState[z].iF = iss - (iss - CurrentState[z].iF) * exp(-dt/tiF);
        NextState[z].iS = iss - (iss - CurrentState[z].iS) * exp(-dt/tiS);
        NextState[z].ap = assp -(assp - CurrentState[z].ap) * exp(-dt/ta);
        NextState[z].iFp = iss - (iss - CurrentState[z].iFp) * exp(-dt/tiFp);
        NextState[z].iSp = iss - (iss - CurrentState[z].iSp) * exp(-dt/tiSp);
        
        double Gto=0.02*c_gto;
        if (celltype==1)
        {
            Gto*=4.0;
        }
        if (celltype==2)
        {
            Gto*=4.0;
        }
        
        double fItop=(1.0/(1.0+KmCaMK/CaMKa));
        Ito = Gto * (CurrentState[z].v - EK) * ((1.0 - fItop) * CurrentState[z].a * i + fItop * CurrentState[z].ap * ip);
        Ito_stop = clock();
        fprintf(foutput, "Ito: %e\n",(double)(Ito_stop - Ito_start) / (CLOCKS_PER_SEC));
        
        /* 4-6. Calculate ICaL, ICaNa, ICaK */
        Ical_start = clock();
        #if TABLE
            static double dss_t[NTABLE], td_t[NTABLE], tff_t[NTABLE], tfs_t[NTABLE], fss_t[NTABLE], tfcaf_t[NTABLE], tfcas_t[NTABLE], Afcaf_t[NTABLE];
        
            iter = 0;
        
            static int first_ica = 1;
            for ( ; first_ica; first_ica = 0){
                for (iter = -1; ++iter<NTABLE; )
                {
                    ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                    tss_d(ee, td_t+iter, dss_t+iter);
                    tss_ff(ee, tff_t+iter, fss_t+iter);
                    tss_fs(ee, tfs_t+iter);
                    tss_fcaf(ee, tfcaf_t+iter);
                    tss_fcas(ee, tfcas_t+iter, Afcaf_t+iter);
                }
            }
        
            double td = td_t[iii]+T1.frac*(td_t[iii+1] - td_t[iii]);
            double dss = dss_t[iii]+T1.frac*(dss_t[iii+1] - dss_t[iii]);
            double tff = tff_t[iii]+T1.frac*(tff_t[iii+1] - tff_t[iii]);
            double tfs = tfs_t[iii]+T1.frac*(tfs_t[iii+1] - tfs_t[iii]);
            double fss = fss_t[iii]+T1.frac*(fss_t[iii+1] - fss_t[iii]);
            double tfcaf = tfcaf_t[iii]+T1.frac*(tfcaf_t[iii+1] - tfcaf_t[iii]);
            double tfcas = tfcas_t[iii]+T1.frac*(tfcas_t[iii+1] - tfcas_t[iii]);
            double Afcaf = Afcaf_t[iii]+T1.frac*(Afcaf_t[iii+1] - Afcaf_t[iii]);
        
        #else
    
            if (iso==1)
            {
                dss=1.0/(1.0+exp((-(CurrentState[z].v+3.940+9.9*c_ICaL_dV_act))/4.230));
            }
            else dss=1.0/(1.0+exp((-(CurrentState[z].v+3.940))/4.230));

            double td=0.6+1.0/(exp(-0.05*(CurrentState[z].v +6.0))+exp(0.09*(CurrentState[z].v +14.0)));
            if (iso==1)
            {
                //fss=1.0/(1.0+exp((CurrentState[z].v + 19.58+2.8)/3.696));
                fss=1.0/(1.0+exp((CurrentState[z].v + 19.58+8.0*c_ICaL_dV_inact)/3.696));
            }
            else fss=1.0/(1.0+exp((CurrentState[z].v + 19.58)/3.696));

            double tff=7.0+1.0/(0.0045*exp(-(CurrentState[z].v + 20.0)/10.0)+0.0045*exp((CurrentState[z].v + 20.0)/10.0));
            double tfs=1000.0+1.0/(0.000035*exp(-(CurrentState[z].v+5.0)/4.0)+0.000035*exp((CurrentState[z].v + 5.0)/6.0));
            double tfcaf=7.0+1.0/(0.04*exp(-(CurrentState[z].v - 4.0)/7.0)+0.04 * exp((CurrentState[z].v - 4.0)/7.0));
            double tfcas=100.0+1.0/(0.00012 * exp(-CurrentState[z].v /3.0) + 0.00012 * exp(CurrentState[z].v /7.0));
            double Afcaf=0.3+0.6/(1.0+exp((CurrentState[z].v - 10.0)/10.0));
	    
        #endif
        
        double Aff=0.6;
        double Afs=1.0-Aff;
        double f = Aff * CurrentState[z].ff + Afs * CurrentState[z].fs;
        double fcass=fss;
        double Afcas=1.0-Afcaf;
        double fca = Afcaf * CurrentState[z].fcaf + Afcas * CurrentState[z].fcas;
        double tjca=75.0;
        double tffp=2.5*tff;
        double fp = Aff * CurrentState[z].ffp + Afs * CurrentState[z].fs;
        double tfcafp=2.5*tfcaf;
        double fcap = Afcaf * CurrentState[z].fcafp + Afcas * CurrentState[z].fcas;
        double Kmn=0.002;
        double k2n=1000.0;
        double km2n = CurrentState[z].jca * 1.0;
        double anca=1.0/(k2n/km2n+pow(1.0+Kmn/CurrentState[z].cass,4.0));
        double PhiCaL;
        
        if (CurrentState[z].cass<0.03)  PhiCaL = 4.0 * vffrt * (CurrentState[z].cass * exp(2.0*vfrt) - 0.341 * cao)/(exp(2.0 * vfrt)-1.0);
        else  PhiCaL = 4.0 * vffrt * (0.03 * exp(2.0*vfrt) - 0.341 * cao)/(exp(2.0 * vfrt)-1.0);
        
        double PhiCaNa=1.0*vffrt*(0.75 * CurrentState[z].nass * exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
        double PhiCaK = 1.0 * vffrt * (0.75 * CurrentState[z].kss * exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
        double zca=2.0;
        
        NextState[z].d = dss - (dss - CurrentState[z].d) * exp(-dt/td);
        NextState[z].ff = fss - (fss - CurrentState[z].ff) * exp(-dt/tff);
        NextState[z].fs = fss - (fss - CurrentState[z].fs) * exp(-dt/tfs);
        NextState[z].fcaf = fcass - (fcass - CurrentState[z].fcaf) * exp(-dt/tfcaf);
        NextState[z].fcas = fcass - (fcass - CurrentState[z].fcas) * exp(-dt/tfcas);
        NextState[z].jca = fcass - (fcass - CurrentState[z].jca) * exp(-dt/tjca);
        NextState[z].ffp = fss - (fss - CurrentState[z].ffp) * exp(-dt/tffp);
        NextState[z].fcafp = fcass - (fcass - CurrentState[z].fcafp) * exp(-dt/tfcafp);
        NextState[z].nca = anca * k2n/km2n - (anca * k2n/km2n - CurrentState[z].nca) * exp(-km2n * dt);
        
        double PCa=0.0001*c_pca;
        if (celltype==1)
        {
            PCa*=1.2;
        }
        if (celltype==2)
        {
            PCa*=2.5;
        }
        double PCap=1.1*PCa;
        double PCaNa=0.00125*PCa;
        double PCaK=3.574e-4*PCa;
        double PCaNap=0.00125*PCap;
        double PCaKp=3.574e-4*PCap;
        double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
        
        ICaL=(1.0-fICaLp) * PCa * PhiCaL * CurrentState[z].d * (f * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fca * CurrentState[z].nca) + fICaLp * PCap * PhiCaL * CurrentState[z].d * (fp * (1.0-CurrentState[z].nca) + CurrentState[z].jca * fcap * CurrentState[z].nca);
        ICaNa=(1.0-fICaLp) * PCaNa * PhiCaNa * CurrentState[z].d * (f * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fca * CurrentState[z].nca)+fICaLp * PCaNap * PhiCaNa * CurrentState[z].d * (fp * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fcap * CurrentState[z].nca);
        ICaK=(1.0-fICaLp) * PCaK * PhiCaK * CurrentState[z].d * (f * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fca * CurrentState[z].nca)+fICaLp * PCaKp * PhiCaK * CurrentState[z].d * (fp * (1.0 - CurrentState[z].nca) + CurrentState[z].jca * fcap * CurrentState[z].nca);
        
        Ical_stop = clock();
        fprintf(foutput,"Ical: %e\n",(double)(Ical_stop - Ical_start) / (CLOCKS_PER_SEC));

        /* 7. Calculate IKr */
        Ikr_start = clock();
        #if TABLE
            double GKr = 0.046;
            #if IKR_GRANDI
                static double xrss_t[NTABLE], txrs_t[NTABLE], rkr_t[NTABLE];
        
                iter = 0;
        
                static int first_ikr = 1;
                for ( ; first_ikr; first_ikr = 0){
                    for (iter = -1; ++iter<NTABLE; )
                    {
                        ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                        tss_xrs(ee, txrs_t+iter, xrss_t+iter, rkr_t+iter);
                    }
                }
        
                double xrss = xrss_t[iii]+T1.frac*(xrss_t[iii+1] - xrss_t[iii]);
                double txrs = txrs_t[iii]+T1.frac*(txrs_t[iii+1] - txrs_t[iii]);
                double rkr  = rkr_t[iii]+T1.frac*(rkr_t[iii+1] - rkr_t[iii]);
        
                NextState[z].xrs = xrss - (xrss - CurrentState[z].xrs) * exp(-dt/txrs);
        
                if (celltype==1)
                {
                    GKr*=1.3;
                }
                if (celltype==2)
                {
                    GKr*=0.8;
                }
                IKr=GKr*sqrt(ko/5.4)*CurrentState[z].xrs*rkr*(CurrentState[z].v - EK);
            #else
                static double xrss_t[NTABLE], txrs_t[NTABLE], rkr_t[NTABLE], txrf_t[NTABLE], Axrf_t[NTABLE];
        
                iter = 0;
        
                static int first_ikr = 1;
                for ( ; first_ikr; first_ikr = 0){
                    for (iter = -1; ++iter<NTABLE; )
                    {
                        ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                        tss_xrs(ee, txrs_t+iter, xrss_t+iter, rkr_t+iter);
                        tss_xrf(ee, txrf_t + iter, Axrf_t + iter);
                    }
                }
        
                double xrss = xrss_t[iii]+T1.frac*(xrss_t[iii+1] - xrss_t[iii]);
                double txrs = txrs_t[iii]+T1.frac*(txrs_t[iii+1] - txrs_t[iii]);
                double rkr  = rkr_t[iii]+T1.frac*(rkr_t[iii+1] - rkr_t[iii]);
                double txrf  = txrf_t[iii]+T1.frac*(txrf_t[iii+1] - txrf_t[iii]);
                double Axrf  = Axrf_t[iii]+T1.frac*(Axrf_t[iii+1] - Axrf_t[iii]);
        
                double Axrs=1.0-Axrf;
                double xr = Axrf * CurrentState[z].xrf + Axrs * CurrentState[z].xrs;
        
                NextState[z].xrf = xrss - (xrss - CurrentState[z].xrf) * exp(-dt/txrf);
                NextState[z].xrs = xrss - (xrss - CurrentState[z].xrs) * exp(-dt/txrs);
        
                if (celltype==1)
                {
                    GKr*=1.3;
                }
                if (celltype==2)
                {
                    GKr*=0.8;
                }
                IKr=GKr*sqrt(ko/5.4)*xr*rkr*(CurrentState[z].v - EK);
            #endif
        #else
            double GKr = 0.046;
            #if IKR_GRANDI
                double xrss =   1.0/(1.0 + exp(-(CurrentState[z].v  + 10.0)/5.0));
                double txrs  = 550.0/(1.0 + exp(-(CurrentState[z].v + 22.0)/9.0)) * 6.0/(1.0 + exp((CurrentState[z].v + 11.0)/9.0)) + 230.0/(1.0 + exp( (CurrentState[z].v + 40.0)/20.0));
                NextState[z].xrs = xrss - (xrss - CurrentState[z].xrs) * exp(-dt/txrs);
                double rkr  =   1.0/(1.0 + exp((CurrentState[z].v + 74.0)/24.0));
        
                if (celltype==1)
                {
                    GKr*=1.3;
                }
                if (celltype==2)
                {
                    GKr*=0.8;
                }
                IKr=GKr*sqrt(ko/5.4)*CurrentState[z].xrs*rkr*(CurrentState[z].v - EK);
        
            #else
                //    double xrss =   1.0/(1.0 + exp(-(CurrentState[z].v  + 10.0)/5.0));
                double xrss=1.0/(1.0+exp((-(CurrentState->v + 8.337))/6.789));
                double txrf=12.98+1.0/(0.3652*exp((CurrentState->v - 31.66)/3.869)+4.123e-5 * exp((-(CurrentState->v -47.78))/20.38));
                //    double txrs  = 550.0/(1.0 + exp(-(CurrentState[z].v + 22.0)/9.0)) * 6.0/(1.0 + exp((CurrentState[z].v + 11.0)/9.0)) + 230.0/(1.0 + exp( (CurrentState[z].v + 40.0)/20.0));
                double txrs=1.865+1.0/(0.06629*exp((CurrentState->v -34.70)/7.355)+1.128e-5 * exp((-(CurrentState->v - 29.74))/25.94));
                double Axrf=1.0/(1.0+exp((CurrentState->v +54.81)/38.21));
                double Axrs=1.0-Axrf;
                NextState[z].xrf = xrss - (xrss - CurrentState[z].xrf) * exp(-dt/txrf);
                NextState[z].xrs = xrss - (xrss - CurrentState[z].xrs) * exp(-dt/txrs);
                double xr = Axrf * CurrentState[z].xrf + Axrs * CurrentState[z].xrs;
                //    double rkr  =   1.0/(1.0 + exp((CurrentState[z].v + 74.0)/24.0));
                double rkr=1.0/(1.0+exp((CurrentState->v + 55.0)/75.0))*1.0/(1.0+exp((CurrentState->v - 10.0)/30.0));
                                                                     
                                                                     
                if (celltype==1)
                {
                    GKr*=1.3;
                }
                if (celltype==2)
                {
                    GKr*=0.8;
                }
                IKr=GKr*sqrt(ko/5.4)*xr*rkr*(CurrentState[z].v - EK);
            #endif
        #endif
        Ikr_stop = clock();
        fprintf(foutput,"Ikr: %e\n",(double)(Ikr_stop - Ikr_start) / (CLOCKS_PER_SEC));
        
        /* 8. Calculate IKs */
        Iks_start = clock();
        #if TABLE
            static double txs1_t[NTABLE], xs1ss_t[NTABLE], txs2_t[NTABLE];
        
            iter = 0;
        
            static int first_iks = 1;
            for ( ; first_iks; first_iks = 0){
                for (iter = -1; ++iter<NTABLE; )
                {
                    ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                    tss_xs1(ee, txs1_t+iter, xs1ss_t+iter);
                    tss_xs2(ee, txs2_t+iter);
                }
            }
        
            double txs1 = txs1_t[iii]+T1.frac*(txs1_t[iii+1] - txs1_t[iii]);
            double xs1ss = xs1ss_t[iii]+T1.frac*(xs1ss_t[iii+1] - xs1ss_t[iii]);
            double txs2 = txs2_t[iii]+T1.frac*(txs2_t[iii+1] - txs2_t[iii]);
        #else
            double xs1ss=1.0/(1.0+exp((-(CurrentState[z].v + 11.60))/8.932));
            //printf("xs1ss: %e\n", xs1ss);
            double txs1=817.3+1.0/(2.326e-4*exp((CurrentState[z].v + 48.28)/17.80)+0.001292*exp((-(CurrentState[z].v + 210.0))/230.0));
            if (iso==1)
            {
                txs1*=5.475*c_txs_1_iso;//0.6;
            }
        
            double txs2=1.0/(0.01*exp((CurrentState[z].v - 50.0)/20.0)+0.0193 * exp((-(CurrentState[z].v + 66.54))/31.0));
        #endif
        double xs2ss=xs1ss;
        
        NextState[z].xs1 = xs1ss - (xs1ss - CurrentState[z].xs1) * exp(-dt/txs1);
        NextState[z].xs2 = xs2ss - (xs2ss - CurrentState[z].xs2) * exp(-dt/txs2);
        
        double KsCa=1.0+0.6/(1.0+pow(3.8e-5/CurrentState[z].cai,1.4));
        double GKs=0.0034;
        double pdgn=1;
        if (iso==1)
        {
            GKs*=3.2*c_IKs_ISO;
            pdgn*=1.;
        }
        
        if (celltype==1)
        {
            GKs*=1.4;
        }
        IKs = pdgn * GKs * KsCa * CurrentState[z].xs1 * CurrentState[z].xs2 * (CurrentState[z].v-EKs);
        
        Iks_stop = clock();
        fprintf(foutput,"Iks: %e\n",(double)(Iks_stop - Iks_start) / (CLOCKS_PER_SEC));
	    
        /* 9. Calculate IK1 */
        Ik1_start = clock();
        #if TABLE 
            static double txk1_t[NTABLE], xk1ss_t[NTABLE], rk1_t[NTABLE];
        
            iter = 0;
        
            static int first_ik1 = 1;
            for ( ; first_ik1; first_ik1 = 0){
                for (iter = -1; ++iter<NTABLE; )
                {
                    ee = ETMIN + (ETMAX - ETMIN)*iter/(NTABLE-1);
                    tss_xk1(ee, txk1_t+iter, xk1ss_t+iter, rk1_t+iter);
                }
            }
        
            double txk1 = txk1_t[iii]+T1.frac*(txk1_t[iii+1] - txk1_t[iii]);
            double xk1ss = xk1ss_t[iii]+T1.frac*(xk1ss_t[iii+1] - xk1ss_t[iii]);
            double rk1 = rk1_t[iii]+T1.frac*(rk1_t[iii+1] - rk1_t[iii]);
        
        #else
            double xk1ss=1.0/(1.0+exp(-(CurrentState[z].v + 2.5538 * ko + 144.59)/(1.5692 * ko + 3.8115)));
            double txk1=122.2/(exp((-(CurrentState[z].v + 127.2))/20.36)+exp((CurrentState[z].v + 236.8)/69.33));
            double rk1=1.0/(1.0+exp((CurrentState[z].v +105.8-2.6*ko)/9.493));
        #endif
        
        NextState[z].xk1 = xk1ss - (xk1ss - CurrentState[z].xk1) * exp(-dt/txk1);
        double GK1=0.1908;
        if (celltype==1)
        {
            GK1*=1.2;
        }
        if (celltype==2)
        {
            GK1*=1.3;
        }
        IK1 = GK1 * sqrt(ko) * rk1 * CurrentState[z].xk1 * (CurrentState[z].v - EK);
        
        Ik1_stop = clock();
        fprintf(foutput,"Ik1: %e\n",(double)(Ik1_stop - Ik1_start) / (CLOCKS_PER_SEC));
        //printf("IK1: %lf\n", IK1);
        
        /* 10. Calculate INaCa */
        Inaca_start = clock();
	    double kna1=15.0;
	    double kna2=5.0;
	    double kna3=88.12;
	    double kasymm=12.5;
	    double wna=6.0e4;
	    double wca=6.0e4;
	    double wnaca=5.0e3;
	    double kcaon=1.5e6;
	    double kcaoff=5.0e3;
	    double qna=0.5224;
	    double qca=0.1670;
	    double hca = exp((qca * CurrentState[z].v * F)/(R * T));
	    double hna = exp((qna * CurrentState[z].v * F)/(R * T));
	    double h1 = 1 + CurrentState[z].nai / kna3 * (1 + hna);
	    double h2=(CurrentState[z].nai * hna)/(kna3*h1);
	    double h3=1.0/h1;
	    double h4 = 1.0 + CurrentState[z].nai/kna1 * (1+ CurrentState[z].nai/kna2);
	    double h5=CurrentState[z].nai * CurrentState[z].nai/(h4*kna1*kna2);
	    double h6=1.0/h4;
	    double h7=1.0+nao/kna3*(1.0+1.0/hna);
	    double h8=nao/(kna3*hna*h7);
	    double h9=1.0/h7;
	    double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
	    double h11=nao*nao/(h10*kna1*kna2);
	    double h12=1.0/h10;
	    double k1=h12*cao*kcaon;
	    double k2=kcaoff;
	    double k3p=h9*wca;
	    double k3pp=h8*wnaca;
	    double k3=k3p+k3pp;
	    double k4p=h3*wca/hca;
	    double k4pp=h2*wnaca;
	    double k4=k4p+k4pp;
	    double k5=kcaoff;
	    double k6 = h6 * CurrentState[z].cai * kcaon;
	    double k7=h5*h2*wna;
	    double k8=h8*h11*wna;
	    double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	    double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	    double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	    double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	    double E1=x1/(x1+x2+x3+x4);
	    double E2=x2/(x1+x2+x3+x4);
	    double E3=x3/(x1+x2+x3+x4);
	    double E4=x4/(x1+x2+x3+x4);
	    double KmCaAct=150.0e-6;
	    double allo=1.0/(1.0+pow(KmCaAct/CurrentState[z].cai,2.0));
	    double zna=1.0;
	    double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	    double JncxCa=E2*k2-E1*k1;
	    double Gncx=0.0008*c_gncx;
	    if (celltype==1)
	    {
		Gncx*=1.1;
	    }
	    if (celltype==2)
	    {
		Gncx*=1.4;
	    }
	    INaCa_i = 0.8 * Gncx * allo * (zna*JncxNa+zca*JncxCa);
	    
	    h1 = 1 + CurrentState[z].nass/kna3*(1 + hna);
	    h2 = (CurrentState[z].nass * hna)/(kna3 * h1);
	    h3 = 1.0/h1;
	    h4 = 1.0 + CurrentState[z].nass/kna1*(1 + CurrentState[z].nass/kna2);
	    h5 = CurrentState[z].nass * CurrentState[z].nass/(h4*kna1*kna2);
	    h6 = 1.0/h4;
	    h7 = 1.0+nao/kna3*(1.0+1.0/hna);
	    h8 = nao/(kna3*hna*h7);
	    h9 = 1.0/h7;
	    h10= kasymm+1.0+nao/kna1*(1+nao/kna2);
	    h11= nao*nao/(h10*kna1*kna2);
	    h12= 1.0/h10;
	    k1 = h12*cao*kcaon;
	    k2 = kcaoff;
	    k3p= h9*wca;
	    k3pp= h8*wnaca;
	    k3 =k3p+k3pp;
	    k4p=h3*wca/hca;
	    k4pp=h2*wnaca;
	    k4=k4p+k4pp;
	    k5=kcaoff;
	    k6 = h6 * CurrentState[z].cass * kcaon;
	    k7=h5*h2*wna;
	    k8=h8*h11*wna;
	    x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	    x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	    x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	    x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	    E1=x1/(x1+x2+x3+x4);
	    E2=x2/(x1+x2+x3+x4);
	    E3=x3/(x1+x2+x3+x4);
	    E4=x4/(x1+x2+x3+x4);
	    KmCaAct=150.0e-6;
	    allo = 1.0/(1.0 + pow(KmCaAct/CurrentState[z].cass,2.0));
	    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	    JncxCa=E2*k2-E1*k1;
	    INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);
	    
	    INaCa=INaCa_i+INaCa_ss;
        
        Inaca_stop = clock();
        fprintf(foutput,"Inaca: %e\n",(double)(Inaca_stop - Inaca_start) / (CLOCKS_PER_SEC));
	    
        /* 11. Calculate INaK*/
        Inak_start = clock();
        
	    double k1p=949.5;
	    double k1m=182.4;
	    double k2p=687.2;
	    double k2m=39.4;
	    k3p=1899.0;
	    double k3m=79300.0;
	    k4p=639.0;
	    double k4m=40.0;
	    double Knai0=9.073;
	    double Knao0=27.78;
            
	    double delta=-0.1550;
	    double Knai = Knai0 * exp((delta * CurrentState[z].v * F)/(3.0 * R * T));
            if (iso==1)
	    {
		   Knai*=0.72*c_Knai;//0.7;
	    }
	    double Knao = Knao0 * exp(((1.0-delta) * CurrentState[z].v * F)/(3.0 * R * T));
	    double Kki = 0.5;
	    double Kko = 0.3582;
	    double MgADP = 0.05;
	    double MgATP = 9.8;
	    double Kmgatp = 1.698e-7;
	    double H = 1.0e-7;
	    double eP = 4.2;
	    double Khp = 1.698e-7;
	    double Knap = 224.0;
	    double Kxkur = 292.0;
	    double P = eP/(1.0+H/Khp + CurrentState[z].nai/Knap + CurrentState[z].ki/Kxkur);
	    double a1 = (k1p*pow(CurrentState[z].nai/Knai,3.0))/(pow(1.0+CurrentState[z].nai/Knai,3.0)+pow(1.0 + CurrentState[z].ki/Kki,2.0)-1.0);
	    double b1 = k1m*MgADP;
	    double a2 = k2p;
	    double b2 = (k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
	    double a3 = (k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
	    double b3 = (k3m*P*H)/(1.0+MgATP/Kmgatp);
	    double a4 = (k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
	    double b4 = (k4m * pow(CurrentState[z].ki/Kki,2.0))/(pow(1.0 + CurrentState[z].nai/Knai,3.0)+pow(1.0 + CurrentState[z].ki/Kki,2.0)-1.0);
	    
	    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
	    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
	    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
	    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
	    E1=x1/(x1+x2+x3+x4);
	    E2=x2/(x1+x2+x3+x4);
	    E3=x3/(x1+x2+x3+x4);
	    E4=x4/(x1+x2+x3+x4);
	    
	    double zk=1.0;
	    double JnakNa=3.0*(E1*a3-E2*b3);
	    double JnakK=2.0*(E4*b1-E3*a1);
	    double Pnak=30*c_pnak;
	    if (celltype==1)
	    {
		Pnak*=0.9;
	    }
	    if (celltype==2)
	    {
		Pnak*=0.7;
	    }
	    INaK = Pnak * (zna * JnakNa + zk * JnakK);
        
        Inak_stop = clock();
        fprintf(foutput,"Inak: %e\n",(double)(Inak_stop - Inak_start) / (CLOCKS_PER_SEC));
	    
        /* 12-14. Calculate INab, ICab, IpCa*/
        Inab_start = clock();
        
	    double xkb=1.0/(1.0+exp(-(CurrentState[z].v - 14.48)/18.34));
	    double GKb = 0.003;
             if (iso==1)
	    {
		  GKb*=3*c_GKb_iso;//2.5;
	    }
             else GKb*=c_gkb;
	    if (celltype==1)
	    {
		GKb*=0.6;
	    }
	    IKb = GKb * xkb * (CurrentState[z].v - EK);
	    
	    double PNab=3.75e-10;
	    INab = PNab * vffrt * (CurrentState[z].nai * exp(vfrt)-nao)/(exp(vfrt)-1.0);
	    
	    double PCab=2.5e-8;
	    ICab=PCab*4.0*vffrt*(CurrentState[z].cai * exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
	    
	    double GpCa=0.0005;
	    IpCa = GpCa * CurrentState[z].cai/(0.0005 + CurrentState[z].cai);
        
        Inab_stop = clock();
        fprintf(foutput,"Inab: %e\n",(double)(Inab_stop - Inab_start) / (CLOCKS_PER_SEC));
        
        /* 15. Calculate IKatp*/
        #if IKATP
            double sigma = 3.8; // channels/mkm^2
            double prob  = 0.91; // open probability in the absence of ATP
        
            double gamma0 = 35.375*pow(ko/5.4, 0.24);
        
            double khMg0 = 0.65/sqrt(ko+5.0);
            double deltaMg = 0.32;
            double Mg_i = 0.5; // mM
            double khMg = khMg0 * exp(-2 * deltaMg * F * CurrentState[z].v / (R*T));
            double f_M = 1/(1 + Mg_i/khMg);
        
            double khNa0 = 25.9; // mmol/L
            double deltaNa = 0.35;
            double khNa = khNa0 * (-2 * deltaNa * F * CurrentState[z].v / (R*T));
            double f_N = 1/(1 + pow(nao/khNa, 2));
        
            double Q10 = 1.3; // temperature coefficient
            double T0 = 36 + 273.15; // K
            double f_T = pow(Q10, (T-T0)/10.0);
        
            double GKatp = gamma0 * f_M * f_N * f_T;
        
            double ADP_i = 0.015; //mmol
            double ATP_i = 5.0; //mmol
            double P_conc = 0.03; //mmol
            double Hill = 1.925;
        
            double kd1 = 1.292e-8;
            double kd2 = 1.328e-1;
            double kd3 = 8.0e3;
            double kd4 = 4.448e-1;
            double koef11 = 1.871e-5;
            double koef1 = 7.533e-9;
            double koef13 = 2.506e-6;
        
            double f_atp = 1 - pow(ATP_i, Hill)/(pow(ATP_i, Hill) + (koef1*koef11*koef13 + ADP_i*koef11*koef13 + koef1*koef13*P_conc + koef1*P_conc*ADP_i)/(koef1*koef11*koef13/kd1 + ADP_i*koef11*koef13/kd2 + koef1*koef13*P_conc/kd3 + koef1*P_conc*ADP_i/kd4));
            IKatp = sigma * GKatp * prob * f_atp * (CurrentState[z].v - EK);
        #endif

	    //----------old: stimulus()---------//
	
	 //   if ((t>(start+n*CL) && t<(start+duration+n*CL-dt)))
             if ((t>(start+n_stim*CL) && t<(start+duration+n_stim*CL-dt)))
    	    {
        	if (Ist==0)
        	{
            		vrest=CurrentState[z].v;
        	}
        	Ist=amp;
    	    }
    	 //   else if (t>(start+duration+n*CL-dt))
             else if (t>(start+duration+n_stim*CL-dt))
    	    {
        	Ist = 0.0;
        //	n = n+1;
               n_stim = n_stim+1;
    	    }
    	    vo = CurrentState[z].v;

	    //----------old: voltage()---------//

	    double I_gap_junc = 0;
	
	    #if CHAIN
            #if IKATP
                if(z<(chain_length-1))	I_gap_junc += -g_gap_junc*(CurrentState[z+1].v-CurrentState[z].v);
                if(z>0)			I_gap_junc += -g_gap_junc*(-CurrentState[z].v + CurrentState[z-1].v);
                if (z==0) NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist+IKatp+I_gap_junc);
                else NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+IKatp+I_gap_junc);
            #else
                if(z<(chain_length-1))	I_gap_junc += -g_gap_junc*(CurrentState[z+1].v-CurrentState[z].v);
                if(z>0)			I_gap_junc += -g_gap_junc*(-CurrentState[z].v + CurrentState[z-1].v);
                if (z==0) NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist+I_gap_junc);
                else NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+I_gap_junc);
            #endif
	    #else
            #if IKATP
                NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+IKatp+Ist);
            #else
                NextState[z].v = CurrentState[z].v-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist);
            #endif
        #endif
        
        
        #if TABLE
        {
            int ii; double dd;
            dd = (CurrentState[z].v - ETMIN)/(ETMAX - ETMIN)*(NTABLE - 1);
            //printf("v: %e\n", NextState[z].v);
            //printf("dd: %e\n", dd);
            ii = (int)dd;
            //printf("ii: %d\n", ii);
            
            T1.n1 = ii;
            if( T1.n1<0 || T1.n1>=NTABLE ) fprintf(stderr, "!!!T1.n1=%d\n", T1.n1);
            T1.frac = dd-ii;
            //printf("T1.frac: %e\n", T1.frac);
        
        }
        #endif //TABLE
        
        //----------old: dVdt_APD()---------//
	   
	    vdot_old = vdot;
	    vdot = (CurrentState[z].v - vo)/dt;
	    if (APD_flag == 0 && CurrentState[z].v > -40 && vdot < vdot_old)
	    {
		vdot_max = vdot_old;
		t_vdot_max = t - dt;
		APD_flag = 1;
	    }
	    if	(APD_flag==1 && CurrentState[z].v < 0.9 * vrest)
	    {
		APD = t - t_vdot_max;
		APD_flag = 0;
	    }

        //------------old: FBS------------//
	    
	    double CaMKb=CaMKo*(1.0 - CurrentState[z].CaMKt)/(1.0 + KmCaM/CurrentState[z].cass);
	    CaMKa = CaMKb + CurrentState[z].CaMKt;
	    NextState[z].CaMKt = CurrentState[z].CaMKt+dt * (aCaMK * CaMKb * (CaMKb + CurrentState[z].CaMKt) - bCaMK * CurrentState[z].CaMKt);
	    
	    JdiffNa = (CurrentState[z].nass - CurrentState[z].nai)/2.0;
	    JdiffK=(CurrentState[z].kss - CurrentState[z].ki)/2.0;
	    Jdiff = (CurrentState[z].cass - CurrentState[z].cai)/0.2;
            
            
	    
	    double bt=4.75;
	    double a_rel=0.5*bt;
            if (iso==1)
	    {
		   a_rel*=1.75*c_a_rel;//1.75;
	    }
	    double Jrel_inf=c_jrelinf*a_rel*(-ICaL)/(1.0+pow(1.5/CurrentState[z].cajsr,8.0));
	    if (celltype==2)
	    {
		Jrel_inf*=1.7;
	    }
	    double tau_rel=bt/(1.0+0.0123/CurrentState[z].cajsr);
             if (iso==1)
	    {
		 tau_rel*=0.5*c_tau_rel;//0.5;
	    }
            
	    if (tau_rel<0.005)
	    {
		tau_rel=0.005;
	    }
            
	    NextState[z].Jrelnp = Jrel_inf - (Jrel_inf - CurrentState[z].Jrelnp) * exp(-dt/tau_rel);
            
	    double btp=1.25*bt;
	    double a_relp=0.5*btp;
             if (iso==1)
	    {
		 a_relp*=1.75*c_a_rel;//1.75;
	    }
	    double Jrel_infp=c_jrelinf*a_relp*(-ICaL)/(1.0+pow(1.5/CurrentState[z].cajsr,8.0));
	    if (celltype==2)
	    {
		Jrel_infp*=1.7;
	    }
	    double tau_relp=btp/(1.0+0.0123/CurrentState[z].cajsr);
             if (iso==1)
	    {
		 tau_relp*=0.5*c_tau_rel;//0.5;
	    }

	    if (tau_relp<0.005)
	    {
		tau_relp=0.005;
	    }
           
            
	    NextState[z].Jrelp = Jrel_infp - (Jrel_infp - CurrentState[z].Jrelp) * exp(-dt/tau_relp);
	    double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
	    Jrel = (1.0-fJrelp) * CurrentState[z].Jrelnp + fJrelp * CurrentState[z].Jrelp;
	    
	    double Jupnp; 
	    double Jupp;
             if (iso==1)
	    {
		 Jupnp = c_jup * 0.004375 * CurrentState[z].cai/(CurrentState[z].cai + 0.00092*0.54*c_SERCA);
	    }
            else   Jupnp =c_jup* 0.004375 * CurrentState[z].cai/(CurrentState[z].cai + 0.00092);

	  
        if (iso==1)
	    {
		 Jupp = c_jup * 2.75 * 0.004375 * CurrentState[z].cai/(CurrentState[z].cai + (0.00092 - 0.00017)*0.54*c_SERCA);
	    }
            else   Jupp = c_jup *2.75 * 0.004375 * CurrentState[z].cai/(CurrentState[z].cai + 0.00092 - 0.00017);
	    if (celltype==1)
	    {
		Jupnp*=1.3;
		Jupp*=1.3;
	    }
	    double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
	    Jleak = 0.0039375 * CurrentState[z].cansr/15.0;
	    Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
	    
	    Jtr = (CurrentState[z].cansr - CurrentState[z].cajsr)/100.0;
	    
	    NextState[z].nai = CurrentState[z].nai +dt * (-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
	    NextState[z].nass = CurrentState[z].nass +dt * (-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);
	    
	    NextState[z].ki = CurrentState[z].ki+ dt * (-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
	    NextState[z].kss = CurrentState[z].kss+ dt * (-(ICaK) * Acap/(F * vss) - JdiffK);
	    
	    double Bcai;
	  //  if (celltype==1)

	 //   {
		Bcai=1.0/(1.0+c_CMDN*cmdnmax*kmcmdn/pow(kmcmdn + CurrentState[z].cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn + CurrentState[z].cai,2.0));
	  //  }
	  //  else
	  //  {
	//	Bcai=1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn + CurrentState[z].cai,2.0)+trpnmax*kmtrpn/pow(kmtrpn + CurrentState[z].cai,2.0));
	 //   }
	    NextState[z].cai = CurrentState[z].cai+ dt * (Bcai * (-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));
	    
	    double Bcass=1.0/(1.0+BSRmax*KmBSR/pow(KmBSR + CurrentState[z].cass,2.0)+BSLmax*KmBSL/pow(KmBSL + CurrentState[z].cass,2.0));
	    NextState[z].cass = CurrentState[z].cass+dt * (Bcass * (-(ICaL - 2.0 * INaCa_ss) * Acap/(2.0 * F * vss) + Jrel * vjsr/vss - Jdiff));
	    
	    NextState[z].cansr = CurrentState[z].cansr+dt * (Jup - Jtr * vjsr/vnsr);
	    
	    double Bcajsr = 1.0/(1.0 + csqnmax * kmcsqn/pow(kmcsqn + CurrentState[z].cajsr,2.0));
	    NextState[z].cajsr = CurrentState[z].cajsr+ dt * (Bcajsr * (Jtr - Jrel));
        
        ical_test_stop = clock();
        fprintf(foutput,"Time: %e\n",(double)(ical_test_stop - ical_test_start) / (CLOCKS_PER_SEC));
        fclose(foutput);
    }
    
}


void open(){
    t_= fopen("tests_ctrl/t.txt","w");
    v_= fopen("tests_ctrl/v.txt","w");
    nai_= fopen("tests_ctrl/nai.txt","w");
    nass_= fopen("tests_ctrl/nass.txt","w");
    ki_= fopen("tests_ctrl/ki.txt","w");
    kss_= fopen("tests_ctrl/kss.txt","w");
    cai_= fopen("tests_ctrl/cai.txt","w");
    cass_= fopen("tests_ctrl/cass.txt","w");
    cansr_= fopen("tests_ctrl/cansr.txt","w");
    cajsr_= fopen("tests_ctrl/cajsr.txt","w");
    Jrel_= fopen("tests_ctrl/Jrel.txt","w");
    CaMKt_= fopen("tests_ctrl/CaMKt.txt","w");
    Jup_= fopen("tests_ctrl/Jup.txt","w");
    Jtr_= fopen("tests_ctrl/Jtr.txt","w");
    Jdiff_= fopen("tests_ctrl/Jdiff.txt","w");
    JdiffNa_= fopen("tests_ctrl/JdiffNa.txt","w");
    JdiffK_= fopen("tests_ctrl/JdiffK.txt","w");
    Jleak_= fopen("tests_ctrl/Jleak.txt","w");
    INa_= fopen("tests_ctrl/INa.txt","w");
    INaL_= fopen("tests_ctrl/INaL.txt","w");
    Ito_= fopen("tests_ctrl/Ito.txt","w");
    ICaL_= fopen("tests_ctrl/ICaL.txt","w");
    ICaNa_= fopen("tests_ctrl/ICaNa.txt","w");
    ICaK_= fopen("tests_ctrl/ICaK.txt","w");
    IKr_= fopen("tests_ctrl/IKr.txt","w");
    IKs_= fopen("tests_ctrl/IKs.txt","w");
    IK1_= fopen("tests_ctrl/IK1.txt","w");
    INaCa_i_= fopen("tests_ctrl/INaCa_i.txt","w");
    INaCa_ss_= fopen("tests_ctrl/INaCa_ss.txt","w");
    INaCa_= fopen("tests_ctrl/INaCa.txt","w");
    INaK_= fopen("tests_ctrl/INaK.txt","w");
    IKb_= fopen("tests_ctrl/IKb.txt","w");
    INab_=fopen("tests_ctrl/INab.txt","w");
    IpCa_=fopen("tests_ctrl/IpCa.txt","w");
    ICab_= fopen("tests_ctrl/ICab.txt","w");
    Ist_= fopen("tests_ctrl/Ist.txt","w");
    dt_= fopen("tests_ctrl/dt.txt","w");
    APD_= fopen("tests_ctrl/APD.txt","w");
    ENa_= fopen("tests_ctrl/ENa.txt","w");
    EK_= fopen("tests_ctrl/EK.txt","w");
    EKs_= fopen("tests_ctrl/EKs.txt","w");


// for ICaL
    f_= fopen("tests_ctrl/f.txt","w");
    d_=fopen("tests_ctrl/d.txt","w");
    nca_=fopen("tests_ctrl/nca.txt","w");
    fca_= fopen("tests_ctrl/fca.txt","w");
    jca_= fopen("tests_ctrl/jca.txt","w");
    fp_= fopen("tests_ctrl/fp.txt","w");
    fcap_= fopen("tests_ctrl/fcap.txt","w");
    tff_= fopen("tests_ctrl/tff.txt","w");
    tfs_= fopen("tests_ctrl/tfs.txt","w");
    tffp_= fopen("tests_ctrl/tffp.txt","w");
    dss_= fopen("tests_ctrl/dss.txt","w");
    fss_= fopen("tests_ctrl/fss.txt","w");
    Ical_V=fopen("tests_ctrl/ICaL_V.txt","w");
//
}

void open_iso(){
    t_= fopen("tests/t.txt","w");
    v_= fopen("tests/v.txt","w");
    nai_= fopen("tests/nai.txt","w");
    nass_= fopen("tests/nass.txt","w");
    ki_= fopen("tests/ki.txt","w");
    kss_= fopen("tests/kss.txt","w");
    cai_= fopen("tests/cai.txt","w");
    cass_= fopen("tests/cass.txt","w");
    cansr_= fopen("tests/cansr.txt","w");
    cajsr_= fopen("tests/cajsr.txt","w");
    Jrel_= fopen("tests/Jrel.txt","w");
    CaMKt_= fopen("tests/CaMKt.txt","w");
    Jup_= fopen("tests/Jup.txt","w");
    Jtr_= fopen("tests/Jtr.txt","w");
    Jdiff_= fopen("tests/Jdiff.txt","w");
    JdiffNa_= fopen("tests/JdiffNa.txt","w");
    JdiffK_= fopen("tests/JdiffK.txt","w");
    Jleak_= fopen("tests/Jleak.txt","w");
    INa_= fopen("tests/INa.txt","w");
    INaL_= fopen("tests/INaL.txt","w");
    Ito_= fopen("tests/Ito.txt","w");
    ICaL_= fopen("tests/ICaL.txt","w");
    ICaNa_= fopen("tests/ICaNa.txt","w");
    ICaK_= fopen("tests/ICaK.txt","w");
    IKr_= fopen("tests/IKr.txt","w");
    IKs_= fopen("tests/IKs.txt","w");
    IK1_= fopen("tests/IK1.txt","w");
    INaCa_i_= fopen("tests/INaCa_i.txt","w");
    INaCa_ss_= fopen("tests/INaCa_ss.txt","w");
    INaCa_= fopen("tests/INaCa.txt","w");
    INaK_= fopen("tests/INaK.txt","w");
    IKb_= fopen("tests/IKb.txt","w");
    INab_=fopen("tests/INab.txt","w");
    IpCa_=fopen("tests/IpCa.txt","w");
    ICab_= fopen("tests/ICab.txt","w");
    Ist_= fopen("tests/Ist.txt","w");
    dt_= fopen("tests/dt.txt","w");
    APD_= fopen("tests/APD.txt","w");
    ENa_= fopen("tests/ENa.txt","w");
    EK_= fopen("tests/EK.txt","w");
    EKs_= fopen("tests/EKs.txt","w");

// for ICaL
    f_= fopen("tests/f.txt","w");
    d_=fopen("tests/d.txt","w");
    nca_=fopen("tests/nca.txt","w");
    fca_= fopen("tests/fca.txt","w");
    jca_= fopen("tests/jca.txt","w");
    fp_= fopen("tests/fp.txt","w");
    fcap_= fopen("tests/fcap.txt","w");
    tff_= fopen("tests/tff.txt","w");
    tfs_= fopen("tests/tds.txt","w");
    tffp_= fopen("tests/tffp.txt","w");
    dss_= fopen("tests/dss.txt","w");
    fss_= fopen("tests/fss.txt","w");
    Ical_V=fopen("tests/ICaL_V.txt","w");
//

}


void close(){
    fclose(t_);
    fclose(v_);
    fclose(nai_);
    fclose(nass_);
    fclose(ki_);
    fclose(kss_);
    fclose(cai_);
    fclose(cass_);
    fclose(cansr_);
    fclose(cajsr_);
    fclose(Jrel_);
    fclose(CaMKt_);
    fclose(Jup_);
    fclose(Jtr_);
    fclose(Jdiff_);
    fclose(JdiffNa_);
    fclose(JdiffK_);
    fclose(Jleak_);
    fclose(INa_);
    fclose(INaL_);
    fclose(Ito_);
    fclose(ICaL_);
    fclose(ICaNa_);
    fclose(ICaK_);
    fclose(IKr_);
    fclose(IKs_);
    fclose(IK1_);
    fclose(INaCa_i_);
    fclose(INaCa_ss_);
    fclose(INaCa_);
    fclose(INaK_);
    fclose(IKb_);
    fclose(INab_);
    fclose(IpCa_);
    fclose(ICab_);
    fclose(Ist_);
    fclose(dt_);
    fclose(APD_);
    fclose(ENa_);
    fclose(EK_);
    fclose(EKs_);

// for ICaL

    fclose(f_);
    fclose(d_);
    fclose(nca_);
    fclose(fca_);
    fclose(jca_);
    fclose(fp_);
    fclose(fcap_);
    fclose(tff_);
    fclose(tfs_);
    fclose(tffp_);
    fclose(dss_);
    fclose(fss_);
    fclose(Ical_V);

//
    
}

//initial values for state variables, there are 41 of them
void stateInitialization(struct State *CurrentState){
    CurrentState->v         = -87.5;
    CurrentState->nai       = 7.;
    CurrentState->nass      = 7.;
    CurrentState->ki        = 145.;
    CurrentState->kss       = 145.;
    CurrentState->cai       = 1.0e-4;
    CurrentState->cass      = 1.0e-4;
    CurrentState->cansr     = 1.2;
    CurrentState->cajsr     = 1.2;
    CurrentState->m         = 0;
    CurrentState->hf        = 1;
    CurrentState->hs        = 1;
    CurrentState->j         = 1;
    CurrentState->hsp       = 1;
    CurrentState->jp        = 1;
    CurrentState->mL        = 0;
    CurrentState->hL        = 1;
    CurrentState->hLp       = 1;
    CurrentState->a         = 0;
    CurrentState->iF        = 1;
    CurrentState->iS        = 1;
    CurrentState->ap        = 0;
    CurrentState->iFp       = 1;
    CurrentState->iSp       = 1;
    CurrentState->d         = 0;
    CurrentState->ff        = 1;
    CurrentState->fs        = 1;
    CurrentState->fcaf      = 1;
    CurrentState->fcas      = 1;
    CurrentState->jca       = 1;
    CurrentState->nca       = 0;
    CurrentState->ffp       = 1;
    CurrentState->fcafp     = 1;
    CurrentState->xrf       = 0;
    CurrentState->xrs       = 0;
    CurrentState->xs1       = 0;
    CurrentState->xs2       = 0;
    CurrentState->xk1       = 1;
    CurrentState->Jrelnp    = 0;
    CurrentState->Jrelp     = 0;
    CurrentState->CaMKt     = 0;
}

void initialization(struct State *CurrentState)
{
    #if READIS
        FILE *fin = fopen( "state_1000_IKR_from_ORd.dat", "r" );
        if(!fin) printf("!!! Cannot open IC file\n");
        fread( CurrentState, sizeof(struct State), 1, fin );
        fclose(fin);
    #else
        stateInitialization(CurrentState);
    #endif

}

float dv_max(struct State *CurrentState, struct State *NextState, int chain_length)

{

    int z;

    float max=0;    

    for (z=0; z<chain_length; z++)

    {

	if(fabs(NextState[z].v-CurrentState[z].v)>max) max=fabs(NextState[z].v-CurrentState[z].v);

    }

    return max;

}

//value holders for state varaibles in the case that the increase in dt was too aggressive, so a smaller one can be taken
double nai0,nass0,ki0,kss0,cai0,cass0,cansr0,cajsr0,m0,hf0,hs0,jO,hsp0,jp0,mL0,hL0,hLp0,a0,iF0,iS0,ap0,iFp0,iSp0,d0,ff0,fs0,fcaf0,fcas0,jca0,nca0,ffp0,fcafp0,xrf0,xrs0,xs10,xs20,xk10,Jrelnp0,Jrelp0,CaMKt0;

int main()
{
    clock_t time_start, time_end;
    time_start = clock();
    
    int chain_length = 1;
    int target_cell = 0;
    int z;
    float dv=0;
    float t_output=0;
    
    //struct State CurrentState;
    
    //memory allocation for the chain:
    if ((CurrentState = (struct State*)calloc( chain_length, sizeof(struct State))) == NULL) {printf ("not enough memory for CurrentState\n"); exit(-1);}
 
    if ((NextState = (struct State*)calloc( chain_length, sizeof(struct State))) == NULL) {printf ("not enough memory for CurrentState\n"); exit(-1);}
    
    for(z=0;z<chain_length;z++)
    {  
        initialization(&NextState[z]);
        CurrentState[z]=NextState[z];
    }
    
    if (iso==1) open_iso();
    else  open();
  
    
    
    while (t<=ft)
    {
            
       #if ADAPTIVE_STEP

/*	if (t>=ft-CL && t<ft-CL+400)

	{

		dt=0.005;

//			Count++; //increase the loop counter



	}

	else*/ if ((t>=(start+n_stim*CL-2)) && ((t<(start+duration+n_stim*CL))||(t<(start+duration+(n_stim-1)*CL+safetime)))) 

	{

		dt=0.005;

	}

	else if((dv>0.08))

	{



		t=t-dt;

		dt=fabs(0.02*dt/dv);

//		for (z=0; z<chain_length; z++)

//		{

//		    CurrentState[z]=NextState[z];

//		}

		RGC(CurrentState,NextState, z, chain_length);

		t=t+dt;

		while((dv_max(CurrentState,NextState, chain_length)>0.08)&&(dt>0.005))

		{

			t=t-dt;

			dt=dt/2.0;

			if(dt<0.005) dt=0.005;	

			RGC(CurrentState,NextState, z, chain_length);

			t=t+dt;

		} 

	}

	else if (dv<0.02)

	{

		dt=fabs(0.08*dt/dv);

		if (dt>0.09)	

		{

			dt=0.09;

		}

	}
#endif

	for (z=0; z<chain_length; z++)

	{

	    CurrentState[z] = NextState[z];

	}

	RGC(CurrentState,NextState, z, chain_length);

	t=t+dt;

            

	dv=dv_max(CurrentState,NextState, chain_length);
        //if (Count%skip==0 && (t>=ft-beatssave*CL) && (t-(ft-beatssave*CL))<=400)//save results ot output file when the sampling interval and time are correct
     //   if ((t>=t_output) && (t>90000))
    
    #if ADAPTIVE_STEP
         if (t>=t_output)
            {   t_output+=0.5;//0.5
                if(t>400000)	printf("%f\n",t);
    #else 
         if (Count%skip==0 && (t>=ft-beatssave*CL) && (t-(ft-beatssave*CL))<=400)//save results ot output file when the sampling interval and time are correct
                {
    #endif
                
                fprintf(v_,"%e\t",t);
                fprintf(nai_,"%e\t",t);
                fprintf(nass_,"%e\t",t);
                fprintf(ki_,"%e\t",t);
                fprintf(kss_,"%e\t",t);
                fprintf(cai_,"%e\t",t);
                fprintf(cass_,"%e\t",t);
                fprintf(cansr_,"%e\t",t);
                fprintf(cajsr_,"%e\t",t);
                fprintf(Jrel_,"%e\t",t);
                fprintf(CaMKt_,"%e\t",t);
                fprintf(Jup_,"%e\t",t);
                fprintf(Jtr_,"%e\t",t);
                fprintf(Jdiff_,"%e\t",t);
                fprintf(JdiffNa_,"%e\t",t);
                fprintf(JdiffK_,"%e\t",t);
                fprintf(Jleak_,"%e\t",t);
                fprintf(INa_,"%e\t",t);
                fprintf(INaL_,"%e\t",t);
                fprintf(Ito_,"%e\t",t);
                fprintf(ICaL_,"%e\t",t);
                fprintf(ICaNa_,"%e\t",t);
                fprintf(ICaK_,"%e\t",t);
                fprintf(IKr_,"%e\t",t);
                fprintf(IKs_,"%e\t",t);
                fprintf(IK1_,"%e\t",t);
                fprintf(INaCa_i_,"%e\t",t);
                fprintf(INaCa_ss_,"%e\t",t);
                fprintf(INaCa_,"%e\t",t);
                fprintf(INaK_,"%e\t",t);
                fprintf(IKb_,"%e\t",t);
                fprintf(INab_,"%e\t",t);

                fprintf(IpCa_,"%e\t",t);
                fprintf(ICab_,"%e\t",t);
                fprintf(Ist_,"%e\t",t);
                fprintf(APD_,"%e\t",t);
                fprintf(ENa_,"%e\t",t);
                fprintf(EK_,"%e\t",t);
                fprintf(EKs_,"%e\t",t);


// for ICaL   
                fprintf(f_,"%e\t",t);
                fprintf(d_,"%e\t",t);
                fprintf(nca_,"%e\t",t);
                fprintf(fca_,"%e\t",t);
                fprintf(jca_,"%e\t",t);
                fprintf(fp_,"%e\t",t);
                fprintf(fcap_,"%e\t",t);
                fprintf(tff_,"%e\t",t);
                fprintf(tfs_,"%e\t",t);
                fprintf(tffp_,"%e\t",t);
                //fprintf(dss_,"%e\t",CurrentState[z].v);
               // fprintf(fss_,"%e\t",CurrentState[z].v);
              
                




		for (z = 0;z < chain_length;z++)
		{
		        fprintf(v_,"%e\t", CurrentState[z].v);
		        /*fprintf(nai_,"%e\t", CurrentState[z].nai);
		        fprintf(nass_,"%e\t", CurrentState[z].nass);
		        fprintf(ki_,"%e\t", CurrentState[z].ki);
		        fprintf(kss_,"%e\t", CurrentState[z].kss);
		        fprintf(cai_,"%e\t", CurrentState[z].cai);
		        fprintf(cass_,"%e\t", CurrentState[z].cass);
		        fprintf(cansr_,"%e\t", CurrentState[z].cansr);
		        fprintf(cajsr_,"%e\t", CurrentState[z].cajsr);
		        fprintf(Jrel_,"%e\t", Jrel);
		        fprintf(CaMKt_,"%e\t", CurrentState[z].CaMKt);
		        fprintf(Jup_,"%e\t", Jup);
		        fprintf(Jtr_,"%e\t", Jtr);
		        fprintf(Jdiff_,"%e\t", Jdiff);
		        fprintf(JdiffNa_,"%e\t", JdiffNa);
		        fprintf(JdiffK_,"%e\t", JdiffK);
		        fprintf(Jleak_,"%e\t", Jleak);
		        fprintf(INa_,"%e\t", INa);
		        fprintf(INaL_,"%e\t", INaL);
		        fprintf(Ito_,"%e\t", Ito);
		        fprintf(ICaL_,"%e\t", ICaL);
		        fprintf(ICaNa_,"%e\t", ICaNa);
		        fprintf(ICaK_,"%e\t", ICaK);
		        fprintf(IKr_,"%e\t", IKr);
		        fprintf(IKs_,"%e\t", IKs);
		        fprintf(IK1_,"%e\t", IK1);
		        fprintf(INaCa_i_,"%e\t", INaCa_i);
		        fprintf(INaCa_ss_,"%e\t", INaCa_ss);
		        fprintf(INaCa_,"%e\t",  INaCa);
		        fprintf(INaK_,"%e\t", INaK);
		        fprintf(IKb_,"%e\t", IKb);
		        fprintf(INab_,"%e\t", INab);
		        fprintf(IpCa_,"%e\t", IpCa);
		        fprintf(ICab_,"%e\t", ICab);
		        fprintf(Ist_,"%e\t",  Ist);
		        fprintf(dt_,"%e\n",dt);
		        fprintf(APD_,"%e\t", APD);
                        fprintf(ENa_,"%e\t",  ENa);
		        fprintf(EK_,"%e\t",EK);
		        fprintf(EKs_,"%e\t", EKs);
                         fprintf(dss_,"%e\t %e\t",CurrentState[z].v,dss);
                         fprintf(fss_,"%e\t %e\t",CurrentState[z].v,fss);
                        if(CurrentState[z].v!=-70) { fprintf(Ical_V,"%e\t %e\t",CurrentState[z].v,ICaL);}*/
// for ICaL   
                /*fprintf(f_,"%e\t",f);
                fprintf(d_,"%e\t",CurrentState[z].d);
                fprintf(nca_,"%e\t",CurrentState[z].nca);
                fprintf(fca_,"%e\t",fca);
                fprintf(jca_,"%e\t",CurrentState[z].jca);
                fprintf(fp_,"%e\t",fp);
                fprintf(fcap_,"%e\t",fcap);
                fprintf(tff_,"%e\t",tff);
                fprintf(tfs_,"%e\t",tfs);
                fprintf(tffp_,"%e\t",tffp);*/
//


		}
		
		fprintf(v_,"\n");
                /*fprintf(nai_,"\n");
                fprintf(nass_,"\n");
                fprintf(ki_,"\n");
                fprintf(kss_,"\n");
                fprintf(cai_,"\n");
                fprintf(cass_,"\n");
                fprintf(cansr_,"\n");
                fprintf(cajsr_,"\n");
                fprintf(Jrel_,"\n");
                fprintf(CaMKt_,"\n");
                fprintf(Jup_,"\n");
                fprintf(Jtr_,"\n");
                fprintf(Jdiff_,"\n");
                fprintf(JdiffNa_,"\n");
                fprintf(JdiffK_,"\n");
                fprintf(Jleak_,"\n");
                fprintf(INa_,"\n");
                fprintf(INaL_,"\n");
                fprintf(Ito_,"\n");
                fprintf(ICaL_,"\n");
                fprintf(ICaNa_,"\n");
                fprintf(ICaK_,"\n");
                fprintf(IKr_,"\n");
                fprintf(IKs_,"\n");
                fprintf(IK1_,"\n");
                fprintf(INaCa_i_,"\n");
                fprintf(INaCa_ss_,"\n");
                fprintf(INaCa_,"\n");
                fprintf(INaK_,"\n");
                fprintf(IKb_,"\n");
                fprintf(INab_,"\n");
                fprintf(IpCa_,"\n");
                fprintf(ICab_,"\n");
                fprintf(Ist_,"\n");
                fprintf(APD_,"\n");
                fprintf(ENa_,"\n");
                fprintf(EK_,"\n");
                fprintf(EKs_,"\n");
// for ICaL   
                fprintf(f_,"\n");
                fprintf(d_,"\n");
                fprintf(nca_,"\n");
                fprintf(fca_,"\n");
                fprintf(jca_,"\n");
                fprintf(fp_,"\n");
                fprintf(fcap_,"\n");
                fprintf(tff_,"\n");
                fprintf(tfs_,"\n");
                fprintf(tffp_,"\n");
                fprintf(dss_,"\n");
                fprintf(fss_,"\n");
                fprintf(Ical_V,"\n");*/
 
            }
  
      int tmp1 = 0;

#if WRITEIS

        if ((t > (np-1)*CL)&&(tmp1 == 0)){

            FILE *fout = fopen( "state_2000_iso_15.dat", "w" );

            fwrite(&NextState[14], sizeof(struct State), 1, fout);

            fclose(fout);

            tmp1 += 1;}

#endif
               Count++;//increase the loop counter

      
      
    }
    close();
        
    time_end = clock();
    printf("CPU_time: %lf\n",(double)(time_end - time_start) / (CLOCKS_PER_SEC));
return 0;
}
