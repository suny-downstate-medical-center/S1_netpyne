/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__StochKv_det
#define _nrn_initial _nrn_initial__StochKv_det
#define nrn_cur _nrn_cur__StochKv_det
#define _nrn_current _nrn_current__StochKv_det
#define nrn_jacob _nrn_jacob__StochKv_det
#define nrn_state _nrn_state__StochKv_det
#define _net_receive _net_receive__StochKv_det 
#define ChkProb ChkProb__StochKv_det 
#define _f_trates _f_trates__StochKv_det 
#define setRNG setRNG__StochKv_det 
#define states states__StochKv_det 
#define trates trates__StochKv_det 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gamma _p[0]
#define eta _p[1]
#define gkbar _p[2]
#define ik _p[3]
#define gk _p[4]
#define N _p[5]
#define N0 _p[6]
#define N1 _p[7]
#define n0_n1 _p[8]
#define n1_n0 _p[9]
#define n _p[10]
#define ek _p[11]
#define scale_dens _p[12]
#define n0_n1_new _p[13]
#define Dn _p[14]
#define v _p[15]
#define _g _p[16]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
#define rng	*_ppvar[3]._pval
#define _p_rng	_ppvar[3]._pval
#define area	*_ppvar[4]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  3;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_BnlDev(void);
 static void _hoc_ChkProb(void);
 static void _hoc_SigmoidRate(void);
 static void _hoc_brand(void);
 static void _hoc_setRNG(void);
 static void _hoc_strap(void);
 static void _hoc_trates(void);
 static void _hoc_urand(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_StochKv_det", _hoc_setdata,
 "BnlDev_StochKv_det", _hoc_BnlDev,
 "ChkProb_StochKv_det", _hoc_ChkProb,
 "SigmoidRate_StochKv_det", _hoc_SigmoidRate,
 "brand_StochKv_det", _hoc_brand,
 "setRNG_StochKv_det", _hoc_setRNG,
 "strap_StochKv_det", _hoc_strap,
 "trates_StochKv_det", _hoc_trates,
 "urand_StochKv_det", _hoc_urand,
 0, 0
};
#define BnlDev BnlDev_StochKv_det
#define SigmoidRate SigmoidRate_StochKv_det
#define brand brand_StochKv_det
#define strap strap_StochKv_det
#define urand urand_StochKv_det
 extern double BnlDev( _threadargsprotocomma_ double , double );
 extern double SigmoidRate( _threadargsprotocomma_ double , double , double , double );
 extern double brand( _threadargsprotocomma_ double , double );
 extern double strap( _threadargsprotocomma_ double );
 extern double urand( _threadargsproto_ );
 
static void _check_trates(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_trates(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[7];
#define _gth 0
#define P_b_StochKv_det _thread1data[0]
#define P_b _thread[_gth]._pval[0]
#define P_a_StochKv_det _thread1data[1]
#define P_a _thread[_gth]._pval[1]
#define Rb Rb_StochKv_det
 double Rb = 0.002;
#define Ra Ra_StochKv_det
 double Ra = 0.02;
#define a_StochKv_det _thread1data[2]
#define a _thread[_gth]._pval[2]
#define b_StochKv_det _thread1data[3]
#define b _thread[_gth]._pval[3]
#define deterministic deterministic_StochKv_det
 double deterministic = 1;
#define ntau_StochKv_det _thread1data[4]
#define ntau _thread[_gth]._pval[4]
#define ninf_StochKv_det _thread1data[5]
#define ninf _thread[_gth]._pval[5]
#define qa qa_StochKv_det
 double qa = 9;
#define q10 q10_StochKv_det
 double q10 = 2.3;
#define tha tha_StochKv_det
 double tha = -40;
#define tadj_StochKv_det _thread1data[6]
#define tadj _thread[_gth]._pval[6]
#define temp temp_StochKv_det
 double temp = 23;
#define usetable usetable_StochKv_det
 double usetable = 1;
#define vmax vmax_StochKv_det
 double vmax = 100;
#define vmin vmin_StochKv_det
 double vmin = -120;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_StochKv_det", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tha_StochKv_det", "mV",
 "Ra_StochKv_det", "/ms",
 "Rb_StochKv_det", "/ms",
 "temp_StochKv_det", "degC",
 "vmin_StochKv_det", "mV",
 "vmax_StochKv_det", "mV",
 "a_StochKv_det", "/ms",
 "b_StochKv_det", "/ms",
 "ntau_StochKv_det", "ms",
 "gamma_StochKv_det", "pS",
 "eta_StochKv_det", "1/um2",
 "gkbar_StochKv_det", "S/cm2",
 "ik_StochKv_det", "mA/cm2",
 "gk_StochKv_det", "S/cm2",
 0,0
};
 static double delta_t = 1;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tha_StochKv_det", &tha_StochKv_det,
 "qa_StochKv_det", &qa_StochKv_det,
 "Ra_StochKv_det", &Ra_StochKv_det,
 "Rb_StochKv_det", &Rb_StochKv_det,
 "temp_StochKv_det", &temp_StochKv_det,
 "q10_StochKv_det", &q10_StochKv_det,
 "deterministic_StochKv_det", &deterministic_StochKv_det,
 "vmin_StochKv_det", &vmin_StochKv_det,
 "vmax_StochKv_det", &vmax_StochKv_det,
 "a_StochKv_det", &a_StochKv_det,
 "b_StochKv_det", &b_StochKv_det,
 "ninf_StochKv_det", &ninf_StochKv_det,
 "ntau_StochKv_det", &ntau_StochKv_det,
 "tadj_StochKv_det", &tadj_StochKv_det,
 "P_a_StochKv_det", &P_a_StochKv_det,
 "P_b_StochKv_det", &P_b_StochKv_det,
 "usetable_StochKv_det", &usetable_StochKv_det,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"StochKv_det",
 "gamma_StochKv_det",
 "eta_StochKv_det",
 "gkbar_StochKv_det",
 0,
 "ik_StochKv_det",
 "gk_StochKv_det",
 "N_StochKv_det",
 "N0_StochKv_det",
 "N1_StochKv_det",
 "n0_n1_StochKv_det",
 "n1_n0_StochKv_det",
 0,
 "n_StochKv_det",
 0,
 "rng_StochKv_det",
 0};
 extern Node* nrn_alloc_node_;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	gamma = 30;
 	eta = 0;
 	gkbar = 0.75;
 	_prop->param = _p;
 	_prop->param_size = 17;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 	_ppvar[4]._pval = &nrn_alloc_node_->_area; /* diam */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 static void bbcore_write(double*, int*, int*, int*, _threadargsproto_);
 extern void hoc_reg_bbcore_write(int, void(*)(double*, int*, int*, int*, _threadargsproto_));
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _StochKv_det_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
   hoc_reg_bbcore_write(_mechtype, bbcore_write);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 17, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "bbcorepointer");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 4, "area");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 StochKv_det /home/fernando/S1_netpyne/sim/mod/StochKv_det.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_ntau;
 static double *_t_ninf;
 static double *_t_a;
 static double *_t_b;
 static double *_t_tadj;
static int _reset;
static char *modelname = "skm95.mod  ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int ChkProb(_threadargsprotocomma_ double);
static int _f_trates(_threadargsprotocomma_ double);
static int setRNG(_threadargsproto_);
static int trates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_trates(_threadargsprotocomma_ double _lv);
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*VERBATIM*/
#include "nrnran123.h"
extern int cvode_active_;
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   trates ( _threadargscomma_ v ) ;
   Dn = a - ( a + b ) * n ;
   if ( deterministic  || dt > 1.0 ) {
     N1 = n * N ;
     }
   else {
     P_a = strap ( _threadargscomma_ a * dt ) ;
     P_b = strap ( _threadargscomma_ b * dt ) ;
     ChkProb ( _threadargscomma_ P_a ) ;
     ChkProb ( _threadargscomma_ P_b ) ;
     n0_n1 = BnlDev ( _threadargscomma_ P_a , N0 ) ;
     n1_n0 = BnlDev ( _threadargscomma_ P_b , N1 ) ;
     N0 = strap ( _threadargscomma_ N0 - n0_n1 + n1_n0 ) ;
     N1 = N - N0 ;
     }
   N0 = N - N1 ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 trates ( _threadargscomma_ v ) ;
 Dn = Dn  / (1. - dt*( ( - ( ( a + b ) )*( 1.0 ) ) )) ;
 if ( deterministic  || dt > 1.0 ) {
   N1 = n * N ;
   }
 else {
   P_a = strap ( _threadargscomma_ a * dt ) ;
   P_b = strap ( _threadargscomma_ b * dt ) ;
   ChkProb ( _threadargscomma_ P_a ) ;
   ChkProb ( _threadargscomma_ P_b ) ;
   n0_n1 = BnlDev ( _threadargscomma_ P_a , N0 ) ;
   n1_n0 = BnlDev ( _threadargscomma_ P_b , N1 ) ;
   N0 = strap ( _threadargscomma_ N0 - n0_n1 + n1_n0 ) ;
   N1 = N - N0 ;
   }
 N0 = N - N1 ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   trates ( _threadargscomma_ v ) ;
    n = n + (1. - exp(dt*(( - ( ( a + b ) )*( 1.0 ) ))))*(- ( a ) / ( ( - ( ( a + b ) )*( 1.0 ) ) ) - n) ;
   if ( deterministic  || dt > 1.0 ) {
     N1 = n * N ;
     }
   else {
     P_a = strap ( _threadargscomma_ a * dt ) ;
     P_b = strap ( _threadargscomma_ b * dt ) ;
     ChkProb ( _threadargscomma_ P_a ) ;
     ChkProb ( _threadargscomma_ P_b ) ;
     n0_n1 = BnlDev ( _threadargscomma_ P_a , N0 ) ;
     n1_n0 = BnlDev ( _threadargscomma_ P_b , N1 ) ;
     N0 = strap ( _threadargscomma_ N0 - n0_n1 + n1_n0 ) ;
     N1 = N - N0 ;
     }
   N0 = N - N1 ;
   }
  return 0;
}
 static double _mfac_trates, _tmin_trates;
  static void _check_trates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_Ra;
  static double _sav_Rb;
  static double _sav_tha;
  static double _sav_qa;
  static double _sav_q10;
  static double _sav_temp;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_Ra != Ra) { _maktable = 1;}
  if (_sav_Rb != Rb) { _maktable = 1;}
  if (_sav_tha != tha) { _maktable = 1;}
  if (_sav_qa != qa) { _maktable = 1;}
  if (_sav_q10 != q10) { _maktable = 1;}
  if (_sav_temp != temp) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_trates)/199.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 200; _x += _dx, _i++) {
    _f_trates(_p, _ppvar, _thread, _nt, _x);
    _t_ntau[_i] = ntau;
    _t_ninf[_i] = ninf;
    _t_a[_i] = a;
    _t_b[_i] = b;
    _t_tadj[_i] = tadj;
   }
   _sav_dt = dt;
   _sav_Ra = Ra;
   _sav_Rb = Rb;
   _sav_tha = tha;
   _sav_qa = qa;
   _sav_q10 = q10;
   _sav_temp = temp;
   _sav_celsius = celsius;
  }
 }

 static int trates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_trates(_p, _ppvar, _thread, _nt);
#endif
 _n_trates(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_trates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  ntau = _xi;
  ninf = _xi;
  a = _xi;
  b = _xi;
  tadj = _xi;
  return;
 }
 if (_xi <= 0.) {
 ntau = _t_ntau[0];
 ninf = _t_ninf[0];
 a = _t_a[0];
 b = _t_b[0];
 tadj = _t_tadj[0];
 return; }
 if (_xi >= 199.) {
 ntau = _t_ntau[199];
 ninf = _t_ninf[199];
 a = _t_a[199];
 b = _t_b[199];
 tadj = _t_tadj[199];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 a = _t_a[_i] + _theta*(_t_a[_i+1] - _t_a[_i]);
 b = _t_b[_i] + _theta*(_t_b[_i+1] - _t_b[_i]);
 tadj = _t_tadj[_i] + _theta*(_t_tadj[_i+1] - _t_tadj[_i]);
 }

 
static int  _f_trates ( _threadargsprotocomma_ double _lv ) {
   tadj = pow( q10 , ( ( celsius - temp ) / ( 10.0 ) ) ) ;
   a = SigmoidRate ( _threadargscomma_ _lv , tha , Ra , qa ) ;
   a = a * tadj ;
   b = SigmoidRate ( _threadargscomma_ - _lv , - tha , Rb , qa ) ;
   b = b * tadj ;
   ntau = 1.0 / ( a + b ) ;
   ninf = a * ntau ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_trates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 trates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double SigmoidRate ( _threadargsprotocomma_ double _lv , double _lth , double _la , double _lq ) {
   double _lSigmoidRate;
  if ( fabs ( _lv - _lth ) > 1e-6 ) {
     _lSigmoidRate = _la * ( _lv - _lth ) / ( 1.0 - exp ( - ( _lv - _lth ) / _lq ) ) ;
      }
   else {
     _lSigmoidRate = _la * _lq ;
     }
   
return _lSigmoidRate;
 }
 
static void _hoc_SigmoidRate(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  SigmoidRate ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double strap ( _threadargsprotocomma_ double _lx ) {
   double _lstrap;
 if ( _lx < 0.0 ) {
     _lstrap = 0.0 ;
     
/*VERBATIM*/
        fprintf (stderr,"skv.mod:strap: negative state");
 }
   else {
     _lstrap = _lx ;
     }
   
return _lstrap;
 }
 
static void _hoc_strap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  strap ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  ChkProb ( _threadargsprotocomma_ double _lp ) {
   if ( _lp < 0.0  || _lp > 1.0 ) {
     
/*VERBATIM*/
    fprintf(stderr, "StochKv.mod:ChkProb: argument not a probability.\n");
 }
    return 0; }
 
static void _hoc_ChkProb(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 ChkProb ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  setRNG ( _threadargsproto_ ) {
   
/*VERBATIM*/
    {
#if !NRNBBCORE
    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
    uint32_t a2 = 0;
    if (*pv) {
        nrnran123_deletestream(*pv);
        *pv = (nrnran123_State*)0;
    } 
    if (ifarg(2)) {
        a2 = (uint32_t)*getarg(2);
    }
    if (ifarg(1)) {
        *pv = nrnran123_newstream((uint32_t)*getarg(1), a2);
    }
#endif
    }
  return 0; }
 
static void _hoc_setRNG(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 setRNG ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
double urand ( _threadargsproto_ ) {
   double _lurand;
 
/*VERBATIM*/
    double value;
    if (_p_rng) {
        /*
        :Supports separate independent but reproducible streams for
        : each instance.
            : distribution MUST be set to Random.uniform(0,1)
        */
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
        //printf("random stream for this simulation = %lf\n",value);
        return value;
    }else{
        //assert(0);
        value = 0.0;
    }
    _lurand = value;
 
return _lurand;
 }
 
static void _hoc_urand(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  urand ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
/*VERBATIM*/
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    if (d) {
        uint32_t* di = ((uint32_t*)d) + *offset;
      // temporary just enough to see how much space is being used
      if (!_p_rng) {
        di[0] = 0; di[1] = 0;
      }else{
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        nrnran123_getids(*pv, di, di+1);
      }
//printf("StochKv.mod %p: bbcore_write offset=%d %d %d\n", _p, *offset, d?di[0]:-1, d?di[1]:-1);
    }
    *offset += 2;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    assert(!_p_rng);
    uint32_t* di = ((uint32_t*)d) + *offset;
        if (di[0] != 0 || di[1] != 0)
        {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
      *pv = nrnran123_newstream(di[0], di[1]);
        }
//printf("StochKv.mod %p: bbcore_read offset=%d %d %d\n", _p, *offset, di[0], di[1]);
    *offset += 2;
}
 
double brand ( _threadargsprotocomma_ double _lP , double _lN ) {
   double _lbrand;
 
/*VERBATIM*/
        /*
        :Supports separate independent but reproducible streams for
        : each instance. However, the corresponding hoc Random
        : distribution MUST be set to Random.uniform(0,1)
        */

        // Should probably be optimized
        double value = 0.0;
        int i;
        for (i = 0; i < _lN; i++) {
           if (urand(_threadargs_) < _lP) {
              value = value + 1;
           }
        }
        return(value);

 _lbrand = value ;
   
return _lbrand;
 }
 
static void _hoc_brand(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  brand ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
/*VERBATIM*/
#define        PI 3.141592654
#define        r_ia     16807
#define        r_im     2147483647
#define        r_am     (1.0/r_im)
#define        r_iq     127773
#define        r_ir     2836
#define        r_ntab   32
#define        r_ndiv   (1+(r_im-1)/r_ntab)
#define        r_eps    1.2e-7
#define        r_rnmx   (1.0-r_eps)
 
/*VERBATIM*/
/* ---------------------------------------------------------------- */
/* gammln - compute natural log of gamma function of xx */
static double
gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}
 
double BnlDev ( _threadargsprotocomma_ double _lppr , double _lnnr ) {
   double _lBnlDev;
 
/*VERBATIM*/
        int j;
        double am,em,g,angle,p,bnl,sq,bt,y;
        double pc,plog,pclog,en,oldg;

        /* prepare to always ignore errors within this routine */

        p=(_lppr <= 0.5 ? _lppr : 1.0-_lppr);
        am=_lnnr*p;
        if (_lnnr < 25) {
            bnl=0.0;
            for (j=1;j<=_lnnr;j++)
                if (urand(_threadargs_) < p) bnl += 1.0;
        }
        else if (am < 1.0) {
            g=exp(-am);
            bt=1.0;
            for (j=0;j<=_lnnr;j++) {
                bt *= urand(_threadargs_);
                if (bt < g) break;
            }
            bnl=(j <= _lnnr ? j : _lnnr);
        }
        else {
            {
                en=_lnnr;
                oldg=gammln(en+1.0);
            }
            {
                pc=1.0-p;
                plog=log(p);
                pclog=log(pc);
            }
            sq=sqrt(2.0*am*pc);
            do {
                do {
                    angle=PI*urand(_threadargs_);
                        angle=PI*urand(_threadargs_);
                    y=tan(angle);
                    em=sq*y+am;
                } while (em < 0.0 || em >= (en+1.0));
                em=floor(em);
                    bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) - 
                    gammln(en-em+1.0)+em*plog+(en-em)*pclog);
            } while (urand(_threadargs_) > bt);
            bnl=em;
        }
        if (p != _lppr) bnl=_lnnr-bnl;

        /* recover error if changed during this routine, thus ignoring
            any errors during this routine */


        return bnl;

 _lBnlDev = bnl ;
   
return _lBnlDev;
 }
 
static void _hoc_BnlDev(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  BnlDev ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_mem_init(Datum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(7, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  n = n0;
 {
   
/*VERBATIM*/
    if (_p_rng) {
        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
    }
    if (cvode_active_ && !deterministic) {
        hoc_execerror("StochKv with deterministic=0", "cannot be used with cvode");
    }
 eta = gkbar / gamma ;
   trates ( _threadargscomma_ v ) ;
   n = ninf ;
   scale_dens = gamma / area ;
   N = floor ( eta * area + 0.5 ) ;
   N1 = n * N ;
   if (  ! deterministic ) {
     N1 = floor ( N1 + 0.5 ) ;
     }
   N0 = N - N1 ;
   n0_n1 = 0.0 ;
   n1_n0 = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_trates(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gk = ( strap ( _threadargscomma_ N1 ) * scale_dens * tadj ) ;
   ik = 1e-4 * gk * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
   _t_ntau = makevector(200*sizeof(double));
   _t_ninf = makevector(200*sizeof(double));
   _t_a = makevector(200*sizeof(double));
   _t_b = makevector(200*sizeof(double));
   _t_tadj = makevector(200*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/fernando/S1_netpyne/sim/mod/StochKv_det.mod";
static const char* nmodl_file_text = 
  "TITLE skm95.mod  \n"
  "\n"
  "COMMENT\n"
  "----------------------------------------------------------------\n"
  "Stochastic version of the K channel mechanism kd3h5.mod by\n"
  "Z. Mainen in Mainen & Sejnowski 95.\n"
  "\n"
  "This represents a potassium channel, with Hodgkin-Huxley like kinetics,\n"
  "based on the gates model, assuming stochastic opening and closing.\n"
  "\n"
  "Kinetic rates based roughly on Sah et al. and Hamill et al. (1991)\n"
  "The main kinetic difference from the standard H-H model (shh.mod) is \n"
  "that the K+ kinetic is different, not n^4, but just n, \n"
  "and the activation curves are different.\n"
  "\n"
  "The rate functions are adapted directly from the Kd3h5.mod file\n"
  "by Zach Mainen.\n"
  "\n"
  "The stochastic model is as following:\n"
  "\n"
  "Potassium\n"
  "\n"
  "       = alpha_n =>      \n"
  "   [N0]             [N1]\n"
  "      <= beta_n =      \n"
  "\n"
  "\n"
  "The model keeps track on the number of channels in each state, and \n"
  "uses a binomial distribution to update these number.\n"
  "\n"
  "Jan 1999, Mickey London, Hebrew University, mikilon@lobster.ls.huji.ac.il\n"
  "        Peter N. Steinmetz, Caltech, peter@klab.caltech.edu\n"
  "14 Sep 99 PNS. Added deterministic flag.\n"
  "19 May 2002 Kamran Diba.  Changed gamma and deterministic from GLOBAL to RANGE.\n"
  "23 Nov 2011 Werner Van Geit @ BBP. Changed the file so that it can use the neuron random number generator. Tuned voltage dependence\n"
  "----------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX StochKv_det\n"
  "    THREADSAFE\n"
  "    USEION k READ ek WRITE ik\n"
  "    RANGE N,eta, gk, gamma, gkbar, ik, N0, N1, n0_n1, n1_n0\n"
  "    GLOBAL ninf, ntau,a,b,P_a,P_b, deterministic\n"
  "    GLOBAL Ra, Rb\n"
  "    GLOBAL vmin, vmax, q10, temp, tadj\n"
  "    BBCOREPOINTER rng\n"
  "    :POINTER rng\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)\n"
  "    (pS) = (picosiemens)\n"
  "    (S) = (siemens)\n"
  "    (um) = (micron)\n"
  "} \n"
  "\n"
  "PARAMETER {\n"
  "    v           (mV)\n"
  "    dt      (ms)\n"
  "    area    (um2)\n"
  "\n"
  "    gamma  =  30          (pS)\n"
  "    eta              (1/um2)\n"
  "    gkbar = .75      (S/cm2)\n"
  "\n"
  "    tha  = -40   (mV)        : v 1/2 for inf\n"
  "    qa   = 9            : inf slope     \n"
  "    Ra   = 0.02 (/ms)       : max act rate\n"
  "    Rb   = 0.002    (/ms)       : max deact rate\n"
  "\n"
  "    celsius (degC)\n"
  "    temp = 23 (degC)   : original temperature for kinetic set\n"
  "    q10 = 2.3               : temperature sensitivity\n"
  "\n"
  "    deterministic = 1   : if non-zero, will use deterministic version\n"
  "    vmin = -120 (mV)    : range to construct tables for\n"
  "    vmax = 100  (mV)\n"
  "} \n"
  "\n"
  "ASSIGNED {\n"
  "    a       (/ms)\n"
  "    b       (/ms)\n"
  "    ik      (mA/cm2)\n"
  "    gk      (S/cm2)\n"
  "    ek      (mV)\n"
  "    ninf        : steady-state value\n"
  "    ntau (ms)   : time constant for relaxation\n"
  "    tadj\n"
  "\n"
  "    N \n"
  "    scale_dens (pS/um2) \n"
  "    P_a     : probability of one channel making alpha transition\n"
  "    P_b     : probability of one channel making beta transition\n"
  "\n"
  "    rng\n"
  "\n"
  "    n0_n1_new\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "    n   \n"
  "}\n"
  "ASSIGNED {\n"
  "    N0 N1       : N states populations\n"
  "    n0_n1 n1_n0 : number of channels moving from one state to the other \n"
  "}\n"
  "\n"
  "COMMENT\n"
  "The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1 \n"
  "for comparison with Pr to decide whether to activate the synapse or not\n"
  "ENDCOMMENT\n"
  "\n"
  "VERBATIM\n"
  "#include \"nrnran123.h\"\n"
  "extern int cvode_active_;\n"
  "ENDVERBATIM\n"
  ": ----------------------------------------------------------------\n"
  ": initialization\n"
  "INITIAL { \n"
  "\n"
  "    VERBATIM\n"
  "    if (_p_rng) {\n"
  "        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);\n"
  "    }\n"
  "    if (cvode_active_ && !deterministic) {\n"
  "        hoc_execerror(\"StochKv with deterministic=0\", \"cannot be used with cvode\");\n"
  "    }\n"
  "    ENDVERBATIM\n"
  "\n"
  "    eta = gkbar / gamma\n"
  "    trates(v)\n"
  "    n = ninf\n"
  "    scale_dens = gamma/area\n"
  "    N = floor(eta*area + 0.5)\n"
  "\n"
  "    N1 = n*N\n"
  "  if (!deterministic) {\n"
  "    N1 = floor(N1 + 0.5)\n"
  "  }\n"
  "    N0 = N-N1       : any round off into non-conducting state\n"
  "\n"
  "    n0_n1 = 0\n"
  "    n1_n0 = 0\n"
  "}\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": Breakpoint for each integration step\n"
  "BREAKPOINT {\n"
  "  SOLVE states METHOD cnexp\n"
  "\n"
  "  gk =  (strap(N1) * scale_dens * tadj)\n"
  "\n"
  "  ik = 1e-4 * gk * (v - ek)\n"
  "} \n"
  "\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": states - updates number of channels in each state\n"
  "DERIVATIVE states {\n"
  "\n"
  "    trates(v)\n"
  "\n"
  "        n' = a - (a + b)*n\n"
  "    if (deterministic || dt > 1) { : ForwardSkip is also deterministic\n"
  "    N1 = n*N\n"
  "    }else{\n"
  "\n"
  "    P_a = strap(a*dt)\n"
  "    P_b = strap(b*dt)\n"
  "\n"
  "    : check that will represent probabilities when used\n"
  "    ChkProb( P_a)\n"
  "    ChkProb( P_b)\n"
  "\n"
  "    : transitions\n"
  "    n0_n1 = BnlDev(P_a, N0)\n"
  "    n1_n0 = BnlDev(P_b, N1)\n"
  "\n"
  "    : move the channels\n"
  "    N0    = strap(N0 - n0_n1 + n1_n0)\n"
  "    N1    = N - N0\n"
  "\n"
  "    }\n"
  "\n"
  "    N0 = N-N1       : any round off into non-conducting state\n"
  "}\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": trates - compute rates, using table if possible\n"
  "PROCEDURE trates(v (mV)) {     \n"
  "    TABLE ntau, ninf, a, b, tadj\n"
  "    DEPEND dt, Ra, Rb, tha, qa, q10, temp, celsius\n"
  "    FROM vmin TO vmax WITH 199\n"
  "\n"
  "    tadj = q10 ^ ((celsius - temp)/(10 (K)))\n"
  "    a = SigmoidRate(v, tha, Ra, qa)\n"
  "    a = a * tadj\n"
  "    b = SigmoidRate(-v, -tha, Rb, qa)\n"
  "    b = b * tadj\n"
  "    ntau = 1/(a+b)\n"
  "    ninf = a*ntau\n"
  "}\n"
  "\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": SigmoidRate - Compute a sigmoid rate function given the \n"
  ": 50% point th, the slope q, and the amplitude a.\n"
  "FUNCTION SigmoidRate(v (mV),th (mV),a (1/ms),q) (1/ms){\n"
  "    UNITSOFF\n"
  "    if (fabs(v-th) > 1e-6 ) {\n"
  "        SigmoidRate = a * (v - th) / (1 - exp(-(v - th)/q))\n"
  "    UNITSON\n"
  "\n"
  "    } else {\n"
  "        SigmoidRate = a * q\n"
  "    }\n"
  "}   \n"
  "\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": sign trap - trap for negative values and replace with zero\n"
  "FUNCTION strap(x) {\n"
  "    if (x < 0) {\n"
  "        strap = 0\n"
  "VERBATIM\n"
  "        fprintf (stderr,\"skv.mod:strap: negative state\");\n"
  "ENDVERBATIM\n"
  "    } else {\n"
  "        strap = x\n"
  "    }\n"
  "}\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": ChkProb - Check that number represents a probability\n"
  "PROCEDURE ChkProb(p) {\n"
  "\n"
  "  if (p < 0.0 || p > 1.0) {\n"
  "    VERBATIM\n"
  "    fprintf(stderr, \"StochKv.mod:ChkProb: argument not a probability.\\n\");\n"
  "    ENDVERBATIM\n"
  "  }\n"
  "\n"
  "}\n"
  "\n"
  "PROCEDURE setRNG() {\n"
  "VERBATIM\n"
  "    {\n"
  "#if !NRNBBCORE\n"
  "    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "    uint32_t a2 = 0;\n"
  "    if (*pv) {\n"
  "        nrnran123_deletestream(*pv);\n"
  "        *pv = (nrnran123_State*)0;\n"
  "    } \n"
  "    if (ifarg(2)) {\n"
  "        a2 = (uint32_t)*getarg(2);\n"
  "    }\n"
  "    if (ifarg(1)) {\n"
  "        *pv = nrnran123_newstream((uint32_t)*getarg(1), a2);\n"
  "    }\n"
  "#endif\n"
  "    }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION urand() {\n"
  "VERBATIM\n"
  "    double value;\n"
  "    if (_p_rng) {\n"
  "        /*\n"
  "        :Supports separate independent but reproducible streams for\n"
  "        : each instance.\n"
  "            : distribution MUST be set to Random.uniform(0,1)\n"
  "        */\n"
  "        value = nrnran123_dblpick((nrnran123_State*)_p_rng);\n"
  "        //printf(\"random stream for this simulation = %lf\\n\",value);\n"
  "        return value;\n"
  "    }else{\n"
  "        //assert(0);\n"
  "        value = 0.0;\n"
  "    }\n"
  "    _lurand = value;\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {\n"
  "    if (d) {\n"
  "        uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "      // temporary just enough to see how much space is being used\n"
  "      if (!_p_rng) {\n"
  "        di[0] = 0; di[1] = 0;\n"
  "      }else{\n"
  "        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "        nrnran123_getids(*pv, di, di+1);\n"
  "      }\n"
  "//printf(\"StochKv.mod %p: bbcore_write offset=%d %d %d\\n\", _p, *offset, d?di[0]:-1, d?di[1]:-1);\n"
  "    }\n"
  "    *offset += 2;\n"
  "}\n"
  "static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {\n"
  "    assert(!_p_rng);\n"
  "    uint32_t* di = ((uint32_t*)d) + *offset;\n"
  "        if (di[0] != 0 || di[1] != 0)\n"
  "        {\n"
  "      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);\n"
  "      *pv = nrnran123_newstream(di[0], di[1]);\n"
  "        }\n"
  "//printf(\"StochKv.mod %p: bbcore_read offset=%d %d %d\\n\", _p, *offset, di[0], di[1]);\n"
  "    *offset += 2;\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  ": Returns random numbers drawn from a binomial distribution\n"
  "FUNCTION brand(P, N) {\n"
  "\n"
  "VERBATIM\n"
  "        /*\n"
  "        :Supports separate independent but reproducible streams for\n"
  "        : each instance. However, the corresponding hoc Random\n"
  "        : distribution MUST be set to Random.uniform(0,1)\n"
  "        */\n"
  "\n"
  "        // Should probably be optimized\n"
  "        double value = 0.0;\n"
  "        int i;\n"
  "        for (i = 0; i < _lN; i++) {\n"
  "           if (urand(_threadargs_) < _lP) {\n"
  "              value = value + 1;\n"
  "           }\n"
  "        }\n"
  "        return(value);\n"
  "\n"
  "ENDVERBATIM\n"
  "\n"
  "        brand = value\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "#define        PI 3.141592654\n"
  "#define        r_ia     16807\n"
  "#define        r_im     2147483647\n"
  "#define        r_am     (1.0/r_im)\n"
  "#define        r_iq     127773\n"
  "#define        r_ir     2836\n"
  "#define        r_ntab   32\n"
  "#define        r_ndiv   (1+(r_im-1)/r_ntab)\n"
  "#define        r_eps    1.2e-7\n"
  "#define        r_rnmx   (1.0-r_eps)\n"
  "ENDVERBATIM\n"
  "\n"
  "VERBATIM\n"
  "/* ---------------------------------------------------------------- */\n"
  "/* gammln - compute natural log of gamma function of xx */\n"
  "static double\n"
  "gammln(double xx)\n"
  "{\n"
  "    double x,tmp,ser;\n"
  "    static double cof[6]={76.18009173,-86.50532033,24.01409822,\n"
  "        -1.231739516,0.120858003e-2,-0.536382e-5};\n"
  "    int j;\n"
  "    x=xx-1.0;\n"
  "    tmp=x+5.5;\n"
  "    tmp -= (x+0.5)*log(tmp);\n"
  "    ser=1.0;\n"
  "    for (j=0;j<=5;j++) {\n"
  "        x += 1.0;\n"
  "        ser += cof[j]/x;\n"
  "    }\n"
  "    return -tmp+log(2.50662827465*ser);\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  "\n"
  ": ----------------------------------------------------------------\n"
  ": BnlDev - draw a uniform deviate from the generator\n"
  "FUNCTION BnlDev (ppr, nnr) {\n"
  "\n"
  "VERBATIM\n"
  "        int j;\n"
  "        double am,em,g,angle,p,bnl,sq,bt,y;\n"
  "        double pc,plog,pclog,en,oldg;\n"
  "\n"
  "        /* prepare to always ignore errors within this routine */\n"
  "\n"
  "        p=(_lppr <= 0.5 ? _lppr : 1.0-_lppr);\n"
  "        am=_lnnr*p;\n"
  "        if (_lnnr < 25) {\n"
  "            bnl=0.0;\n"
  "            for (j=1;j<=_lnnr;j++)\n"
  "                if (urand(_threadargs_) < p) bnl += 1.0;\n"
  "        }\n"
  "        else if (am < 1.0) {\n"
  "            g=exp(-am);\n"
  "            bt=1.0;\n"
  "            for (j=0;j<=_lnnr;j++) {\n"
  "                bt *= urand(_threadargs_);\n"
  "                if (bt < g) break;\n"
  "            }\n"
  "            bnl=(j <= _lnnr ? j : _lnnr);\n"
  "        }\n"
  "        else {\n"
  "            {\n"
  "                en=_lnnr;\n"
  "                oldg=gammln(en+1.0);\n"
  "            }\n"
  "            {\n"
  "                pc=1.0-p;\n"
  "                plog=log(p);\n"
  "                pclog=log(pc);\n"
  "            }\n"
  "            sq=sqrt(2.0*am*pc);\n"
  "            do {\n"
  "                do {\n"
  "                    angle=PI*urand(_threadargs_);\n"
  "                        angle=PI*urand(_threadargs_);\n"
  "                    y=tan(angle);\n"
  "                    em=sq*y+am;\n"
  "                } while (em < 0.0 || em >= (en+1.0));\n"
  "                em=floor(em);\n"
  "                    bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) - \n"
  "                    gammln(en-em+1.0)+em*plog+(en-em)*pclog);\n"
  "            } while (urand(_threadargs_) > bt);\n"
  "            bnl=em;\n"
  "        }\n"
  "        if (p != _lppr) bnl=_lnnr-bnl;\n"
  "\n"
  "        /* recover error if changed during this routine, thus ignoring\n"
  "            any errors during this routine */\n"
  "\n"
  "\n"
  "        return bnl;\n"
  "\n"
  "    ENDVERBATIM\n"
  "    BnlDev = bnl\n"
  "}  \n"
  "\n"
  ;
#endif
