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
 
#define nrn_init _nrn_init__StochKv_deterministic
#define _nrn_initial _nrn_initial__StochKv_deterministic
#define nrn_cur _nrn_cur__StochKv_deterministic
#define _nrn_current _nrn_current__StochKv_deterministic
#define nrn_jacob _nrn_jacob__StochKv_deterministic
#define nrn_state _nrn_state__StochKv_deterministic
#define _net_receive _net_receive__StochKv_deterministic 
#define rates rates__StochKv_deterministic 
#define states states__StochKv_deterministic 
 
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
#define gmax _p[0]
#define conductance _p[1]
#define q10ConductanceScaling_q10Factor _p[2]
#define q10ConductanceScaling_experimentalTemp _p[3]
#define n_instances _p[4]
#define n_reverseRate_rate _p[5]
#define n_reverseRate_midpoint _p[6]
#define n_reverseRate_scale _p[7]
#define n_forwardRate_rate _p[8]
#define n_forwardRate_midpoint _p[9]
#define n_forwardRate_scale _p[10]
#define n_q10Settings_q10Factor _p[11]
#define n_q10Settings_experimentalTemp _p[12]
#define n_q10Settings_TENDEGREES _p[13]
#define gion _p[14]
#define q10ConductanceScaling_factor _p[15]
#define n_reverseRate_x _p[16]
#define n_reverseRate_r _p[17]
#define n_forwardRate_x _p[18]
#define n_forwardRate_r _p[19]
#define n_q10Settings_q10 _p[20]
#define n_rateScale _p[21]
#define n_alpha _p[22]
#define n_beta _p[23]
#define n_fcond _p[24]
#define n_inf _p[25]
#define n_tau _p[26]
#define conductanceScale _p[27]
#define fopenHHrates _p[28]
#define fopenHHtauInf _p[29]
#define fopenHHratesTau _p[30]
#define fopenHHratesInf _p[31]
#define fopenHHratesTauInf _p[32]
#define fopen _p[33]
#define g _p[34]
#define n_q _p[35]
#define temperature _p[36]
#define ek _p[37]
#define ik _p[38]
#define rate_n_q _p[39]
#define Dn_q _p[40]
#define v _p[41]
#define _g _p[42]
#define _ion_ik	*_ppvar[0]._pval
#define _ion_dikdv	*_ppvar[1]._pval
 
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
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
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
 "setdata_StochKv_deterministic", _hoc_setdata,
 "rates_StochKv_deterministic", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_StochKv_deterministic", "S/cm2",
 "conductance_StochKv_deterministic", "uS",
 "q10ConductanceScaling_experimentalTemp_StochKv_deterministic", "K",
 "n_reverseRate_rate_StochKv_deterministic", "kHz",
 "n_reverseRate_midpoint_StochKv_deterministic", "mV",
 "n_reverseRate_scale_StochKv_deterministic", "mV",
 "n_forwardRate_rate_StochKv_deterministic", "kHz",
 "n_forwardRate_midpoint_StochKv_deterministic", "mV",
 "n_forwardRate_scale_StochKv_deterministic", "mV",
 "n_q10Settings_experimentalTemp_StochKv_deterministic", "K",
 "n_q10Settings_TENDEGREES_StochKv_deterministic", "K",
 "gion_StochKv_deterministic", "S/cm2",
 "n_reverseRate_r_StochKv_deterministic", "kHz",
 "n_forwardRate_r_StochKv_deterministic", "kHz",
 "n_alpha_StochKv_deterministic", "kHz",
 "n_beta_StochKv_deterministic", "kHz",
 "n_tau_StochKv_deterministic", "ms",
 "g_StochKv_deterministic", "uS",
 0,0
};
 static double delta_t = 0.01;
 static double n_q0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"StochKv_deterministic",
 "gmax_StochKv_deterministic",
 "conductance_StochKv_deterministic",
 "q10ConductanceScaling_q10Factor_StochKv_deterministic",
 "q10ConductanceScaling_experimentalTemp_StochKv_deterministic",
 "n_instances_StochKv_deterministic",
 "n_reverseRate_rate_StochKv_deterministic",
 "n_reverseRate_midpoint_StochKv_deterministic",
 "n_reverseRate_scale_StochKv_deterministic",
 "n_forwardRate_rate_StochKv_deterministic",
 "n_forwardRate_midpoint_StochKv_deterministic",
 "n_forwardRate_scale_StochKv_deterministic",
 "n_q10Settings_q10Factor_StochKv_deterministic",
 "n_q10Settings_experimentalTemp_StochKv_deterministic",
 "n_q10Settings_TENDEGREES_StochKv_deterministic",
 0,
 "gion_StochKv_deterministic",
 "q10ConductanceScaling_factor_StochKv_deterministic",
 "n_reverseRate_x_StochKv_deterministic",
 "n_reverseRate_r_StochKv_deterministic",
 "n_forwardRate_x_StochKv_deterministic",
 "n_forwardRate_r_StochKv_deterministic",
 "n_q10Settings_q10_StochKv_deterministic",
 "n_rateScale_StochKv_deterministic",
 "n_alpha_StochKv_deterministic",
 "n_beta_StochKv_deterministic",
 "n_fcond_StochKv_deterministic",
 "n_inf_StochKv_deterministic",
 "n_tau_StochKv_deterministic",
 "conductanceScale_StochKv_deterministic",
 "fopenHHrates_StochKv_deterministic",
 "fopenHHtauInf_StochKv_deterministic",
 "fopenHHratesTau_StochKv_deterministic",
 "fopenHHratesInf_StochKv_deterministic",
 "fopenHHratesTauInf_StochKv_deterministic",
 "fopen_StochKv_deterministic",
 "g_StochKv_deterministic",
 0,
 "n_q_StochKv_deterministic",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 43, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	conductance = 1e-05;
 	q10ConductanceScaling_q10Factor = 2.3;
 	q10ConductanceScaling_experimentalTemp = 296.15;
 	n_instances = 1;
 	n_reverseRate_rate = 0.018;
 	n_reverseRate_midpoint = -40;
 	n_reverseRate_scale = -9;
 	n_forwardRate_rate = 0.18;
 	n_forwardRate_midpoint = -40;
 	n_forwardRate_scale = 9;
 	n_q10Settings_q10Factor = 2.3;
 	n_q10Settings_experimentalTemp = 296.15;
 	n_q10Settings_TENDEGREES = 10;
 	_prop->param = _p;
 	_prop->param_size = 43;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _StochKv_deterministic_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", 1.0);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 43, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 StochKv_deterministic /home/fernando/S1_netpyne/sim/x86_64/StochKv_deterministic.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Mod file for component: Component(id=StochKv_deterministic type=ionChannelHH)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dn_q = rate_n_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargs_ ) ;
 Dn_q = Dn_q  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargs_ ) ;
    n_q = n_q - dt*(- ( rate_n_q ) ) ;
   }
  return 0;
}
 
static int  rates ( _threadargsproto_ ) {
   q10ConductanceScaling_factor = pow( q10ConductanceScaling_q10Factor , ( ( temperature - q10ConductanceScaling_experimentalTemp ) / 10.0 ) ) ;
   n_reverseRate_x = ( v - n_reverseRate_midpoint ) / n_reverseRate_scale ;
   if ( n_reverseRate_x  != 0.0 ) {
     n_reverseRate_r = n_reverseRate_rate * n_reverseRate_x / ( 1.0 - exp ( 0.0 - n_reverseRate_x ) ) ;
     }
   else if ( n_reverseRate_x  == 0.0 ) {
     n_reverseRate_r = n_reverseRate_rate ;
     }
   n_forwardRate_x = ( v - n_forwardRate_midpoint ) / n_forwardRate_scale ;
   if ( n_forwardRate_x  != 0.0 ) {
     n_forwardRate_r = n_forwardRate_rate * n_forwardRate_x / ( 1.0 - exp ( 0.0 - n_forwardRate_x ) ) ;
     }
   else if ( n_forwardRate_x  == 0.0 ) {
     n_forwardRate_r = n_forwardRate_rate ;
     }
   n_q10Settings_q10 = pow( n_q10Settings_q10Factor , ( ( temperature - n_q10Settings_experimentalTemp ) / n_q10Settings_TENDEGREES ) ) ;
   n_rateScale = n_q10Settings_q10 ;
   n_alpha = n_forwardRate_r ;
   n_beta = n_reverseRate_r ;
   n_fcond = pow( n_q , n_instances ) ;
   n_inf = n_alpha / ( n_alpha + n_beta ) ;
   n_tau = 1.0 / ( ( n_alpha + n_beta ) * n_rateScale ) ;
   rate_n_q = ( n_inf - n_q ) / n_tau ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt );
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
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  n_q = n_q0;
 {
   ek = - 85.0 ;
   temperature = celsius + 273.15 ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   n_q = n_inf ;
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
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   conductanceScale = q10ConductanceScaling_factor ;
   fopenHHrates = n_fcond ;
   fopenHHtauInf = 1.0 ;
   fopenHHratesTau = 1.0 ;
   fopenHHratesInf = 1.0 ;
   fopenHHratesTauInf = 1.0 ;
   fopen = conductanceScale * fopenHHrates * fopenHHtauInf * fopenHHratesTau * fopenHHratesInf * fopenHHratesTauInf ;
   g = conductance * fopen ;
   gion = gmax * fopen ;
   ik = gion * ( v - ek ) ;
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
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n_q) - _p;  _dlist1[0] = &(Dn_q) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/fernando/S1_netpyne/sim/mod/StochKv_deterministic.mod";
static const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=StochKv_deterministic type=ionChannelHH)\n"
  "\n"
  "COMMENT\n"
  "\n"
  "    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)\n"
  "         org.neuroml.export  v1.4.2\n"
  "         org.neuroml.model   v1.4.2\n"
  "         jLEMS               v0.9.7.3\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX StochKv_deterministic\n"
  "    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!\n"
  "    \n"
  "    RANGE gion                           \n"
  "    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!\n"
  "    RANGE conductance                       : parameter\n"
  "    \n"
  "    RANGE g                                 : exposure\n"
  "    \n"
  "    RANGE fopen                             : exposure\n"
  "    RANGE q10ConductanceScaling_q10Factor   : parameter\n"
  "    RANGE q10ConductanceScaling_experimentalTemp: parameter\n"
  "    \n"
  "    RANGE q10ConductanceScaling_factor      : exposure\n"
  "    RANGE n_instances                       : parameter\n"
  "    \n"
  "    RANGE n_alpha                           : exposure\n"
  "    \n"
  "    RANGE n_beta                            : exposure\n"
  "    \n"
  "    RANGE n_tau                             : exposure\n"
  "    \n"
  "    RANGE n_inf                             : exposure\n"
  "    \n"
  "    RANGE n_rateScale                       : exposure\n"
  "    \n"
  "    RANGE n_fcond                           : exposure\n"
  "    RANGE n_reverseRate_rate                : parameter\n"
  "    RANGE n_reverseRate_midpoint            : parameter\n"
  "    RANGE n_reverseRate_scale               : parameter\n"
  "    \n"
  "    RANGE n_reverseRate_r                   : exposure\n"
  "    RANGE n_forwardRate_rate                : parameter\n"
  "    RANGE n_forwardRate_midpoint            : parameter\n"
  "    RANGE n_forwardRate_scale               : parameter\n"
  "    \n"
  "    RANGE n_forwardRate_r                   : exposure\n"
  "    RANGE n_q10Settings_q10Factor           : parameter\n"
  "    RANGE n_q10Settings_experimentalTemp    : parameter\n"
  "    RANGE n_q10Settings_TENDEGREES          : parameter\n"
  "    \n"
  "    RANGE n_q10Settings_q10                 : exposure\n"
  "    RANGE n_reverseRate_x                   : derived variable\n"
  "    RANGE n_forwardRate_x                   : derived variable\n"
  "    RANGE conductanceScale                  : derived variable\n"
  "    RANGE fopenHHrates                      : derived variable\n"
  "    RANGE fopenHHtauInf                     : derived variable\n"
  "    RANGE fopenHHratesTau                   : derived variable\n"
  "    RANGE fopenHHratesInf                   : derived variable\n"
  "    RANGE fopenHHratesTauInf                : derived variable\n"
  "    \n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    \n"
  "    (nA) = (nanoamp)\n"
  "    (uA) = (microamp)\n"
  "    (mA) = (milliamp)\n"
  "    (A) = (amp)\n"
  "    (mV) = (millivolt)\n"
  "    (mS) = (millisiemens)\n"
  "    (uS) = (microsiemens)\n"
  "    (molar) = (1/liter)\n"
  "    (kHz) = (kilohertz)\n"
  "    (mM) = (millimolar)\n"
  "    (um) = (micrometer)\n"
  "    (umol) = (micromole)\n"
  "    (S) = (siemens)\n"
  "    \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    \n"
  "    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!\n"
  "    \n"
  "    conductance = 1.0E-5 (uS)\n"
  "    q10ConductanceScaling_q10Factor = 2.3 \n"
  "    q10ConductanceScaling_experimentalTemp = 296.15 (K)\n"
  "    n_instances = 1 \n"
  "    n_reverseRate_rate = 0.018000001 (kHz)\n"
  "    n_reverseRate_midpoint = -40 (mV)\n"
  "    n_reverseRate_scale = -9 (mV)\n"
  "    n_forwardRate_rate = 0.18 (kHz)\n"
  "    n_forwardRate_midpoint = -40 (mV)\n"
  "    n_forwardRate_scale = 9 (mV)\n"
  "    n_q10Settings_q10Factor = 2.3 \n"
  "    n_q10Settings_experimentalTemp = 296.15 (K)\n"
  "    n_q10Settings_TENDEGREES = 10 (K)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    \n"
  "    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel\n"
  "    v (mV)\n"
  "    celsius (degC)\n"
  "    temperature (K)\n"
  "    ek (mV)\n"
  "    ik (mA/cm2)\n"
  "    \n"
  "    \n"
  "    q10ConductanceScaling_factor           : derived variable\n"
  "    \n"
  "    n_reverseRate_x                        : derived variable\n"
  "    \n"
  "    n_reverseRate_r (kHz)                  : conditional derived var...\n"
  "    \n"
  "    n_forwardRate_x                        : derived variable\n"
  "    \n"
  "    n_forwardRate_r (kHz)                  : conditional derived var...\n"
  "    \n"
  "    n_q10Settings_q10                      : derived variable\n"
  "    \n"
  "    n_rateScale                            : derived variable\n"
  "    \n"
  "    n_alpha (kHz)                          : derived variable\n"
  "    \n"
  "    n_beta (kHz)                           : derived variable\n"
  "    \n"
  "    n_fcond                                : derived variable\n"
  "    \n"
  "    n_inf                                  : derived variable\n"
  "    \n"
  "    n_tau (ms)                             : derived variable\n"
  "    \n"
  "    conductanceScale                       : derived variable\n"
  "    \n"
  "    fopenHHrates                           : derived variable\n"
  "    \n"
  "    fopenHHtauInf                          : derived variable\n"
  "    \n"
  "    fopenHHratesTau                        : derived variable\n"
  "    \n"
  "    fopenHHratesInf                        : derived variable\n"
  "    \n"
  "    fopenHHratesTauInf                     : derived variable\n"
  "    \n"
  "    fopen                                  : derived variable\n"
  "    \n"
  "    g (uS)                                 : derived variable\n"
  "    rate_n_q (/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    n_q \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    ek = -85.0\n"
  "    \n"
  "    temperature = celsius + 273.15\n"
  "    \n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    n_q = n_inf\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=StochKv_deterministic type=ionChannelHH), from conductanceScaling; Component(id=null type=q10ConductanceScaling)\n"
  "    conductanceScale = q10ConductanceScaling_factor ? multiply applied to all instances of factor in: <conductanceScaling> ([Component(id=null type=q10ConductanceScaling)]) c2 ([Component(id=null type=annotation), Component(id=null type=notes), Component(id=null type=q10ConductanceScaling), Component(id=n type=gateHHrates)]) ? path based\n"
  "    \n"
  "    ? DerivedVariable is based on path: gatesHHrates[*]/fcond, on: Component(id=StochKv_deterministic type=ionChannelHH), from gatesHHrates; Component(id=n type=gateHHrates)\n"
  "    fopenHHrates = n_fcond ? multiply applied to all instances of fcond in: <gatesHHrates> ([Component(id=n type=gateHHrates)]) c2 ([Component(id=null type=annotation), Component(id=null type=notes), Component(id=null type=q10ConductanceScaling), Component(id=n type=gateHHrates)]) ? path based\n"
  "    \n"
  "    ? DerivedVariable is based on path: gatesHHtauInf[*]/fcond, on: Component(id=StochKv_deterministic type=ionChannelHH), from gatesHHtauInf; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    fopenHHtauInf = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gatesHHratesTau[*]/fcond, on: Component(id=StochKv_deterministic type=ionChannelHH), from gatesHHratesTau; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    fopenHHratesTau = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gatesHHratesInf[*]/fcond, on: Component(id=StochKv_deterministic type=ionChannelHH), from gatesHHratesInf; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    fopenHHratesInf = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gatesHHratesTauInf[*]/fcond, on: Component(id=StochKv_deterministic type=ionChannelHH), from gatesHHratesTauInf; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    fopenHHratesTauInf = 1 \n"
  "    \n"
  "    fopen = conductanceScale  *  fopenHHrates  *  fopenHHtauInf  *  fopenHHratesTau  *  fopenHHratesInf  *  fopenHHratesTauInf ? evaluable\n"
  "    g = conductance  *  fopen ? evaluable\n"
  "    gion = gmax * fopen \n"
  "    \n"
  "    ik = gion * (v - ek)\n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    n_q' = rate_n_q \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    \n"
  "    q10ConductanceScaling_factor = q10ConductanceScaling_q10Factor ^((temperature -  q10ConductanceScaling_experimentalTemp )/10) ? evaluable\n"
  "    n_reverseRate_x = (v -  n_reverseRate_midpoint ) /  n_reverseRate_scale ? evaluable\n"
  "    if (n_reverseRate_x  != 0)  { \n"
  "        n_reverseRate_r = n_reverseRate_rate  *  n_reverseRate_x  / (1 - exp(0 -  n_reverseRate_x )) ? evaluable cdv\n"
  "    } else if (n_reverseRate_x  == 0)  { \n"
  "        n_reverseRate_r = n_reverseRate_rate ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    n_forwardRate_x = (v -  n_forwardRate_midpoint ) /  n_forwardRate_scale ? evaluable\n"
  "    if (n_forwardRate_x  != 0)  { \n"
  "        n_forwardRate_r = n_forwardRate_rate  *  n_forwardRate_x  / (1 - exp(0 -  n_forwardRate_x )) ? evaluable cdv\n"
  "    } else if (n_forwardRate_x  == 0)  { \n"
  "        n_forwardRate_r = n_forwardRate_rate ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    n_q10Settings_q10 = n_q10Settings_q10Factor ^((temperature -  n_q10Settings_experimentalTemp )/ n_q10Settings_TENDEGREES ) ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=n type=gateHHrates), from q10Settings; Component(id=null type=q10ExpTemp)\n"
  "    n_rateScale = n_q10Settings_q10 ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10ExpTemp)]) c2 ([Component(id=null type=HHExpLinearRate), Component(id=null type=HHExpLinearRate), Component(id=null type=q10ExpTemp)]) ? path based\n"
  "    \n"
  "    ? DerivedVariable is based on path: forwardRate/r, on: Component(id=n type=gateHHrates), from forwardRate; Component(id=null type=HHExpLinearRate)\n"
  "    n_alpha = n_forwardRate_r ? path based\n"
  "    \n"
  "    ? DerivedVariable is based on path: reverseRate/r, on: Component(id=n type=gateHHrates), from reverseRate; Component(id=null type=HHExpLinearRate)\n"
  "    n_beta = n_reverseRate_r ? path based\n"
  "    \n"
  "    n_fcond = n_q ^ n_instances ? evaluable\n"
  "    n_inf = n_alpha /( n_alpha + n_beta ) ? evaluable\n"
  "    n_tau = 1/(( n_alpha + n_beta ) *  n_rateScale ) ? evaluable\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    rate_n_q = ( n_inf  -  n_q ) /  n_tau ? Note units of all quantities used here need to be consistent!\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "}\n"
  "\n"
  ;
#endif
