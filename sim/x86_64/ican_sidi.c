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
 
#define nrn_init _nrn_init__ican
#define _nrn_initial _nrn_initial__ican
#define nrn_cur _nrn_cur__ican
#define _nrn_current _nrn_current__ican
#define nrn_jacob _nrn_jacob__ican
#define nrn_state _nrn_state__ican
#define _net_receive _net_receive__ican 
#define evaluate_fct evaluate_fct__ican 
#define iassign iassign__ican 
#define states states__ican 
 
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
#define gbar _p[0]
#define i _p[1]
#define m_inf _p[2]
#define tau_m _p[3]
#define m _p[4]
#define cai _p[5]
#define Dm _p[6]
#define ina _p[7]
#define tadj _p[8]
#define g _p[9]
#define v _p[10]
#define _g _p[11]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
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
 static void _hoc_evaluate_fct(void);
 static void _hoc_iassign(void);
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
 "setdata_ican", _hoc_setdata,
 "evaluate_fct_ican", _hoc_evaluate_fct,
 "iassign_ican", _hoc_iassign,
 0, 0
};
 /* declare global and static user variables */
#define beta beta_ican
 double beta = 0.0001;
#define cac cac_ican
 double cac = 0.0004;
#define erev erev_ican
 double erev = -20;
#define taumin taumin_ican
 double taumin = 0.1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "erev_ican", "mV",
 "beta_ican", "1/ms",
 "cac_ican", "mM",
 "taumin_ican", "ms",
 "gbar_ican", "mho/cm2",
 "i_ican", "mA/cm2",
 "tau_m_ican", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "erev_ican", &erev_ican,
 "beta_ican", &beta_ican,
 "cac_ican", &cac_ican,
 "taumin_ican", &taumin_ican,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ican",
 "gbar_ican",
 0,
 "i_ican",
 "m_inf_ican",
 "tau_m_ican",
 0,
 "m_ican",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gbar = 0.0001;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 prop_ion = need_memb(_na_sym);
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
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

 void _ican_sidi_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	ion_reg("na", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ican /home/fernando/S1_netpyne/sim/mod/ican_sidi.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Slow Ca-dependent cation current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(_threadargsprotocomma_ double, double);
static int iassign(_threadargsproto_);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
static int  iassign ( _threadargsproto_ ) {
   g = gbar * m * m ;
   i = g * ( v - erev ) ;
   ina = 0.7 * i ;
    return 0; }
 
static void _hoc_iassign(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 iassign ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   evaluate_fct ( _threadargscomma_ v , cai ) ;
   Dm = ( m_inf - m ) / tau_m ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 evaluate_fct ( _threadargscomma_ v , cai ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_m )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   evaluate_fct ( _threadargscomma_ v , cai ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_m)))*(- ( ( ( m_inf ) ) / tau_m ) / ( ( ( ( - 1.0 ) ) ) / tau_m ) - m) ;
   }
  return 0;
}
 
static int  evaluate_fct ( _threadargsprotocomma_ double _lv , double _lcai ) {
   double _lalpha2 ;
 _lalpha2 = beta * pow( ( _lcai / cac ) , 2.0 ) ;
   tau_m = 1.0 / ( _lalpha2 + beta ) / tadj ;
   m_inf = _lalpha2 / ( _lalpha2 + beta ) ;
   if ( tau_m < taumin ) {
     tau_m = taumin ;
     }
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
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
  cai = _ion_cai;
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
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  m = m0;
 {
   tadj = pow( 3.0 , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v , cai ) ;
   m = m_inf ;
   iassign ( _threadargs_ ) ;
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
  cai = _ion_cai;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   iassign ( _threadargs_ ) ;
   }
 _current += i;
 _current += ina;

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
  cai = _ion_cai;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
  cai = _ion_cai;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/fernando/S1_netpyne/sim/mod/ican_sidi.mod";
static const char* nmodl_file_text = 
  "TITLE Slow Ca-dependent cation current\n"
  ": from\n"
  ":  https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=144089&file=/PFCcell/mechanism/ican.mod\n"
  ":\n"
  ":   Ca++ dependent nonspecific cation current ICAN\n"
  ":   Differential equations\n"
  ":\n"
  ":   Model based on a first order kinetic scheme\n"
  ":\n"
  ":       + n cai <->     (alpha,beta)\n"
  ":\n"
  ":   Following this model, the activation fct will be half-activated at \n"
  ":   a concentration of Cai = (beta/alpha)^(1/n) = cac (parameter)\n"
  ":\n"
  ":   The mod file is here written for the case n=2 (2 binding sites)\n"
  ":   ---------------------------------------------\n"
  ":\n"
  ":   Kinetics based on: Partridge & Swandulla, TINS 11: 69-72, 1988.\n"
  ":\n"
  ":   This current has the following properties:\n"
  ":      - inward current (non specific for cations Na, K, Ca, ...)\n"
  ":      - activated by intracellular calcium\n"
  ":      - NOT voltage dependent\n"
  ":\n"
  ":   A minimal value for the time constant has been added\n"
  ":\n"
  ":   Ref: Destexhe et al., J. Neurophysiology 72: 803-818, 1994.\n"
  ":   See also:  http://www.cnl.salk.edu/~alain , http://cns.fmed.ulaval.ca\n"
  ":\n"
  "\n"
  ": Updated by Kiki Sidiropoulou (2010) so that dADP has slow inactivation kinetics and it\n"
  ": is activated after 5 spikes\n"
  "\n"
  ": Updated by Sam Neymotin (2016) to avoid using n ion and get rid of 'mystart' rule; also\n"
  ": make sure that INITIAL block assigns currents\n"
  "\n"
  "NEURON {\n"
  "  SUFFIX ican\n"
  "  NONSPECIFIC_CURRENT i\n"
  "  USEION ca READ cai\n"
  "  USEION na WRITE ina\n"
  "  RANGE gbar, m_inf, tau_m\n"
  "  GLOBAL beta, cac, taumin\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "  (mA) = (milliamp)\n"
  "  (mV) = (millivolt)\n"
  "  (molar) = (1/liter)\n"
  "  (mM) = (millimolar)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  v               (mV)\n"
  "  celsius = 36    (degC)\n"
  "  erev      = -20   (mV)            	: reversal potential\n"
  "  cai     	     (mM)           	: initial [Ca]i\n"
  "  gbar    = 0.0001 (mho/cm2)\n"
  "  beta    = 0.0001 (1/ms) 	 	: backward rate constant\n"
  "  cac     = 0.0004 (mM)\n"
  "  : middle point of activation fct, for ip3 as somacar, for current injection\n"
  "  taumin  = 0.1   (ms)            	: minimal value of time constant\n"
  "}\n"
  "\n"
  "STATE {\n"
  "  m\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  i      (mA/cm2)\n"
  "  ina     (mA/cm2)\n"
  "  m_inf\n"
  "  tau_m   (ms)\n"
  "  tadj\n"
  "  g (mho/cm2)\n"
  "}\n"
  "\n"
  "PROCEDURE iassign () {\n"
  "  g = gbar * m * m\n"
  "  i = g * (v - erev) \n"
  "  ina = 0.7 * i\n"
  "}\n"
  "\n"
  "BREAKPOINT { \n"
  "  SOLVE states METHOD cnexp  \n"
  "  iassign()\n"
  "}\n"
  "\n"
  "DERIVATIVE states { \n"
  "  evaluate_fct(v,cai)  \n"
  "  m' = (m_inf - m) / tau_m\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "INITIAL {\n"
  "  :  activation kinetics are assumed to be at 22 deg. C\n"
  "  :  Q10 is assumed to be 3\n"
  "  tadj = 3.0 ^ ((celsius-22.0)/10)  \n"
  "  evaluate_fct(v,cai)\n"
  "  m = m_inf\n"
  "  iassign()\n"
  "}\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2  \n"
  "  alpha2 = beta * (cai/cac)^2  \n"
  "  tau_m = 1 / (alpha2 + beta) / tadj\n"
  "  m_inf = alpha2 / (alpha2 + beta)  \n"
  "  if(tau_m < taumin) { tau_m = taumin } : min value of time constant\n"
  "}\n"
  "UNITSON\n"
  ;
#endif
