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
 
#define nrn_init _nrn_init__ch_KCaS
#define _nrn_initial _nrn_initial__ch_KCaS
#define nrn_cur _nrn_cur__ch_KCaS
#define _nrn_current _nrn_current__ch_KCaS
#define nrn_jacob _nrn_jacob__ch_KCaS
#define nrn_state _nrn_state__ch_KCaS
#define _net_receive _net_receive__ch_KCaS 
#define rate rate__ch_KCaS 
#define state state__ch_KCaS 
 
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
#define ik _p[1]
#define g _p[2]
#define qinf _p[3]
#define qtau _p[4]
#define myi _p[5]
#define q _p[6]
#define ek _p[7]
#define cai _p[8]
#define Dq _p[9]
#define qexp _p[10]
#define v _p[11]
#define _g _p[12]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
#define _ion_cai	*_ppvar[3]._pval
 
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
 static void _hoc_rate(void);
 static void _hoc_state(void);
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
 "setdata_ch_KCaS", _hoc_setdata,
 "rate_ch_KCaS", _hoc_rate,
 "state_ch_KCaS", _hoc_state,
 0, 0
};
 #define _zq10 _thread[0]._pval[0]
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_ch_KCaS", "mho/cm2",
 "ik_ch_KCaS", "mA/cm2",
 "g_ch_KCaS", "mho/cm2",
 "qtau_ch_KCaS", "ms",
 "myi_ch_KCaS", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double q0 = 0;
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ch_KCaS",
 "gmax_ch_KCaS",
 0,
 "ik_ch_KCaS",
 "g_ch_KCaS",
 "qinf_ch_KCaS",
 "qtau_ch_KCaS",
 "myi_ch_KCaS",
 0,
 "q_ch_KCaS",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* cai */
 
}
 static void _initlists();
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ch_KCaS_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", 1.0);
 	ion_reg("ca", 2.0);
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ch_KCaS /home/fernando/S1_netpyne/sim/mod/ch_KCaS.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 /*Top LOCAL _zq10 */
static int _reset;
static char *modelname = "calcium-activated potassium channel (non-voltage-dependent)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(_threadargsprotocomma_ double);
static int state(_threadargsproto_);
 
/*VERBATIM*/
#include <stdlib.h> /* 	Include this library so that the following
						(innocuous) warning does not appear:
						 In function '_thread_cleanup':
						 warning: incompatible implicit declaration of 
						          built-in function 'free'  */
 
static int  state ( _threadargsproto_ ) {
   rate ( _threadargscomma_ cai ) ;
   q = q + ( qinf - q ) * qexp ;
    return 0; }
 
static void _hoc_state(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 state ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
static int  rate ( _threadargsprotocomma_ double _lcai ) {
   double _lalpha , _lbeta , _ltinc ;
 _zq10 = pow( 3.0 , ( ( celsius - 34.0 ) / 10.0 ) ) ;
   _lalpha = 1.25e1 * _lcai * _lcai ;
   _lbeta = 0.00025 ;
   qtau = 1.0 / ( _lalpha + _lbeta ) / _zq10 ;
   qinf = _lalpha * qtau ;
   _ltinc = - dt ;
   qexp = 1.0 - exp ( _ltinc / qtau ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rate ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("ch_KCaS", "cannot be used with CVODE"); return 0;}
 
static void _thread_mem_init(Datum* _thread) {
   _thread[0]._pval = (double*)ecalloc(1, sizeof(double));
 }
 
static void _thread_cleanup(Datum* _thread) {
   free((void*)(_thread[0]._pval));
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  q = q0;
 {
   q = qinf ;
   rate ( _threadargscomma_ cai ) ;
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
  ek = _ion_ek;
  cai = _ion_cai;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gmax * q * q ;
   ik = g * ( v - ek ) ;
   myi = ik ;
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
  cai = _ion_cai;
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
  cai = _ion_cai;
 {  { state(_p, _ppvar, _thread, _nt); }
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/fernando/S1_netpyne/sim/mod/ch_KCaS.mod";
static const char* nmodl_file_text = 
  "TITLE calcium-activated potassium channel (non-voltage-dependent)\n"
  "\n"
  "COMMENT\n"
  "Ca2+ activated K+ channel (not voltage dependent)\n"
  "\n"
  "From:  original said for granule cells, but used in all the cell types\n"
  "\n"
  "Updates:\n"
  "2014 December (Marianne Bezaire): documented\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX ch_KCaS\n"
  "	USEION k READ ek WRITE ik VALENCE 1\n"
  "	USEION ca READ cai VALENCE 2\n"
  "	RANGE g, gmax, qinf, qtau, ik\n"
  "	RANGE myi\n"
  "    THREADSAFE\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "#include <stdlib.h> /* 	Include this library so that the following\n"
  "						(innocuous) warning does not appear:\n"
  "						 In function '_thread_cleanup':\n"
  "						 warning: incompatible implicit declaration of \n"
  "						          built-in function 'free'  */\n"
  "ENDVERBATIM\n"
  "\n"
  "UNITS {\n"
  "        (molar) = (1/liter)\n"
  "        (mM)    = (millimolar)\n"
  "	(mA)	= (milliamp)\n"
  "	(mV)	= (millivolt)\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "PARAMETER {\n"
  "      celsius (degC) : temperature - set in hoc; default is 6.3\n"
  "	v		(mV)\n"
  "	dt		(ms)\n"
  "	gmax  (mho/cm2)\n"
  "	ek	(mV)\n"
  "	cai (mM)\n"
  "}\n"
  "\n"
  "STATE { q }\n"
  "\n"
  "ASSIGNED {\n"
  "	ik (mA/cm2) \n"
  "	g (mho/cm2) \n"
  "	qinf \n"
  "	qtau (ms) \n"
  "	qexp\n"
  "	myi (mA/cm2)\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {          :Computes i=g*q^2*(v-ek)\n"
  "	SOLVE state\n"
  "    g = gmax * q*q\n"
  "	ik = g * (v-ek)\n"
  "	myi = ik\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  ": verbatim blocks are not thread safe (perhaps related, this mechanism cannot be used with cvode)\n"
  "INITIAL {\n"
  "	q=qinf\n"
  "	rate(cai)\n"
  "}\n"
  "\n"
  "PROCEDURE state() {  :Computes state variable q at current v and dt.\n"
  "	:cai = ncai + lcai + tcai\n"
  "	rate(cai)\n"
  "	q = q + (qinf-q) * qexp\n"
  "}\n"
  "\n"
  "LOCAL q10\n"
  "PROCEDURE rate(cai) {  :Computes rate and other constants at current v.\n"
  "	LOCAL alpha, beta, tinc\n"
  "	q10 = 3^((celsius - 34)/10) : set to 1 for the cutsuridis model?\n"
  "		:\"q\" activation system\n"
  "alpha = 1.25e1 * cai * cai\n"
  "beta = 0.00025 \n"
  "\n"
  "	qtau = 1 /(alpha + beta)/q10\n"
  "	qinf = alpha * qtau\n"
  "	tinc = -dt\n"
  "	qexp = 1 - exp(tinc/qtau)\n"
  "}\n"
  "\n"
  "UNITSON\n"
  "\n"
  ;
#endif
