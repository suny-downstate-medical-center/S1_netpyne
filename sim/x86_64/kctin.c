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
 
#define nrn_init _nrn_init__kctin
#define _nrn_initial _nrn_initial__kctin
#define nrn_cur _nrn_cur__kctin
#define _nrn_current _nrn_current__kctin
#define nrn_jacob _nrn_jacob__kctin
#define nrn_state _nrn_state__kctin
#define _net_receive _net_receive__kctin 
#define kin kin__kctin 
#define rates rates__kctin 
 
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
#define gkcbar _p[0]
#define ik _p[1]
#define cst _p[2]
#define ost _p[3]
#define ist _p[4]
#define cai _p[5]
#define ek _p[6]
#define k1 _p[7]
#define k2 _p[8]
#define k3 _p[9]
#define k4 _p[10]
#define q10 _p[11]
#define Dcst _p[12]
#define Dost _p[13]
#define Dist _p[14]
#define v _p[15]
#define _g _p[16]
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
 static void _hoc_alpha(void);
 static void _hoc_alp(void);
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
 "setdata_kctin", _hoc_setdata,
 "alpha_kctin", _hoc_alpha,
 "alp_kctin", _hoc_alp,
 "rates_kctin", _hoc_rates,
 0, 0
};
#define alpha alpha_kctin
#define alp alp_kctin
 extern double alpha( _threadargsprotocomma_ double , double , double , double , double );
 extern double alp( _threadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gkcbar_kctin", "S/cm2",
 "ik_kctin", "mA/cm2",
 0,0
};
 static double cst0 = 0;
 static double delta_t = 0.01;
 static double ist0 = 0;
 static double ost0 = 0;
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
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"kctin",
 "gkcbar_kctin",
 0,
 "ik_kctin",
 0,
 "cst_kctin",
 "ost_kctin",
 "ist_kctin",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	gkcbar = 0;
 	_prop->param = _p;
 	_prop->param_size = 17;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
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
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _kctin_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 3);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 17, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kctin /home/fernando/S1_netpyne/sim/mod/kctin.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Kct current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double, double);
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3]; static double *_temp1;
 static int kin();
 
static int kin (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<3;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ v , cai ) ;
   /* ~ cst <-> ost ( k3 , k4 )*/
 f_flux =  k3 * cst ;
 b_flux =  k4 * ost ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  k3 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 2 ,1)  -= _term;
 _term =  k4 ;
 _MATELM1( 1 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ ost <-> ist ( k1 , 0.0 )*/
 f_flux =  k1 * ost ;
 b_flux =  0.0 * ist ;
 _RHS1( 2) -= (f_flux - b_flux);
 
 _term =  k1 ;
 _MATELM1( 2 ,2)  += _term;
 _term =  0.0 ;
 _MATELM1( 2 ,0)  -= _term;
 /*REACTION*/
  /* ~ ist <-> cst ( k2 , 0.0 )*/
 f_flux =  k2 * ist ;
 b_flux =  0.0 * cst ;
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  k2 ;
 _MATELM1( 1 ,0)  -= _term;
 _term =  0.0 ;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
   /* cst + ost + ist = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= ist ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= ost ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= cst ;
 /*CONSERVATION*/
   } return _reset;
 }
 
static int  rates ( _threadargsprotocomma_ double _lv , double _lcai ) {
   k1 = alp ( _threadargscomma_ 0.01 , _lv , - 10.0 , 1.0 ) ;
   k2 = alp ( _threadargscomma_ 0.1 , _lv , - 120.0 , - 10.0 ) ;
   k3 = alpha ( _threadargscomma_ 0.001 , 1.0 , _lv , - 20.0 , 7.0 ) * 1.0e8 * pow( ( _lcai * 1.0 ) , 3.0 ) ;
   k4 = alp ( _threadargscomma_ 0.2 , _lv , - 44.0 , - 5.0 ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double alpha ( _threadargsprotocomma_ double _ltmin , double _ltmax , double _lv , double _lvhalf , double _lk ) {
   double _lalpha;
 _lalpha = 1.0 / ( _ltmin + 1.0 / ( 1.0 / ( _ltmax - _ltmin ) + exp ( ( _lv - _lvhalf ) / _lk ) * 1.0 ) ) ;
   
return _lalpha;
 }
 
static void _hoc_alpha(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alpha ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) , *getarg(5) );
 hoc_retpushx(_r);
}
 
double alp ( _threadargsprotocomma_ double _ltmin , double _lv , double _lvhalf , double _lk ) {
   double _lalp;
 _lalp = 1.0 / ( _ltmin + exp ( - ( _lv - _lvhalf ) / _lk ) * 1.0 ) ;
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alp ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<3;_i++) _p[_dlist1[_i]] = 0.0;}
 rates ( _threadargscomma_ v , cai ) ;
 /* ~ cst <-> ost ( k3 , k4 )*/
 f_flux =  k3 * cst ;
 b_flux =  k4 * ost ;
 Dcst -= (f_flux - b_flux);
 Dost += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ost <-> ist ( k1 , 0.0 )*/
 f_flux =  k1 * ost ;
 b_flux =  0.0 * ist ;
 Dost -= (f_flux - b_flux);
 Dist += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ ist <-> cst ( k2 , 0.0 )*/
 f_flux =  k2 * ist ;
 b_flux =  0.0 * cst ;
 Dist -= (f_flux - b_flux);
 Dcst += (f_flux - b_flux);
 
 /*REACTION*/
   /* cst + ost + ist = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<3;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ v , cai ) ;
 /* ~ cst <-> ost ( k3 , k4 )*/
 _term =  k3 ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 2 ,1)  -= _term;
 _term =  k4 ;
 _MATELM1( 1 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ ost <-> ist ( k1 , 0.0 )*/
 _term =  k1 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 0 ,2)  -= _term;
 _term =  0.0 ;
 _MATELM1( 2 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ ist <-> cst ( k2 , 0.0 )*/
 _term =  k2 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 1 ,0)  -= _term;
 _term =  0.0 ;
 _MATELM1( 0 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
   /* cst + ost + ist = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 3;}
 
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
  cai = _ion_cai;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 3, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
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
  cst = cst0;
  ist = ist0;
  ost = ost0;
 {
    _ss_sparse_thread(&_thread[_spth1]._pvoid, 3, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
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
   ik = gkcbar * ost * ( v - ek ) ;
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
 {  sparse_thread(&_thread[_spth1]._pvoid, 3, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(ist) - _p;  _dlist1[0] = &(Dist) - _p;
 _slist1[1] = &(cst) - _p;  _dlist1[1] = &(Dcst) - _p;
 _slist1[2] = &(ost) - _p;  _dlist1[2] = &(Dost) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/fernando/S1_netpyne/sim/mod/kctin.mod";
static const char* nmodl_file_text = 
  "\n"
  "TITLE Kct current\n"
  "\n"
  "COMMENT Equations from \n"
  "		  Shao L.R., Halvorsrud R., Borg-Graham L., Storm J.F. The role of BK-type Ca2_-dependent K+ channels in spike broadening during repetitive firing in rat hippocampal pyramidal cells J.Physiology (1999),521:135-146 \n"
  "		  \n"
  "		  The Krasnow Institute\n"
  "		  George Mason University\n"
  "\n"
  "Copyright	  Maciej Lazarewicz, 2001\n"
  "		  (mlazarew@gmu.edu)\n"
  "		  All rights reserved.\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	:SUFFIX  mykca\n"
  "	SUFFIX kctin\n"
  "	USEION k READ ek WRITE ik\n"
  "	USEION ca READ cai   \n"
  "	RANGE  gkcbar,ik\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mM)	= (millimolar)\n"
  "	(S)  	= (siemens)\n"
  "	(mA) 	= (milliamp)\n"
  "	(mV) 	= (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "        gkcbar	= 0	 (S/cm2)\n"
  "	:gkbar	= 1.0e-3 (S/cm2)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "        v       (mV)\n"
  "        cai	(mM)     :26/10/04   htan se sxolio\n"
  "	celsius		 (degC)\n"
  "	ek		(mV)\n"
  "	ik	(mA/cm2)\n"
  "	k1	(/ms)\n"
  "	k2	(/ms)\n"
  "	k3	(/ms)\n"
  "	k4	(/ms)\n"
  "	q10	(1)\n"
  "}\n"
  "\n"
  "STATE { cst ost ist }\n"
  "\n"
  "BREAKPOINT { \n"
  "	SOLVE kin METHOD sparse\n"
  "	ik = gkcbar * ost *( v - ek ) \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	SOLVE kin STEADYSTATE sparse\n"
  "}\n"
  "\n"
  "KINETIC kin {\n"
  ":	rates(v, cani)\n"
  "	rates(v, cai)\n"
  "	~cst<->ost  (k3,k4)\n"
  "	~ost<->ist  (k1,0.0)\n"
  "	~ist<->cst  (k2,0.0)\n"
  "	CONSERVE cst+ost+ist=1\n"
  "}\n"
  "\n"
  ":change feb8th for pfc\n"
  ":PROCEDURE rates( v(mV), cani(mM)) {\n"
  "PROCEDURE rates( v(mV), cai(mM)) {\n"
  ":	 k1=alp( 0.1, v,  -10.0,   1.0 ) : original\n"
  "	 k1=alp( 0.01, v,  -10.0,   1.0 ) :increases the current\n"
  "	 k2=alp( 0.1, v, -120.0, -10.0 ) :original (0.1, -120, -10)\n"
  ":	 k3=alpha( 0.001, 1.0, v, -20.0, 7.0 ) *1.0e8* ( cai*1.0(/mM) )^3  :original\n"
  "	 k3=alpha( 0.001, 1.0, v, -20.0, 7.0 ) *1.0e8* (cai*1.0(/mM) )^3\n"
  "	 :k3 changes the attenuation\n"
  "	 k4=alp( 0.2, v, -44.0,  -5.0 ) :original\n"
  ":	 k4=alp( 0.2, v, -44.0,  -5.0 )\n"
  "}\n"
  "\n"
  "FUNCTION alpha( tmin(ms), tmax(ms), v(mV), vhalf(mV), k(mV) )(/ms){\n"
  "        alpha = 1.0 / ( tmin + 1.0 / ( 1.0 / (tmax-tmin) + exp((v-vhalf)/k)*1.0(/ms) ) )\n"
  "}\n"
  "\n"
  "FUNCTION alp( tmin(ms), v(mV), vhalf(mV), k(mV) )(/ms){\n"
  "        alp = 1.0 / ( tmin + exp( -(v-vhalf) / k )*1.0(ms) )\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
