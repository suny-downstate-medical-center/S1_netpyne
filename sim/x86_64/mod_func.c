#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _ar_traub_reg(void);
extern void _cadad_reg(void);
extern void _CaDynamics_E2_reg(void);
extern void _cadyn_reg(void);
extern void _cagk_reg(void);
extern void _Ca_HVA_reg(void);
extern void _cal_mh_reg(void);
extern void _cal_mig_reg(void);
extern void _Ca_LVAst_reg(void);
extern void _Ca_reg(void);
extern void _cancr_reg(void);
extern void _canin_reg(void);
extern void _can_mig_reg(void);
extern void _catcb_reg(void);
extern void _cat_mig_reg(void);
extern void _cat_traub_reg(void);
extern void _ch_CavL_reg(void);
extern void _ch_CavN_reg(void);
extern void _ch_KCaS_reg(void);
extern void _ch_Kdrfastngf_reg(void);
extern void _ch_KvAngf_reg(void);
extern void _ch_KvCaB_reg(void);
extern void _ch_leak_reg(void);
extern void _ch_Navngf_reg(void);
extern void _DetAMPANMDA_reg(void);
extern void _DetGABAAB_reg(void);
extern void _gabab_reg(void);
extern void _h_BS_reg(void);
extern void _HCN1_reg(void);
extern void _HH2_reg(void);
extern void _h_harnett_reg(void);
extern void _hin_reg(void);
extern void _h_kole_reg(void);
extern void _h_migliore_reg(void);
extern void _htc_reg(void);
extern void _ican_sidi_reg(void);
extern void _iccr_reg(void);
extern void _IC_reg(void);
extern void _iconc_Ca_reg(void);
extern void _Ih_reg(void);
extern void _ikscr_reg(void);
extern void _IKsin_reg(void);
extern void _Im_reg(void);
extern void _IT2_reg(void);
extern void _IT_reg(void);
extern void _kap_BS_reg(void);
extern void _kapcb_reg(void);
extern void _kapin_reg(void);
extern void _kBK_reg(void);
extern void _kca_reg(void);
extern void _kctin_reg(void);
extern void _kdmc_BS_reg(void);
extern void _kdr_BS_reg(void);
extern void _kdrcr_reg(void);
extern void _kdrin_reg(void);
extern void _KdShu2007_reg(void);
extern void _kl_reg(void);
extern void _km_reg(void);
extern void _K_Pst_reg(void);
extern void _K_Tst_reg(void);
extern void _kv_reg(void);
extern void _MyExp2SynBB_reg(void);
extern void _my_exp2syn_reg(void);
extern void _MyExp2SynNMDABB_reg(void);
extern void _nafcr_reg(void);
extern void _nafx_reg(void);
extern void _Nap_Et2_reg(void);
extern void _nap_sidi_reg(void);
extern void _NaTa_t_reg(void);
extern void _NaTs2_t_reg(void);
extern void _nax_BS_reg(void);
extern void _naz_reg(void);
extern void _Nca_reg(void);
extern void _ProbAMPANMDA_EMS_reg(void);
extern void _ProbGABAAB_EMS_reg(void);
extern void _savedist_reg(void);
extern void _SK_E2_reg(void);
extern void _SKv3_1_reg(void);
extern void _StochKv_deterministic_reg(void);
extern void _StochKv_det_reg(void);
extern void _StochKv_reg(void);
extern void _tia_reg(void);
extern void _vecstim_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"mod//ar_traub.mod\"");
    fprintf(stderr," \"mod//cadad.mod\"");
    fprintf(stderr," \"mod//CaDynamics_E2.mod\"");
    fprintf(stderr," \"mod//cadyn.mod\"");
    fprintf(stderr," \"mod//cagk.mod\"");
    fprintf(stderr," \"mod//Ca_HVA.mod\"");
    fprintf(stderr," \"mod//cal_mh.mod\"");
    fprintf(stderr," \"mod//cal_mig.mod\"");
    fprintf(stderr," \"mod//Ca_LVAst.mod\"");
    fprintf(stderr," \"mod//Ca.mod\"");
    fprintf(stderr," \"mod//cancr.mod\"");
    fprintf(stderr," \"mod//canin.mod\"");
    fprintf(stderr," \"mod//can_mig.mod\"");
    fprintf(stderr," \"mod//catcb.mod\"");
    fprintf(stderr," \"mod//cat_mig.mod\"");
    fprintf(stderr," \"mod//cat_traub.mod\"");
    fprintf(stderr," \"mod//ch_CavL.mod\"");
    fprintf(stderr," \"mod//ch_CavN.mod\"");
    fprintf(stderr," \"mod//ch_KCaS.mod\"");
    fprintf(stderr," \"mod//ch_Kdrfastngf.mod\"");
    fprintf(stderr," \"mod//ch_KvAngf.mod\"");
    fprintf(stderr," \"mod//ch_KvCaB.mod\"");
    fprintf(stderr," \"mod//ch_leak.mod\"");
    fprintf(stderr," \"mod//ch_Navngf.mod\"");
    fprintf(stderr," \"mod//DetAMPANMDA.mod\"");
    fprintf(stderr," \"mod//DetGABAAB.mod\"");
    fprintf(stderr," \"mod//gabab.mod\"");
    fprintf(stderr," \"mod//h_BS.mod\"");
    fprintf(stderr," \"mod//HCN1.mod\"");
    fprintf(stderr," \"mod//HH2.mod\"");
    fprintf(stderr," \"mod//h_harnett.mod\"");
    fprintf(stderr," \"mod//hin.mod\"");
    fprintf(stderr," \"mod//h_kole.mod\"");
    fprintf(stderr," \"mod//h_migliore.mod\"");
    fprintf(stderr," \"mod//htc.mod\"");
    fprintf(stderr," \"mod//ican_sidi.mod\"");
    fprintf(stderr," \"mod//iccr.mod\"");
    fprintf(stderr," \"mod//IC.mod\"");
    fprintf(stderr," \"mod//iconc_Ca.mod\"");
    fprintf(stderr," \"mod//Ih.mod\"");
    fprintf(stderr," \"mod//ikscr.mod\"");
    fprintf(stderr," \"mod//IKsin.mod\"");
    fprintf(stderr," \"mod//Im.mod\"");
    fprintf(stderr," \"mod//IT2.mod\"");
    fprintf(stderr," \"mod//IT.mod\"");
    fprintf(stderr," \"mod//kap_BS.mod\"");
    fprintf(stderr," \"mod//kapcb.mod\"");
    fprintf(stderr," \"mod//kapin.mod\"");
    fprintf(stderr," \"mod//kBK.mod\"");
    fprintf(stderr," \"mod//kca.mod\"");
    fprintf(stderr," \"mod//kctin.mod\"");
    fprintf(stderr," \"mod//kdmc_BS.mod\"");
    fprintf(stderr," \"mod//kdr_BS.mod\"");
    fprintf(stderr," \"mod//kdrcr.mod\"");
    fprintf(stderr," \"mod//kdrin.mod\"");
    fprintf(stderr," \"mod//KdShu2007.mod\"");
    fprintf(stderr," \"mod//kl.mod\"");
    fprintf(stderr," \"mod//km.mod\"");
    fprintf(stderr," \"mod//K_Pst.mod\"");
    fprintf(stderr," \"mod//K_Tst.mod\"");
    fprintf(stderr," \"mod//kv.mod\"");
    fprintf(stderr," \"mod//MyExp2SynBB.mod\"");
    fprintf(stderr," \"mod//my_exp2syn.mod\"");
    fprintf(stderr," \"mod//MyExp2SynNMDABB.mod\"");
    fprintf(stderr," \"mod//nafcr.mod\"");
    fprintf(stderr," \"mod//nafx.mod\"");
    fprintf(stderr," \"mod//Nap_Et2.mod\"");
    fprintf(stderr," \"mod//nap_sidi.mod\"");
    fprintf(stderr," \"mod//NaTa_t.mod\"");
    fprintf(stderr," \"mod//NaTs2_t.mod\"");
    fprintf(stderr," \"mod//nax_BS.mod\"");
    fprintf(stderr," \"mod//naz.mod\"");
    fprintf(stderr," \"mod//Nca.mod\"");
    fprintf(stderr," \"mod//ProbAMPANMDA_EMS.mod\"");
    fprintf(stderr," \"mod//ProbGABAAB_EMS.mod\"");
    fprintf(stderr," \"mod//savedist.mod\"");
    fprintf(stderr," \"mod//SK_E2.mod\"");
    fprintf(stderr," \"mod//SKv3_1.mod\"");
    fprintf(stderr," \"mod//StochKv_deterministic.mod\"");
    fprintf(stderr," \"mod//StochKv_det.mod\"");
    fprintf(stderr," \"mod//StochKv.mod\"");
    fprintf(stderr," \"mod//tia.mod\"");
    fprintf(stderr," \"mod//vecstim.mod\"");
    fprintf(stderr, "\n");
  }
  _ar_traub_reg();
  _cadad_reg();
  _CaDynamics_E2_reg();
  _cadyn_reg();
  _cagk_reg();
  _Ca_HVA_reg();
  _cal_mh_reg();
  _cal_mig_reg();
  _Ca_LVAst_reg();
  _Ca_reg();
  _cancr_reg();
  _canin_reg();
  _can_mig_reg();
  _catcb_reg();
  _cat_mig_reg();
  _cat_traub_reg();
  _ch_CavL_reg();
  _ch_CavN_reg();
  _ch_KCaS_reg();
  _ch_Kdrfastngf_reg();
  _ch_KvAngf_reg();
  _ch_KvCaB_reg();
  _ch_leak_reg();
  _ch_Navngf_reg();
  _DetAMPANMDA_reg();
  _DetGABAAB_reg();
  _gabab_reg();
  _h_BS_reg();
  _HCN1_reg();
  _HH2_reg();
  _h_harnett_reg();
  _hin_reg();
  _h_kole_reg();
  _h_migliore_reg();
  _htc_reg();
  _ican_sidi_reg();
  _iccr_reg();
  _IC_reg();
  _iconc_Ca_reg();
  _Ih_reg();
  _ikscr_reg();
  _IKsin_reg();
  _Im_reg();
  _IT2_reg();
  _IT_reg();
  _kap_BS_reg();
  _kapcb_reg();
  _kapin_reg();
  _kBK_reg();
  _kca_reg();
  _kctin_reg();
  _kdmc_BS_reg();
  _kdr_BS_reg();
  _kdrcr_reg();
  _kdrin_reg();
  _KdShu2007_reg();
  _kl_reg();
  _km_reg();
  _K_Pst_reg();
  _K_Tst_reg();
  _kv_reg();
  _MyExp2SynBB_reg();
  _my_exp2syn_reg();
  _MyExp2SynNMDABB_reg();
  _nafcr_reg();
  _nafx_reg();
  _Nap_Et2_reg();
  _nap_sidi_reg();
  _NaTa_t_reg();
  _NaTs2_t_reg();
  _nax_BS_reg();
  _naz_reg();
  _Nca_reg();
  _ProbAMPANMDA_EMS_reg();
  _ProbGABAAB_EMS_reg();
  _savedist_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
  _StochKv_deterministic_reg();
  _StochKv_det_reg();
  _StochKv_reg();
  _tia_reg();
  _vecstim_reg();
}
