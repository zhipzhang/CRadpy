#include "nanobind/nanobind.h"
#include "nanobind/nb_defs.h"
#include "constants.h"
namespace nb = nanobind;

NB_MODULE(CRadpy_constants, m)
{
    m.attr("TeV_to_erg") = constants::TeV_to_erg;
    m.attr("GeV_to_erg") = constants::GeV_to_erg;
    m.attr("eV_to_erg") = constants::eV_to_erg;
    m.attr("kb") = constants::kb;
    m.attr("sigma_T") = constants::sigma_T;
    m.attr("e_radius") = constants::e_radius;
    m.attr("pc_to_cm") = constants::pc_to_cm;
    m.attr("AU_to_cm") = constants::AU_to_cm;
    m.attr("m_p_g") = constants::m_p_g;
    m.attr("m_p") = constants::m_p;
    m.attr("m_e_g") = constants::m_e_g;
    m.attr("m_e") = constants::m_e;
    m.attr("m_pi") = constants::m_pi;
    m.attr("yr_to_sec") = constants::yr_to_sec;
    m.attr("pi") = constants::pi;
    m.attr("mSol") = constants::mSol;
    m.attr("c_speed") = constants::c_speed;
    m.attr("el_charge") = constants::el_charge;
    m.attr("eRadius") = constants::eRadius;
    m.attr("hp") = constants::hp;
    m.attr("fineStructConst") = constants::fineStructConst;
    m.attr("h_to_sec") = constants::h_to_sec;
    m.attr("pc_to_lyr") = constants::pc_to_lyr;
    m.attr("ln10") = constants::ln10;

}
