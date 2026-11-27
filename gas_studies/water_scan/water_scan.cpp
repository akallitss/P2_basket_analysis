// water_scan.cpp
// -----------------------------------------------------------------------------
// Magboltz (via Garfield++ MediumMagboltz::RunMagboltz) drift-field scan for the
// det2 drift-velocity audit: v_d, D_L, D_T, Townsend alpha and attachment eta
// versus E (20 .. 1500 V/cm) for Ar/iC4H10 95/5 with contaminants (H2O, air)
// or a different isobutane fraction.
//
// Usage:
//   water_scan <out.csv> <ncoll_1e7> <T_K> <p_Torr> <gasA> <fracA> [<gasB> <fracB> ...]
//     (1..6 components, fractions in percent, normalised by Magboltz)
// Output CSV (one row per field, flushed per row):
//   E_Vcm,vz_cm_ns,vz_err_pct,dl_sqrtcm,dl_err_pct,dt_sqrtcm,dt_err_pct,
//   alpha_cm,alpha_err_pct,eta_cm,eta_err_pct
// Env WS_FIELDS="20,30,..." overrides the field list.
// -----------------------------------------------------------------------------
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;

static bool SetMixture(MediumMagboltz& g, const std::vector<std::string>& n,
                       const std::vector<double>& f) {
  switch (n.size()) {
    case 1: return g.SetComposition(n[0], f[0]);
    case 2: return g.SetComposition(n[0], f[0], n[1], f[1]);
    case 3: return g.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2]);
    case 4: return g.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2], n[3], f[3]);
    case 5: return g.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2], n[3], f[3], n[4], f[4]);
    case 6: return g.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2], n[3], f[3], n[4], f[4], n[5], f[5]);
    default: return false;
  }
}

int main(int argc, char** argv) {
  if (argc < 7 || (argc % 2) == 0) {
    std::cerr << "usage: " << argv[0]
              << " <out.csv> <ncoll_1e7> <T_K> <p_Torr> <gasA> <fracA> [...]\n";
    return 1;
  }
  const std::string outfile = argv[1];
  const int ncoll = std::atoi(argv[2]);
  const double tk = std::atof(argv[3]);
  const double ptorr = std::atof(argv[4]);
  std::vector<std::string> names;
  std::vector<double> fracs;
  for (int i = 5; i + 1 < argc; i += 2) {
    names.push_back(argv[i]);
    fracs.push_back(std::atof(argv[i + 1]));
  }
  std::vector<double> fields = {20.,  30.,  40.,  50.,  62.5, 75.,  87.5, 100.,
                                125., 150., 187.5, 225., 250., 312.5, 375.,
                                437.5, 500., 562.5, 625., 750., 875., 1000.,
                                1150., 1312.5, 1500.};
  if (const char* s = std::getenv("WS_FIELDS")) {
    fields.clear();
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) fields.push_back(std::atof(tok.c_str()));
  }

  MediumMagboltz gas;
  if (!SetMixture(gas, names, fracs)) {
    std::cerr << "SetComposition failed\n";
    return 2;
  }
  gas.SetTemperature(tk);
  gas.SetPressure(ptorr);
  gas.EnableThermalMotion();
  gas.PrintGas();

  std::ofstream out(outfile.c_str());
  out << "E_Vcm,vz_cm_ns,vz_err_pct,dl_sqrtcm,dl_err_pct,dt_sqrtcm,dt_err_pct,"
         "alpha_cm,alpha_err_pct,eta_cm,eta_err_pct\n";
  for (const double e : fields) {
    double vx, vy, vz, wv, wr, dl, dt, alpha, eta, rion, ratt, lor;
    double vxe, vye, vze, wve, wre, dle, dte, alphae, etae, rione, ratte, lore;
    double alphatof;
    std::array<double, 6> dten;
    // Garfield >= 2024 signature (with TOF outputs, unused here)
    gas.RunMagboltz(e, 0., 0., ncoll, false, vx, vy, vz, wv, wr, dl, dt, alpha,
                    eta, rion, ratt, lor, vxe, vye, vze, wve, wre, dle, dte,
                    alphae, etae, rione, ratte, lore, alphatof, dten);
    out << e << "," << vz << "," << vze << "," << dl << "," << dle << "," << dt
        << "," << dte << "," << alpha << "," << alphae << "," << eta << ","
        << etae << "\n";
    out.flush();
    std::cout << "E=" << e << " vz=" << vz << " cm/ns (" << vze << "%) dl=" << dl
              << " dt=" << dt << " eta=" << eta << " /cm" << std::endl;
  }
  return 0;
}
