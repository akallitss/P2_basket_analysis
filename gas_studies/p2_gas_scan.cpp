// p2_gas_scan.cpp
// -----------------------------------------------------------------------------
// Magboltz (via Garfield++) transport scan for the P2 Micromegas, used to map
// operating high voltages between gas mixtures.
//
// For a gas mixture given on the command line it scans two field regions:
//   region "amp"   : high field in the 150 um amplification gap.  We scan the
//                    field that corresponds to a mesh voltage V_mesh = E*d_amp,
//                    d_amp = 0.015 cm, over V_mesh = 250..700 V.  The Townsend
//                    coefficient (alpha) and attachment (eta) give the gas gain
//                    G = exp[(alpha-eta) * d_amp].
//   region "drift" : low field in the 3 mm conversion gap.  Drift velocity vs
//                    E; the drift HV difference across the gap is dV = E*d_drift,
//                    d_drift = 0.30 cm.
//
// Output CSV columns (one file per gas):
//   region,E_Vcm,vz_cm_ns,vz_err_pct,dl,dl_err,dt,dt_err,alpha_cm,eta_cm,
//   alphatof_cm   (alphatof = effective Townsend (alpha-eta) via time-of-flight)
// The Python side (analyze.py) converts E<->voltage using the geometry above.
//
// Usage:
//   p2_gas_scan <out.csv> <gasA> <fracA> [<gasB> <fracB> ...]   (up to 6 gases)
// Example (NSW gas):
//   p2_gas_scan ar_co2_iso_95_3_2.csv ar 95 co2 3 ic4h10 2
//
// Built and run on lxplus by run_lxplus.sh against an LCG view (Garfield++).
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;

namespace {
const double kDAmp = 0.015;    // amplification gap [cm] = 150 um
const double kDDrift = 0.300;  // drift gap [cm] = 3 mm

// Magboltz SetComposition only has fixed-arity overloads; dispatch by size.
bool SetMixture(MediumMagboltz& gas, const std::vector<std::string>& n,
                const std::vector<double>& f) {
  switch (n.size()) {
    case 1: gas.SetComposition(n[0], f[0]); return true;
    case 2: gas.SetComposition(n[0], f[0], n[1], f[1]); return true;
    case 3: gas.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2]); return true;
    case 4: gas.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2],
                               n[3], f[3]); return true;
    case 5: gas.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2],
                               n[3], f[3], n[4], f[4]); return true;
    case 6: gas.SetComposition(n[0], f[0], n[1], f[1], n[2], f[2],
                               n[3], f[3], n[4], f[4], n[5], f[5]); return true;
    default: return false;
  }
}

// One Magboltz point -> append a CSV row.
void ScanPoint(MediumMagboltz& gas, std::ofstream& out, const std::string& reg,
               double e, int ncoll) {
  double vx, vy, vz, dl, dt, alpha, eta, lor;
  double vxerr, vyerr, vzerr, dlerr, dterr, alphaerr, etaerr, lorerr, alphatof;
  std::array<double, 6> difftens;
  gas.RunMagboltz(e, 0., 0., ncoll, false, vx, vy, vz, dl, dt, alpha, eta, lor,
                  vxerr, vyerr, vzerr, dlerr, dterr, alphaerr, etaerr, lorerr,
                  alphatof, difftens);
  out << reg << "," << e << "," << vz << "," << vzerr << "," << dl << ","
      << dlerr << "," << dt << "," << dterr << "," << alpha << "," << eta << ","
      << alphatof << "\n";
  out.flush();
  std::cout << "  [" << reg << "] E=" << e << " V/cm  vz=" << vz
            << " cm/ns  alpha=" << alpha << "  eta=" << eta
            << "  a_tof=" << alphatof << std::endl;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 4 || (argc % 2) != 0) {
    std::cerr << "usage: " << argv[0]
              << " <out.csv> <gasA> <fracA> [<gasB> <fracB> ...]\n";
    return 1;
  }
  const std::string outfile = argv[1];
  std::vector<std::string> names;
  std::vector<double> fracs;
  for (int i = 2; i + 1 < argc; i += 2) {
    names.push_back(argv[i]);
    fracs.push_back(std::atof(argv[i + 1]));
  }

  MediumMagboltz gas;
  if (!SetMixture(gas, names, fracs)) {
    std::cerr << "error: 1..6 gas components required\n";
    return 1;
  }
  gas.SetTemperature(293.15);  // 20 C
  gas.SetPressure(760.);       // 1 atm
  gas.EnableThermalMotion();

  std::cout << "Gas:";
  for (size_t i = 0; i < names.size(); ++i)
    std::cout << " " << names[i] << " " << fracs[i] << "%";
  std::cout << "  -> " << outfile << std::endl;

  std::ofstream out(outfile.c_str());
  out << "region,E_Vcm,vz_cm_ns,vz_err_pct,dl,dl_err,dt,dt_err,alpha_cm,eta_cm,"
         "alphatof_cm\n";

  // --- drift region FIRST: low-field grid over the 3 mm gap (fast) ---------
  // Deliberately before the amplification scan: these points cost ~15 s each,
  // so the drift-velocity / diffusion answer (the timing-resolution driver) is
  // on disk within a couple of minutes even if the slow Townsend scan below is
  // still running. The grid is expressed as the mesh->drift HV difference dV
  // actually dialled at the beam, dV = E * 3 mm; it brackets the present
  // working point (mesh 410 / drift 610 -> dV = 200 V -> E = 667 V/cm).
  const int ncoll_drift = 3;
  const std::vector<double> driftDV = {30.,  60.,  100., 150., 200., 250.,
                                       300., 400., 500., 600., 800.};
  for (const double dv : driftDV) ScanPoint(gas, out, "drift", dv / kDDrift, ncoll_drift);

  // --- amplification region: scan by mesh voltage across the 150 um gap ----
  // High-field Townsend (steady-state method) is the slow part of Magboltz
  // (minutes/point), so keep the grid coarse: the gain curve is smooth and
  // monotonic and analyze.py interpolates. Points bracket the 410 V working
  // point and extend to 780 V so a heavily quenched CF4 mixture (which needs
  // markedly more mesh voltage for the same gain) still has an equal-gain match.
  const int ncoll_amp = 2;  // x 1e7 collisions (Townsend needs statistics)
  const std::vector<double> ampV = {350., 420., 490., 560., 630., 700., 780.};
  for (const double vmesh : ampV) ScanPoint(gas, out, "amp", vmesh / kDAmp, ncoll_amp);

  std::cout << "done -> " << outfile << " (d_amp=" << kDAmp
            << " cm, d_drift=" << kDDrift << " cm)" << std::endl;
  return 0;
}
