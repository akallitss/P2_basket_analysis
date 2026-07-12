// P2 timing simulation inputs, computed with Garfield++ on lxplus:
//   mode "gas"  : Magboltz scan of Ar/iC4H10 95:5 -> gas_table.csv
//                 (E [V/cm], drift velocity, longitudinal/transverse diffusion)
//   mode "heed" : HEED primary-ionisation clusters for cosmic muons crossing
//                 a 3 mm gap at several zenith angles -> heed_clusters.csv
// Both CSVs are consumed by cosmic_bench_analysis/14_timing_simulation.py.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "Garfield/ComponentConstant.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "gas";

  MediumMagboltz gas;
  gas.SetComposition("ar", 95., "ic4h10", 5.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);

  if (mode == "gas") {
    const std::vector<double> efields = {100,  150,  200,  300,  400,
                                         500,  567,  700,  850,  1000,
                                         1250, 1500, 2000, 2500, 3000};
    const int ncoll = 5;  // x 1e7 collisions per point
    std::ofstream out("gas_table.csv");
    out << "E_Vcm,vz,vz_err,dl,dl_err,dt,dt_err,alpha,eta\n";
    for (const double e : efields) {
      double vx, vy, vz, wv, wr, dl, dt, alpha, eta, riontof, ratttof, lor;
      double vxerr, vyerr, vzerr, wverr, wrerr, dlerr, dterr;
      double alphaerr, etaerr, riontoferr, ratttoferr, lorerr, alphatof;
      std::array<double, 6> difftens;
      gas.RunMagboltz(e, 0., 0., ncoll, false, vx, vy, vz, wv, wr, dl, dt,
                      alpha, eta, riontof, ratttof, lor, vxerr, vyerr, vzerr,
                      wverr, wrerr, dlerr, dterr, alphaerr, etaerr, riontoferr,
                      ratttoferr, lorerr, alphatof, difftens);
      out << e << "," << vz << "," << vzerr << "," << dl << "," << dlerr << ","
          << dt << "," << dterr << "," << alpha << "," << eta << "\n";
      out.flush();
      std::cout << "E=" << e << " V/cm  vz=" << vz << " (+-" << vzerr
                << "%)  dl=" << dl << " dt=" << dt << std::endl;
    }
    return 0;
  }

  // ---- mode "heed": primary clusters in a 3 mm slab ----------------------
  GeometrySimple geo;
  // half-lengths [cm]; gap spans z = 0 .. 0.3 (mesh at z = 0)
  SolidBox box(0., 0., 0.15, 10., 10., 0.15);
  geo.AddSolid(&box, &gas);
  ComponentConstant cmp;
  cmp.SetGeometry(&geo);
  cmp.SetElectricField(0., 0., 100.);  // irrelevant for HEED clustering
  Sensor sensor;
  sensor.AddComponent(&cmp);

  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle("mu-");
  track.SetMomentum(4.e9);  // 4 GeV/c, typical cosmic muon

  const std::vector<double> angles = {0., 10., 20., 30., 40., 50.};
  const int ntracks = 4000;
  std::ofstream out("heed_clusters.csv");
  out << "track,angle_deg,z_cm,ne\n";
  for (const double ang : angles) {
    const double rad = ang * M_PI / 180.;
    for (int i = 0; i < ntracks; ++i) {
      // enter at the top face, direction tilted by the zenith angle
      const double x0 = -0.3 * std::tan(rad);
      if (!track.NewTrack(x0, 0., 0.2999, 0., std::sin(rad), 0., -std::cos(rad)))
        continue;
      double xc, yc, zc, tc, ec, extra;
      int ne, ni;
      while (track.GetCluster(xc, yc, zc, tc, ne, ni, ec, extra)) {
        out << i << "," << ang << "," << zc << "," << ne << "\n";
      }
    }
    std::cout << "HEED angle " << ang << " deg done" << std::endl;
  }
  return 0;
}
