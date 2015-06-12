#include <faunus/faunus.h>
/**
 * @brief  Multiple polymers/peptides/proteins in bulk. Coarse-grained constant-pH, NVT, 
 *         MC simulation with implicit solvent and salt. 
 *
 * @note   Includes new non-bonded routine to handle the special case where two potentials 
 *         of the same type are required for different types of interactions.
 *         For example: Two square well potentials, one applied to intra-protein interactions 
 *         while the other handles inter-protein interactions.
 *
 * @author Joao Henriques
 * @date   2015/06/12 
 */
//
using namespace Faunus;
typedef Space<Geometry::Cuboid> Tspace;
//
#ifdef SPECIAL
namespace Faunus {
  namespace Energy {
    template<class Tspace, class Tpairpot, class Tpairpot_g2g>
    class NonbondedSpecial : public Nonbonded<Tspace, Tpairpot> {
    private:
      typedef Nonbonded<Tspace, Tpairpot> base;
      Nonbonded<Tspace, Tpairpot_g2g> nb;
    public:
      NonbondedSpecial(InputMap &in) : base(in), nb(in) {
	base::pairpot = Tpairpot(in, "squarewell");
	base::name = "Nonbonded N" + textio::squared + " - " + base::pairpot.name + " (intra)";
	nb.pairpot = Tpairpot_g2g(in, "squarewell_g2g");
	nb.name = "Nonbonded N" + textio::squared + " - " + nb.pairpot.name + " (inter)";
      }
      double g2g(const typename base::Tpvec &p, Group &g1, Group &g2) FOVERRIDE {
	assert(&g1 != &g2 && "Groups must be different");
	return nb.g2g(p, g1, g2);
      }
      double i2g(const typename base::Tpvec &p, Group &g, int i) FOVERRIDE {
	if (&g != base::spc->findGroup(i))
	  return nb.i2g(p, g, i);
	return base::i2g(p, g, i);
      }
      double i2all(typename base::Tpvec &p, int i) FOVERRIDE {
	double u = 0;
	for (auto g : base::spc->groupList())
	  u += i2g(p, *g, i);
	return u;
      }
      string info() { return base::info() + nb.info(); }
    };
  }
}
#endif
//
int main() {
  cout << textio::splash();
  //
  ::atom.includefile("param.json");
  InputMap mcp("simulation.input");  
  Tspace spc(mcp);
  //
#ifdef SPECIAL
  auto pot = Energy::Nonbonded<Tspace, Potential::DebyeHuckel>(mcp)
           + Energy::Nonbonded<Tspace, Potential::LennardJonesR12>(mcp)
           + Energy::NonbondedSpecial<Tspace, Potential::SquareWellHydrophobic, Potential::SquareWellHydrophobic>(mcp)
           + Energy::EquilibriumEnergy<Tspace>(mcp)
           + Energy::Bonded<Tspace>();
#else
  auto pot = Energy::Nonbonded<Tspace, Potential::DebyeHuckel>(mcp)
           + Energy::Nonbonded<Tspace, Potential::LennardJonesR12>(mcp)
           + Energy::Nonbonded<Tspace, Potential::SquareWellHydrophobic>(mcp)
           + Energy::EquilibriumEnergy<Tspace>(mcp)
           + Energy::Bonded<Tspace>();
#endif
  //
  auto bonded = &pot.second;
  //
  vector<Group> pepvec(mcp.get("peptide_N", 0));
  string file = mcp.get<string>("molecule", "");
  double req  = mcp.get<double>("harmonic_eqdist", 0);
  double k    = mcp.get<double>("harmonic_forceconst", 0);
  for (auto &element : pepvec) {
    Tspace::ParticleVector v;
    FormatAAM::load(file, v);
    Geometry::FindSpace().find(spc.geo, spc.p, v);
    element = spc.insert(v);
    element.name = "peptide";
    spc.enroll(element);
    for (int i = element.front(); i < element.back(); i++)
      bonded->add(i, i+1, Potential::Harmonic(k, req) - Potential::LennardJonesR12(mcp, "r12repex"));
  }
  //
  Group peptides(pepvec.front().front(), pepvec.back().back());
  //
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::TranslateRotate<Tspace> gmv(mcp, pot, spc);
  Move::CrankShaft<Tspace> crank(mcp, pot, spc);
  Move::Pivot<Tspace> pivot(mcp, pot, spc);
  Move::Reptation<Tspace> rep(mcp, pot, spc);
  Move::SwapMove<Tspace> tit(mcp, pot, spc, pot.first.second);
  //
  Analysis::PolymerShape shape;
  Average<double> rg2;
  Average<double> re;
  //
  Analysis::LineDistribution<> p2p_dist(0.2);
  Analysis::ChargeMultipole mp;
  //
  spc.load("simulation.state");
  //
  FormatGRO gro;
  gro.len=spc.geo.len.x();
  //
  EnergyDrift sys;
  sys.init(Energy::systemEnergy(spc, pot, spc.p));
  //
  std::ofstream rg_step("rg_step.dat");
  std::ofstream re_step("re_step.dat");
  //
  cout << atom.info()
       << spc.info() 
       << pot.info() 
       << "\n  ------------------------------------------------------------------------------------\n"
       << "  Louis: Hey, keep your f*cking mouth shut, all right? I mean it not one f*cking word!\n"
       << "  Melanie: Okay, Louis...\n"                                       
       << "  [Louis pulls a gun and shoots Melanie twice]\n"
       << "  ------------------------------------------------------------------------------------\n\n";
  //
  MCLoop loop(mcp);
  while (loop[0]) {
    while (loop[1]) {
      int k, i = slp_global.rand() % 6;
      switch (i) {
      // Particle translation.
      case 0:
	{
      	mv.setGroup(peptides);
	sys += mv.move(peptides.size());
	break;
	}
      // Group translation/rotation.
      case 1:
	k = pepvec.size();
	while (k-- > 0) { // postfix decrements after use
	  gmv.setGroup( pepvec[slp_global.rand() % pepvec.size()] );
	  sys += gmv.move();
	}
        break;
      // Crankshaft move.
      case 2:
	k = pepvec.size();
	while (k-- > 0) {
	  crank.setGroup( pepvec[slp_global.rand() % pepvec.size()] );
	  sys += crank.move();
	}
        break;
      // Pivot move. 
      case 3:
	k = pepvec.size();
	while (k-- > 0) {
	  pivot.setGroup( pepvec[slp_global.rand() % pepvec.size()] );
	  sys += pivot.move();
	}
	break;
      // Reptation move. 
      case 4:
      	k = pepvec.size();
	while (k-- > 0) {
	  rep.setGroup( pepvec[slp_global.rand() % pepvec.size()] );
	  sys+=rep.move();
	}
	break;	
      // Titration move.
      case 5:
	sys += tit.move();
	mp.sample(pepvec, spc);
	break;
      }
      // (Average) individual shape and Rg/Ree vs. step.
      for (auto &element : pepvec) {
	element.setMassCenter(spc);
	shape.sample(element, spc);
	Point pt = shape.vectorgyrationRadiusSquared(element, spc);
	double val = pt.x() + pt.y() + pt.z();
	rg2.add(val);
	auto i = element.front();
	auto j = element.back();
	double dist = spc.geo.dist(spc.p[i], spc.p[j]);	
	re.add(dist);
      }
      // Peptide-peptide center of mass distance distribution.
      for (auto i = pepvec.begin(); i != pepvec.end()-1; i++)
	for (auto j = i+1; j != pepvec.end(); j++) {
	  double dist = spc.geo.dist(i->cm, j->cm);
	  p2p_dist(dist)++;	    
	}
      // End of loop[1].  
    }
    rg_step << loop.innerCount() << " " << sqrt(rg2.avg()) << endl;
    re_step << loop.innerCount() << " " << re.avg() << endl;
    sys.checkDrift(Energy::systemEnergy(spc, pot, spc.p));
    cout << loop.timing();
    gro.save("simulation.gro", spc.p, "append");
    // End of loop[0].
  }
  p2p_dist.save("p2p_dist.dat");
  spc.save("simulation.state");  
  FormatPQR::save("simulation.pqr", spc.p);
  //
  cout << sys.info() 
       << mv.info() 
       << gmv.info() 
       << pivot.info() 
       << crank.info()
       << rep.info()
       << tit.info()
       << mp.info()
       << shape.info()
       << loop.info();
}
