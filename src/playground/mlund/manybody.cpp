#include <faunus/faunus.h>
#include <tclap/CmdLine.h>

using namespace Faunus;
using namespace TCLAP;

typedef Geometry::Cuboid Tgeometry;
typedef Potential::CombinedPairPotential<Potential::SquareWellHydrophobic, Potential::LennardJones> SRpot;
typedef Potential::CoulombSR<Tgeometry, Potential::DebyeHuckel, SRpot> Tpairpot;
//typedef Potential::CoulombSR<Tgeometry, Potential::DebyeHuckel, Potential::HardSphere> Tpairpot;

int main(int argc, char** argv) {
  string inputfile,istate,ostate;
  try {
    cout << textio::splash();
    CmdLine cmd("NPT Monte Carlo simulation of rigid bodies in continuum", ' ', "0.1");
    ValueArg<string> inputArg("i","inputfile","InputMap key/value file",true,"","inputfile");
    ValueArg<string> istateArg("c","instate","Name of input statefile",false,"state","instate");
    ValueArg<string> ostateArg("o","outstate","Name of output statefile",false,"state","outstate");
    cmd.add( inputArg );
    cmd.add( istateArg );
    cmd.add( ostateArg );
    cmd.parse( argc, argv );
    inputfile = inputArg.getValue();
    istate = istateArg.getValue();
    ostate = ostateArg.getValue();
  }
  catch (ArgException &e)  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  InputMap mcp(inputfile);
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  Space spc( pot.getGeometry() );

  // Add molecules
  int N1 = mcp.get("molecule_N1",0);
  int N2 = mcp.get("molecule_N2",0);
  string file;
  vector<GroupMolecular> pol(N1+N2);
  for (int i=0; i<N1+N2; i++) {
    GroupMolecular g;
    if (i>=N1)
      file = mcp.get<string>("molecule_file2", "");
    else
      file = mcp.get<string>("molecule_file1", "");
    aam.load(file);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.p);        // find empty spot in particle vector
    pol[i] = spc.insert( aam.p );          // insert into space
    pol[i].name=file;
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  spc.load(istate);

  Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::SwapMoveMSR tit(mcp,pot,spc);
  return 0;
  Analysis::RadialDistribution<float,int> rdf(0.25);
  Analysis::ChargeMultipole mpol;

  double utot=pot.external() + pot.g_internal(spc.p, allpol);
  for (auto &g : pol)
    utot += pot.g_external(spc.p, g);
  sys.init( utot );

  cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 3;
      switch (i) {
        case 0:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo->dist(i->cm,j->cm) )++;
          break;
        case 1:
          sys+=iso.move();
          break;
        case 2:
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
      }
      if ( slp_global.runtest(0.0001) ) {
        xtc.setbox( nonbonded->pair.geo.len );
        xtc.save("traj.xtc", spc);
      }
    } // end of micro loop

    double utot=pot.external() + pot.g_internal(spc.p, allpol);
    for (auto &g : pol)
      utot += pot.g_external(spc.p, g);
    sys.checkDrift( utot );

    cout << loop.timing();

  } // end of macro loop

  /*
  iso.test(test);
  gmv.test(test);
  sys.test(test);
  */

  cout << loop.info() << sys.info() << gmv.info() << iso.info() << tit.info()
    << mpol.info();

  rdf.save("rdf_p2p.dat");
  pqr.save("confout.pqr", spc.p);
  top.save("mytopol.top", spc);
  spc.save(ostate);
}