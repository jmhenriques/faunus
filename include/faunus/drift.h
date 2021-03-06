#ifndef FAUNUS_DRIFT_H
#define FAUNUS_DRIFT_H

#include <faunus/common.h>
#include <faunus/average.h>

namespace Faunus {
  class EnergyDrift {
    private:
      double delta;
      double initial;
      Average<double> avg;
    public:
      double drift;
      EnergyDrift();
      void init(const double&);
      double current() const;
      EnergyDrift& operator+=(const double&);
      EnergyDrift& operator()(const std::pair<double,double>&);
      double weight;
      bool rejection;
      double checkDrift(const double&);
      string info();
      void test(UnitTest&);
  };

  EnergyDrift::EnergyDrift() {
    delta=initial=drift=0;
  }

  void EnergyDrift::init(const double &u0) {
    avg.reset();
    delta=drift=0;
    initial=u0;
    avg+=u0;
  }

  double EnergyDrift::current() const {
    return initial + delta;
  }

  EnergyDrift& EnergyDrift::operator+=(const double &du) {
    delta+=du;
    avg+=current();
    return *this;
  }

  /**
   * @brief Compute acceptance probability of a move 
   *
   * @details Useful when sampling with waste-recycling method.                                              
   * @param du Energy change pair. First member =du in case of acceptance and
   * =0 in case of rejection. Second member =du in both cases.  
   *
   * [More info](http://dx.doi.org/10.1007/3-540-35273-2_4)
   */
  EnergyDrift& EnergyDrift::operator()(const std::pair<double,double> &du) {
    delta+=du.first;
    avg+=current();
    weight=std::exp(-du.second);
    if (du.first==du.second) rejection=false;
    else rejection=true;
    return *this;
  }

  double EnergyDrift::checkDrift(const double &snapshot) {
    drift = snapshot-current();
    return drift;
  }

  string EnergyDrift::info() {
    using namespace Faunus::textio;
    std::ostringstream o;
    o << header("System Energy and Drift");
    if (avg.cnt>0) {
      char w=25;
      o << textio::pad(SUB,w, "Average") << avg.avg() << kT << ", "
        << sigma << "=" << avg.stdev() << endl
        << textio::pad(SUB,w, "Initial energy") << initial << kT << endl
        << textio::pad(SUB,w, "Initial + changes") << current() << kT << endl;
      o.precision(4);
      o << pad(SUB,w, "Total energy drift") << drift << kT
        << " (" << drift/current()*100.
        << percent << ")" << endl;
    }
    return o.str();
  }

  void EnergyDrift::test(UnitTest &t) {
    //t("initialEnergy", initial, 1e-3);
    t("energyAverage", avg.avg() );
    t("relativeEnergyDrift", std::abs(drift/current()), 10.0 ); // allow 200% deviation    
  }
}//namespace
#endif
