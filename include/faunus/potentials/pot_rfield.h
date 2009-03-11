#ifndef FAU_POT_RFIELD_H
#define FAU_POT_RFIELD_H

#include "faunus/point.h"
#include "faunus/legendre.h"
#include "faunus/potentials/base.h"

namespace Faunus {
  /*! \brief Charge interactions in- and outside a spherical cavity
   *  \author mikaek lund
   *  \date 2005-2006 Canberra
   *  \todo The sketch could be more clear.
   *  \note Woodward+Svensson, J.Phys.Chem, 1991, 95, p7471.
   *  \image html reactionfield.png 
   *
   *  Calculate the pair interaction between charges locate in-
   *  and/or outside a spherical dielectric discontinuity.
   *  
   *  Faunus::inputfile keywords:
   *  \li "rfield_epso" - Dielectric constant outside sphere
   *  \li "rfield_epsi" - Dielectric constant inside sphere
   *  \li "rfield_cavity" - Radius of sphere (origo = [0,0,0]).
   *  \li "rfield_bjerrum" - the Bjerrum length OUTSIDE the sphere
   */
  class pot_rfield {
    private:
      legendre l;           // legendre polynomium 
      int steps;            // # of steps in summation (precision)
    public:
      double eo, ei;        // dielectric constants (solvent, inside sphere)
      double f;             //!< Factor to convert to kT
      double a;             // cavity radius
      pot_rfield( inputfile &in ) {
        a  = in.getflt("rfield_cavity");
        eo = in.getflt("rfield_epso",80.);
        ei = in.getflt("rfield_epsi",1.);
        f  = in.getflt("rfield_bjerrum",7.0025) * eo;
        steps=in.getint("rfield_steps",50); 
        l.resize(steps);
      };
      inline double pairpot(const particle &p1, const particle &p2) {
        return (p2.charge!=0) ? p2.charge * phi(p1, p2) : 0;
      }

      /*! \brief Self energy (interaction with image charge)
       *  \return \f$ \frac{1}{2} q \phi \f$
       */
      inline double selfenergy(const particle &p) {
        return 0.5 * p.charge * phi(p,p,true);
      }

      /*! \brief Calculate potential in point "p" from charge in "p0".
       *  \note Multiply w. lB*eo (=f) to get kT.
       *  \param p0 Charged particle
       *  \param p  Calculate potential in this point
       *  \param self Set to true to calculate self-energy of p0
       */
      double phi(const particle &p0, const point &p, bool self=false) {
        if (p0.charge==0)
          return 0;
        double sum = 0;
        double r   = p.len();
        double r0  = p0.len();
        double dn;
        double costheta = p.dot(p0) / (r * r0);
        l.eval(costheta); 

        //charge outside, potential outside
        if (r0>a && r>a) {
          for (int n=1; n<steps; n+=1) {
            dn   = double(n);
            sum += pow(a*a/(r*r0),n+1) * l.p[n]
              / (  (ei+eo*(1+1/dn))   );
          }
          sum *= (eo-ei) / (eo * a );
          if (self==false)
            sum += 1/(eo*p.dist(p0));
        }
        //charge inside, potential inside
        if (r0<a && r<a) {
          for (int n=0; n<steps; n+=1) {
            dn = double(n);
            sum += pow(r*r0/(a*a),n) * l.p[n]
              / ( eo-ei*(dn/(dn+1)) )  ;
          }
          sum *= (ei-eo) / (ei*a);
          if (self==false)
            sum += 1/(ei*p.dist(p0));
        }

        //charge inside, potential outside
        if (r0<a && r>a) {
          for (int n=0; n<steps; n+=1) {
            dn=double(n);
            sum += (2*dn+1) * pow(r0/r,n) * l.p[n]
              / (  ( dn*ei + eo*(dn+1) )  );
          };
          sum = sum / r;
        };

        //charge outside, potential inside (as above...)
        if (r0>a && r<a) {
          for (int n=1; n<steps; n+=1) {
            dn = double(n);
            sum += pow(r/r0,n) * l.p[n]
              / ( (ei+eo*(1+1/dn))   ) ;
          };
          sum *= (eo-ei) / (eo * r0) ;
          sum += 1/(eo*p.dist(p0));
        };
        return sum * p0.charge;
      }

      string info() {
        std::ostringstream o;
        o << endl;
        o << "#  Reaction field pair-potential:" << endl
          << "#    Cavity radius                = " << a << endl
          << "#    Dielectric constant (in out) = " << ei << " " << eo << endl
          << "#    Legendre steps               = " << steps << endl;
        return o.str();
      }
  };
}//namespace
#endif