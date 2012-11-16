#ifndef FAUNUS_POINT_H
#define FAUNUS_POINT_H

#ifndef SWIG
#include "faunus/common.h"
#endif

namespace Faunus {
  /*!
   * \brief Cartesian coordinates
   * \author Mikael Lund
   * \date 2002-2007
   */
  class Point {
    public:
      typedef double Tcoord;                       //!< Floating point type for Point coordinates
      Tcoord x,y,z;                                //!< Cartesian coordinates
      Point();                                     //!< Constructor, zero data.
      Point(Tcoord,Tcoord,Tcoord);                 //!< Constructor, set vector
      virtual ~Point();                         
      void clear();                                //!< Zero all data.
      Tcoord len() const;                          //!< Get scalar
      Tcoord dot(const Point&) const;              //!< Angle with another point
      void ranunit(RandomBase&);                   //!< Generate a random unit vector

      virtual void rotate(Geometry::VectorRotate&);//!< Rotate around vector
      virtual void translate(const Geometry::Geometrybase&, const Point&);//!< Translate along a vector
      virtual void scale(const Geometry::Geometrybase&, double);          //!< NPT volume scaling

      Point operator-() const;                     //!< Sign reversal
      const Point operator*(Tcoord) const;         //!< Scale vector
      bool operator==(const Point&) const;         //!< Equality operator
      Point operator+(const Point&) const;         //!< Add two vectors
      Point operator-(const Point&) const;         //!< Subtract vector
      Point& operator+=(const Point&);             //!< Vector addition
      Point& operator*=(const Tcoord);             //!< Scaling vector
      Point& operator<<(std::istream&);            //!< Read from stream
      friend std::ostream &operator<<(std::ostream&, const Point&);//!< Write to stream

      /*
         enum vecFmt {XYZ,XYZQ,XYZQI,ALL};
         virtual void toVector(std::vector<double>&, vecFmt) const; //!< Write data to vector
         virtual int fromVector(const std::vector<double>&, int, vecFmt); //!< Read data from vector
         */
    private:
      inline int anint(Tcoord) const;              //!< Ala fortran function
  };

  inline int Point::anint(Tcoord a) const { return int(a>0 ? a+.5 : a-.5); }

  /*!
   * \brief Class for particles
   * \author Mikael Lund
   * \date 2002-2007
   *
   * Example\n
   * \code
   * std::vector<PointParticle> p(2);
   * p[0].radius = 2.0;
   * p[1].z = 10;
   * std::cout << p[0];
   * \endcode
   */
  class PointParticle : public Point {
    public:
      typedef unsigned char Tid;
      typedef double Tradius;
      typedef double Tcharge;
      typedef float Tmw;
      typedef bool Thydrophobic;

      Tcharge charge;                           //!< Charge number
      Tradius radius;                           //!< Radius
      Tmw mw;                                   //!< Molecular weight
      Tid id;                                   //!< Particle identifier
      Thydrophobic hydrophobic;                 //!< Hydrophobic flag

      PointParticle();
      double volume() const;                    //!< Return volume
      void deactivate();                        //!< Deactivate for use w. faster energy loops
      void clear();                             //!< Zero all data
      PointParticle& operator=(const Point&);   //!< Copy coordinates from a point
      PointParticle& operator=(const AtomData&);//!< Copy data from AtomData
      PointParticle& operator<<(std::istream&); //!< Copy data from stream
      friend std::ostream &operator<<(std::ostream&, const PointParticle&);//!< Write to stream
  };

  /*!
   * \brief Sphero-cylindrical particle
   * \author ...
   * \date ...
   *
   * detailed information here...
   */
  class CigarParticle : public PointParticle {
    public:
      Point dir;
      Point patchdir, patchsides[2], chdir;
      double patchangle, length;

      void rotate(Geometry::VectorRotate&);

      CigarParticle operator+(const Point&) const;
      CigarParticle& operator=(const Point&);
      CigarParticle& operator=(const AtomData&);
      CigarParticle& operator=(const PointParticle&);
      CigarParticle& operator<<(std::istream&);
      friend std::ostream &operator<<(std::ostream &, const CigarParticle&); //!< Output information
  };

#ifdef HYPERSPHERE
  /*!
   * \brief Hypersphere particle
   * \author Martin Trulsson
   * \date Lund, 2009
   * \warning Unfinished - need to transfer from jurassic branch
   */
  class Hyperpoint : public PointParticle {
    public:
      double z1,z2,z3,z4;                     //!< Reduced Coordinates on hypersphere
      translate(const Geometry::Geometrybase&, const Point&);
      friend std::ostream &operator<<(std::ostream&, Hyperpoint&);
      Hyperpoint &operator<<(std::istream&);

      void clear() {
        z1=z2=z3=0;
        z4=1;
      }

      Hyperpoint() { clear(); }

      /*!
       * \brief Squared distance between two points.
       * \return \f[ r^2 = z_1z_1' + z_2z_2' + z_3z_3' + z_4z_4' \f]
       */
      inline double sqdist(const Hyperpoint &p) const {
        return z1*p.z1+z2*p.z2+z3*p.z3+z4*p.z4;
      }

      /*!
       * \brief Geodesic distance between two hyperpoints
       * \return \f[ r_{\mbox{\scriptsize{geod}}} = \arccos{ (r^2) } \f]
       */
      inline double geodesic(const Hyperpoint &p) const {
        return std::acos(sqdist(p));
      }
  };
#endif

}//namespace
#endif
