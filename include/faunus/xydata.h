#ifndef FAU_XYDATE_h
#define FAU_XYDATA_h
#include "faunus/common.h"
#include "faunus/average.h"

namespace Faunus {
  /*!
   * \brief Class to handle a dataset of xy values with equidistant x-values.
   * \author Mikael Lund
   * \date March 2007
   */
  template <class T=double>
    class xydata {
      private:
        T invres;
        vector<T> y;                  //!< y value vector
        unsigned short x2i(T);        //!< xvalue->index
        T i2x(int);                   //!< index->xvalue
        T i2y(int);                   //!< index->yvalue
      public:
        string comment;               //!< Arbitrary comment
        T xmin;                       //!< Minimum x-value
        T xmax;                       //!< Maximum x-value
        T res;                        //!< Data resolution
        void setup(T,T,T);            //!< Reset boundaries
        xydata(T=0.,T=0.,T=0.1);      //!< Constructor
        void add(T,T);                //!< Add a point
        bool add(vector<T> &, vector<T> &); //!< Add a high res. dataset
        inline T x2y(T);              //!< xvalue->yvalue
        void list();                  //!< Print info
    };

  template <class T>
    unsigned short xydata<T>::x2i(T xin) {
      return static_cast<unsigned short>((xin-xmin)*invres+.5);
    }
  template <class T> T xydata<T>::i2x(int i) { return i*res+xmin; }
  template <class T> T xydata<T>::i2y(int i) { return y[i]; }
  template <class T> T xydata<T>::x2y(T x) { return y[ x2i(x) ]; }

  // re-set xmin, xmax and data resolution
  template <class T>
    void xydata<T>::setup(T min, T max, T resolution) {
      xmin=min; //minimum xvalue
      xmax=max; //maximum xvalue
      res=resolution; //xsteps
      invres=1./res;   //for speed
      y.resize( int( (xmax-xmin)/res+0.5) ); //adjust y vector
    }

  /*!
   *  \param min Minimum x-value
   *  \param max Maximum x-value
   *  \param resolution Data resolution
   */
  template <class T>
    xydata<T>::xydata(T min, T max, T resolution) {
      setup(min,max,resolution);
    }

  template <class T>
    void xydata<T>::add(T xin, T yin) {
      unsigned int i=x2i(xin);
      if (i>=y.size())
        y.resize(i+1);
      y[i] = yin;
    }

  /*! Load a high resolution dataset and average it
   * to the current resolution.
   */
  template <class T>
    bool xydata<T>::add(vector<T> &xin, vector<T> &yin) {
      T d0,d;
      unsigned int i=0, n=xin.size()-1;
      average<T> meanx, meany;
      setup( xin[0], xin[n], res); 
      d0=std::abs( xin[1] - xin[0] );
      if (d0>res)
        return false;
      d=d0;
      while (i<n) {
        while (d<res && i<n) {
          meanx+=xin.at(i);              // Add data to average
          meany+=yin.at(i);
          //cout << xin[i] << endl;
          i++; 
          d+=std::abs(xin.at(i)-xin.at(i-1));
        }
        add( meanx.avg(), meany.avg() );
        //cout << meanx.avg() << ", " << meany.avg() << " " <<i2x(x2i(meanx.avg())) << endl;
        meanx.reset();
        meany.reset();
        d=0;
      }
      xmax=i2x( y.size()-1 ) ;
      return true;
    }

  template <class T>
    void xydata<T>::list() {
      for (unsigned int i=0; i<y.size(); i++)
        std::cout << i2x(i) << " " << i2y(i) << std::endl;
    }
}
#endif