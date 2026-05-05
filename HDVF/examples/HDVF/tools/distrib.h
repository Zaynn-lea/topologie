//
//  distrib.hpp
//  test1
//
//  Created by Bac Alexandra on 21/03/2022.
//

#ifndef distrib_hpp
#define distrib_hpp

#include <iostream>
using namespace std ;

template <typename T>
class distrib {
    vector<T> _data ;
    T _mean_sum, _min, _max, _stdev ;
    size_t _n_bins;
    std::vector<size_t> _hist;
public:
    distrib(size_t nbins = 10) : _mean_sum(0.), _min(0.), _max(0.), _stdev(0.), _n_bins(nbins), _hist(nbins)  { } ;
    distrib(const distrib<T>& d) : _data(d._data), _mean_sum(d._mean_sum), _min(d._min), _max(d._max), _stdev(d._stdev), _n_bins(d._n_bins), _hist(d._hist) {}
    distrib(vector<T> &vec, size_t nbins = 10) : _data(vec), _mean_sum(0.), _min(0.), _max(0.), _stdev(0.), _n_bins(nbins), _hist(nbins)
    {
        if (_data.size()>0)
        {
            _mean_sum = _data.at(0) ;
            _min = _data.at(0) ;
            _max = _data.at(0) ;
            for (int i=1; i<_data.size(); ++i)
            {
                _mean_sum += _data.at(i) ;
                if (_data.at(i) < _min)
                    _min = _data.at(i) ;
                if (_data.at(i) > _max)
                    _max = _data.at(i) ;
            }
        }
        update_stdev();
        update_hist();
    } ;
    
    template <typename TT>
    friend std::ostream& operator<<(std::ostream& out, distrib<TT>& dist);
    
    void update_stdev ()
    {
        T mean = _mean_sum / _data.size(), tmp ;
        for (int i=0; i<_data.size(); ++i)
        {
            tmp = _data.at(i) - mean ;
            _stdev += tmp*tmp ;
        }
        _stdev = sqrt(_stdev / _data.size()) ;
    } ;

    void update_hist () {
        double delta(double(_max-_min)/_n_bins);
        for (T x : _data) {
            size_t i;
            if (x == _max)
                i = _n_bins-1;
            else
                i = floor((x-_min)/delta);
            ++_hist.at(i);
        }
    }
    // Accesseurs
    T get_mean () const { return _mean_sum/_data.size() ; } ;
    T get_stdev () { update_stdev() ; return _stdev ; } ;
    T get_min (bool with_data=false, size_t* pindex=NULL) const {
        if (with_data) {
            typename std::vector<T>::const_iterator it = std::find(_data.cbegin(), _data.cend(), _min);
            *pindex = it-_data.cbegin();
        }
        return _min ;
    }

    T get_max (bool with_data=false, size_t* pindex=NULL) const {
        if (with_data) {
            typename std::vector<T>::const_iterator it = std::find(_data.cbegin(), _data.cend(), _max);
            *pindex = it-_data.cbegin();
        }
        return _max ;
    }

    void add_data (T dat)
    {
        if (_data.size() == 0)
        {
            _min = dat ;
            _max = dat ;
        }
        else
        {
            if (dat < _min)
                _min = dat ;
            if (dat > _max)
                _max = dat ;
        }
        _data.push_back(dat) ; _mean_sum += dat ;
    }
    void clear ()
    { _data.clear() ; _min = 0. ; _max = 0. ; _mean_sum = 0. ; _stdev = 0. ;}
};

template <typename T>
std::ostream& operator<<(std::ostream& out, distrib<T>& dist) {
    out << "K min : " << dist.get_min() << " - max : " << dist.get_max() << std::endl ;
    out << "K mean : " << dist.get_mean() << " - sigma : " << dist.get_stdev() << std::endl;

    dist.update_hist();
    out << "edges = [";
    for (int i = 0; i<=dist._n_bins; ++i)
        out << dist._min + i*(double(dist._max-dist._min)/dist._n_bins) << " " ;
    out << "];" << std::endl;

    out << "hist = [";
    for (T x : dist._hist)
        out << x << " " ;
    out << "];" << std::endl;
    
    return out;
}

#endif /* distrib_hpp */
