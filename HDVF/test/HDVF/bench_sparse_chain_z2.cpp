#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/OSM/OSM.h>
#include <CGAL/Z2.h>
#include <CGAL/Timer.h>
#include <CGAL/OSM/Sparse_chain_Z2Z.h>
#include <CGAL/OSM/Sparse_chain.h>
#include "../../examples/HDVF/tools/distrib.h"

typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::OSM::Sparse_chain_z2<CGAL::OSM::COLUMN> Column_chain_z2;
typedef CGAL::OSM::Sparse_chain_z2<CGAL::OSM::ROW> Row_chain_z2 ;
typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::ROW> Row_chain ;

typedef CGAL::Timer Timer;

struct stat_data {
    std::vector<size_t> labels;
    std::vector<double> tcreate, tcreate_z2, diff_create;
    std::vector<double> tsum, tsum_z2, diff_sum;
    std::vector<double> tget, tget_z2, diff_get;
};

std::ostream& operator<< (std::ostream& out, const stat_data& stats) {
    out << "labels = [";
    for (size_t label : stats.labels)
        out << label << " ";
    out << "];" << std::endl;

    out << "tcreate = [";
    for (double t : stats.tcreate)
        out << t << " ";
    out << "];" << std::endl;

    out << "tcreate_z2 = [";
    for (double t : stats.tcreate_z2)
        out << t << " ";
    out << "];" << std::endl;

    out << "tsum = [";
    for (double t : stats.tsum)
        out << t << " ";
    out << "];" << std::endl;

    out << "tsum_z2 = [";
    for (double t : stats.tsum_z2)
        out << t << " ";
    out << "];" << std::endl;

    out << "tget = [";
    for (double t : stats.tget)
        out << t << " ";
    out << "];" << std::endl;

    out << "tget_z2 = [";
    for (double t : stats.tget_z2)
        out << t << " ";
    out << "];" << std::endl;
    return out;
}

std::ostream& operator<< (std::ostream& out, const std::unordered_set<size_t>& chain) {
    out << "[";
    for (size_t i : chain)
        out << i << " ";
    out << "];" << std::endl;
    return out;
}

// Benchmark
// -> create chains of size N with M% coefs
// -> sum chains of size N with M% coefs
// -> get coef in chain of size N with M% coefs

std::vector<std::unordered_set<size_t> > rand_data (int nb_chains, int N, int M, std::mt19937 &rng)
{
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,N-1);
    std::vector<std::unordered_set<size_t> > res(nb_chains);
    const int N_vals(int(double(N)*M/100));

    for (int i=0; i<nb_chains; ++i) {
        while (res.at(i).size() < N_vals) {
            int j = dist(rng) ; // Random index
            res.at(i).insert(j);
        }
    }
    return res ;
}

std::pair<std::vector<Column_chain>, std::vector<Column_chain_z2> > generate_chains (const std::vector<std::unordered_set<size_t> >& data, int N, Timer& timer, Timer& timer_z2) {
    int nb_chains(data.size());
    std::vector<Column_chain_z2> chains_z2(nb_chains);
    std::vector<Column_chain> chains(nb_chains);

    // Generate standard chains
    timer.start();
    for (int i=0; i<nb_chains; ++i) {
        chains.at(i) = Column_chain(N);
        for (std::unordered_set<size_t>::iterator it = data.at(i).begin(); it != data.at(i).end(); ++it) {
            chains.at(i).set_coefficient(*it,1);
        }
    }
    timer.stop();

    // Generate z2 chains
    timer_z2.start();
    for (int i=0; i<nb_chains; ++i) {
        chains_z2.at(i) = Column_chain_z2(N);
        for (std::unordered_set<size_t>::iterator it = data.at(i).begin(); it != data.at(i).end(); ++it) {
            chains_z2.at(i).set_coefficient(*it,1);
        }
    }
    timer_z2.stop();

    // Compute the distribution of times

    distrib<double> d(50), d_z2(50);
    std::vector<Column_chain_z2> chains_tmp_z2(nb_chains);
    std::vector<Column_chain> chains_tmp(nb_chains);
    Timer t, t_z2;

    // Generate standard chains
    for (int i=0; i<nb_chains; ++i) {
        t.start();
        chains_tmp.at(i) = Column_chain(N);
        for (std::unordered_set<size_t>::iterator it = data.at(i).begin(); it != data.at(i).end(); ++it) {
            chains_tmp.at(i).set_coefficient(*it,1);
        }
        t.stop();
        d.add_data(t.time());
    }

    // Generate z2 chains
    for (int i=0; i<nb_chains; ++i) {
        t_z2.start();
        chains_tmp_z2.at(i) = Column_chain_z2(N);
        for (std::unordered_set<size_t>::iterator it = data.at(i).begin(); it != data.at(i).end(); ++it) {
            chains_tmp_z2.at(i).set_coefficient(*it,1);
        }
        t_z2.stop();
        d_z2.add_data(t_z2.time());
    }

    std::cerr << "---" << "Create Sparse_chains: " <<  d ;
    std::cerr << "Create Sparse_chains_z2: " <<  d_z2 ;

//    size_t id_max, id_min;
//    d.get_min(true,&id_min);
//    d.get_max(true, &id_max);
//    std::cerr << "       " << "min for: " << data.at(id_min) ;
//    std::cerr << "       " << "max for: " << data.at(id_max) ;

    return std::pair<std::vector<Column_chain>, std::vector<Column_chain_z2> >({chains, chains_z2});
}

void sum_chains(std::pair<const std::vector<Column_chain>, const std::vector<Column_chain_z2>> data, Timer& timer, Timer& timer_z2) {
    std::vector<Column_chain> chains(data.first);
    std::vector<Column_chain_z2> chains_z2(data.second);
    size_t nb_chains(chains.size());

    // Sum standard chains
    timer.start();
    for (int i=0; i<nb_chains-1; ++i) {
        chains.at(i) += chains.at(i+1);
    }
    timer.stop();

    // Sum z2 chains
    timer_z2.start();
    for (int i=0; i<nb_chains-1; ++i) {
        chains_z2.at(i) += chains_z2.at(i+1);
    }
    timer_z2.stop();

    // Compute the distribution of times
    std::vector<Column_chain> chains_tmp(data.first);
    std::vector<Column_chain_z2> chains_tmp_z2(data.second);

    distrib<double> d, d_z2;
    Timer t, t_z2;

    // Sum standard chains
    for (int i=0; i<nb_chains-1; ++i) {
        t.start();
        chains_tmp.at(i) += chains_tmp.at(i+1);
        t.stop();
        d.add_data(t.time());
    }

    // Sum z2 chains
    for (int i=0; i<nb_chains-1; ++i) {
        t_z2.start();
        chains_tmp_z2.at(i) += chains_tmp_z2.at(i+1);
        t_z2.stop();
        d_z2.add_data(t_z2.time());
    }

    std::cerr << "---" << "Sum Sparse_chains: " <<  d ;
    std::cerr << "Sum Sparse_chains_z2: " <<  d_z2 ;
}

void get_coef_chains(std::pair<const std::vector<Column_chain>, const std::vector<Column_chain_z2>> data, int N, Timer& timer, Timer& timer_z2, std::mt19937 &rng) {
    std::vector<Column_chain> chains(data.first);
    std::vector<Column_chain_z2> chains_z2(data.second);
    size_t n_chain(chains.at(0).dimension());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,n_chain-1);

    // Indices to get
    std::vector<size_t> indices(N);
    for (int i=0; i<N; ++i)
        indices.at(i) = dist(rng);

    size_t nb_chains(chains.size());

    // Get coef in standard chains
    int tmp;
    timer.start();
    for (int i=0; i<nb_chains-1; ++i) {
        for (int j=0; j<N; ++j) {
            tmp = chains.at(i).get_coefficient(j);
        }
    }
    timer.stop();

    // Get coef in z2 chains
    Coefficient_ring tmp2;
    timer_z2.start();
    for (int i=0; i<nb_chains-1; ++i) {
        for (int j=0; j<N; ++j) {
            tmp2 = chains_z2.at(i).get_coefficient(j);
        }
    }
    timer_z2.stop();

    // Compute the distribution of times
    std::vector<Column_chain> chains_tmp(data.first);
    std::vector<Column_chain_z2> chains_tmp_z2(data.second);
    distrib<double> d, d_z2;
    Timer t, t_z2;

    // Get coef in standard chains
    for (int i=0; i<nb_chains-1; ++i) {
        t.start();
        for (int j=0; j<N; ++j) {
            tmp = chains.at(i).get_coefficient(j);
        }
        t.stop();
        d.add_data(t.time());
    }

    // Get coef in z2 chains
    for (int i=0; i<nb_chains-1; ++i) {
        t_z2.start();
        for (int j=0; j<N; ++j) {
            tmp2 = chains_z2.at(i).get_coefficient(j);
        }
        t_z2.stop();
        d_z2.add_data(t_z2.time());
    }
    std::cerr << "---" << "Get Sparse_chains: " <<  d ;
    std::cerr << "Get Sparse_chains_z2: " <<  d_z2 ;
}

void  fct(int nb_chains, int N, int M, std::mt19937 &rng, stat_data& stats)
{
    Timer tcreate, tsum, tget;
    Timer tcreate_z2, tsum_z2, tget_z2;

    std::cerr << "----------------------" << std::endl;
    std::cerr << "Chains of size " << N << " with " << M << "% coefs" << std::endl;

    // Generate random data for nb_chains chains of size N with M% coefficients
//    std::cerr << "==> Data generation" << std::endl;
    std::vector<std::unordered_set<size_t> > data(rand_data(nb_chains, N, M, rng));
//    std::cerr << "<== END Data generation" << std::endl;

    // Test chains generation
//    std::cerr << "==> Chains generation" << std::endl;
    std::pair<std::vector<Column_chain>, std::vector<Column_chain_z2> > p(generate_chains(data, N, tcreate, tcreate_z2));
//    std::cerr << "<== END Chains generation" << std::endl;

    // Test sum
//    std::cerr << "==> Chains sums" << std::endl;
    sum_chains(p, tsum, tsum_z2);
//    std::cerr << "<== END Chains sums" << std::endl;

    // Test get
//    std::cerr << "==> Chains get" << std::endl;
    get_coef_chains(p, 100, tget, tget_z2, rng);
//    std::cerr << "<== END Chains get" << std::endl;


    stats.labels.push_back(N);

    stats.tcreate.push_back(tcreate.time());
    stats.tcreate_z2.push_back(tcreate_z2.time());
    stats.diff_create.push_back(-(tcreate.time()-tcreate_z2.time())/tcreate.time()*100);

    stats.tsum.push_back(tsum.time());
    stats.tsum_z2.push_back(tsum_z2.time());
    stats.diff_sum.push_back(-(tsum.time()-tsum_z2.time())/tsum.time()*100);

    stats.tget.push_back(tget.time());
    stats.tget_z2.push_back(tget_z2.time());
    stats.diff_get.push_back(-(tget.time()-tget_z2.time())/tget.time()*100);
}

int main(int , char* argv[])
{
    std::random_device dev;
    //    std::mt19937 rng(dev()); // Random
    std::mt19937 rng(0); // Pseudo-random (for tests)

    stat_data stats;

    fct(10000, 100, 10, rng, stats);
    fct(10000, 500, 10, rng, stats);
    fct(10000, 1000, 10, rng, stats);
    fct(10000, 5000, 10, rng, stats);
    fct(10000, 10000, 10, rng, stats);

    std::cerr << "----------------------" << std::endl;
    std::cerr << stats;
    return 0;
}
