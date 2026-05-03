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

typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::OSM::Sparse_chain_z2<CGAL::OSM::COLUMN> Column_chain_z2;
typedef CGAL::OSM::Sparse_chain_z2<CGAL::OSM::ROW> Row_chain_z2 ;
typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain<int, CGAL::OSM::ROW> Row_chain ;

typedef CGAL::Timer Timer;

// Benchmark
// -> create chains of size N with M% coefs
// -> sum chains of size N with M% coefs
// -> get coef in chain of size N with M% coefs

std::vector<std::set<size_t> > rand_data (int nb_chains, int N, int M, std::mt19937 &rng)
{
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,N-1);
    std::vector<std::set<size_t> > res(nb_chains);
    const int N_vals(int(double(N)*M/100));

    for (int i=0; i<nb_chains; ++i) {
        int nb_vals = 0 ;
        while (nb_vals < N_vals)
        {
            int j = dist(rng) ; // Random column
//            if(res.at(i).find(j) != res.at(i).end()) // j not in the ith-chain
            {
                res.at(i).insert(j);
                ++nb_vals;
            }
        }
    }
    return res ;
}

std::pair<std::vector<Column_chain>, std::vector<Column_chain_z2> > generate_chains(const std::vector<std::set<size_t> >& data, int N, Timer& timer, Timer& timer_z2) {
    int nb_chains(data.size());
    std::vector<Column_chain_z2> chains_z2(nb_chains);
    std::vector<Column_chain> chains(nb_chains);

    // Generate standard chains
    timer.start();
    for (int i=0; i<nb_chains; ++i) {
        chains.at(i) = Column_chain(N);
        for (std::set<size_t>::iterator it = data.at(i).begin(); it != data.at(i).end(); ++it) {
            chains.at(i).set_coefficient(*it,1);
        }
    }
    timer.stop();

    // Generate z2 chains
    timer_z2.start();
    for (int i=0; i<nb_chains; ++i) {
        chains_z2.at(i) = Column_chain_z2(N);
        for (std::set<size_t>::iterator it = data.at(i).begin(); it != data.at(i).end(); ++it) {
            chains_z2.at(i).set_coefficient(*it,1);
        }
    }
    timer_z2.stop();
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
}

void  fct(int nb_chains, int N, int M, std::mt19937 &rng)
{
    Timer tcreate, tsum, tget;
    Timer tcreate_z2, tsum_z2, tget_z2;

    // Generate random data for nb_chains chains of size N with M% coefficients
//    std::cerr << "==> Data generation" << std::endl;
    std::vector<std::set<size_t> > data(rand_data(nb_chains, N, M, rng));
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

    std::cerr << nb_chains << " chains created of size " << N << " with " << M << "% coefs" << std::endl;
    std::cerr << "Sparse_chain\t: " << tcreate.time() << " sec." << std::endl;
    std::cerr << "Sparse_chain_z2\t: " << tcreate_z2.time() << " sec." << std::endl;
    std::cerr << "Diff:\t\t" << -(tcreate.time()-tcreate_z2.time())/tcreate.time()*100 << " %" << std::endl;

    std::cerr << nb_chains << " chains sumed of size " << N << " with " << M << "% coefs" << std::endl;
    std::cerr << "Sparse_chain\t: " << tsum.time() << " sec." << std::endl;
    std::cerr << "Sparse_chain_z2\t: " << tsum_z2.time() << " sec." << std::endl;
    std::cerr << "Diff:\t\t" << -(tsum.time()-tsum_z2.time())/tsum.time()*100 << " %" << std::endl;

    std::cerr << nb_chains << " chains of size " << N << " with " << M << "% coefs:  " << 100 << " gets" << std::endl;
    std::cerr << "Sparse_chain\t: " << tget.time() << " sec." << std::endl;
    std::cerr << "Sparse_chain_z2\t: " << tget_z2.time() << " sec." << std::endl;
    std::cerr << "Diff:\t\t" << -(tget.time()-tget_z2.time())/tget.time()*100 << " %" << std::endl;
}

int main(int , char* argv[])
{
    std::random_device dev;
    //    std::mt19937 rng(dev()); // Random
    std::mt19937 rng(0); // Pseudo-random (for tests)

    fct(10000, 100, 10, rng);
    fct(10000, 500, 10, rng);
    fct(10000, 1000, 10, rng);
    fct(10000, 5000, 10, rng);
    fct(10000, 10000, 10, rng);
    return 0;
}
