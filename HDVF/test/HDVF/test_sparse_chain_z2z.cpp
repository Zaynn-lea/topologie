#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/OSM/OSM.h>
#include <CGAL/Z2.h>
#include <CGAL/OSM/Sparse_chain_Z2Z.h>


typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::OSM::Sparse_chain_z2<CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain_z2<CGAL::OSM::ROW> Row_chain ;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::COLUMN, CGAL::OSM::Sparse_chain_z2_core> Column_matrix;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::ROW, CGAL::OSM::Sparse_chain_z2_core> Row_matrix;


/*
// temporarly :
typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::OSM::Sparse_chain_z2_core<Coefficient_ring, CGAL::OSM::COLUMN> Column_chain;
typedef CGAL::OSM::Sparse_chain_z2_core<Coefficient_ring, CGAL::OSM::ROW> Row_chain ;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::COLUMN, CGAL::OSM::Sparse_chain_z2_core> Column_matrix;
typedef CGAL::OSM::Sparse_matrix<Coefficient_ring, CGAL::OSM::ROW, CGAL::OSM::Sparse_chain_z2_core> Row_matrix;
*/

int main(int argc, char **argv)
{
    CGAL::Z2 one(1), zero(0);
    Column_chain c1(5), c2(5), c3(5);
    Row_chain r1(5), r2(5), r3(5) ;

    c1.set_coefficient(0, one);
    c1.set_coefficient(2, one);
    c2.set_coefficient(0, one);
    c2.set_coefficient(1, one);
    c3.set_coefficient(2, one);
    c3.set_coefficient(0, one);

    r1.set_coefficient(0, one);
    r1.set_coefficient(2, one);
    r2.set_coefficient(0, one);
    r2.set_coefficient(1, one);
    r3.set_coefficient(2, one);
    r3.set_coefficient(0, one);

    std::cerr << "-- Test Sparse_chain_z2z operator==" << std::endl;

    {
        std::cerr << "----> Column chains" << std::endl;
        bool comp_true(c1==c3) ;
        std::cerr << "Compare similar column chains: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(c1==c2) ;
        std::cerr << "Compare different column chains: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        bool comp_true(r1==r3) ;
        std::cerr << "Compare similar row chains: " << comp_true << std::endl ;
        assert(comp_true) ;
        bool comp_false(r1==r2) ;
        std::cerr << "Compare different row chains: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }

    std::cerr << "-- Test Sparse_chain_z2z operator+" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5), res2(4);
        res.set_coefficient(2, one);
        res.set_coefficient(1, one);
        bool comp_true(c1+c2==res) ;
        std::cerr << "Test c1+c2: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(2, one);
        res2.set_coefficient(1, one);
        bool comp_false(c1+c2==res2) ;
        std::cerr << "Compare c1+c2 with invalid size chain: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5), res2(4);
        res.set_coefficient(2, one);
        res.set_coefficient(1, one);
        bool comp_true(r1+r2==res) ;
        std::cerr << "Test r1+r2: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(2, one);
        res2.set_coefficient(1, one);
        bool comp_false(r1+r2==res2) ;
        std::cerr << "Compare r1+r2 with invalid size chain: " << comp_false << std::endl ;
        assert(!comp_false) ;
    }

    std::cerr << "-- Test Sparse_chain_z2z operator-" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5);
//        res.set_coefficient(0, one);
        res.set_coefficient(2, one);
        res.set_coefficient(1, one);
        bool comp_true(c1-c2==res) ;
        std::cerr << "Test c1-c2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5);
 //       res.set_coefficient(0, one);
        res.set_coefficient(2, one);
        res.set_coefficient(1, one);
        bool comp_true(r1-r2==res) ;
        std::cerr << "Test r1-r2: " << comp_true << std::endl ;
        assert(comp_true) ;
    }

    std::cerr << "-- Test Sparse_chain_z2z operator* (coefficient*chain)" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(c1);
        bool comp_true(one*c1==res) ;
        std::cerr << "Test 1*c1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5);
        res.set_coefficient(0, one);
        res.set_coefficient(2, one);
        bool comp_true(one*r1==res) ;
        std::cerr << "Test 3*r1: " << comp_true << std::endl ;
        assert(comp_true) ;
    }

    std::cerr << "-- Test Sparse_chain_z2z operator* (chain*coefficient)" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(5);
        res.set_coefficient(0, one);
        res.set_coefficient(2, one);
        bool comp_true(c1*one==res) ;
        std::cerr << "Test c1*3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }
    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(5);
        res.set_coefficient(0, one);
        res.set_coefficient(2, one);
        bool comp_true(r1*one==res) ;
        std::cerr << "Test r1*3: " << comp_true << std::endl ;
        assert(comp_true) ;
    }

    std::cerr << "-- Test Sparse_chain_z2z operator/" << std::endl;
    {
        std::cerr << "----> Column chains" << std::endl;
        Column_chain res(c1), res2(5);
        res.set_coefficient(1, one);
        bool comp_true(c1==(res/1)) ;
        std::cerr << "Test c1 == res/1: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(0, one);
        bool comp_true2(res2==(res/std::vector<size_t>({1,2}))) ;
        std::cerr << "Test res2 == res/{1,2}: " << comp_true2 << std::endl ;
        assert(comp_true2) ;
    }

    {
        std::cerr << "----> Row chains" << std::endl;
        Row_chain res(r1), res2(5);
        res.set_coefficient(1, one);
        bool comp_true(r1==(res/1)) ;
        std::cerr << "Test r1 == res/1: " << comp_true << std::endl ;
        assert(comp_true) ;
        res2.set_coefficient(0, one);
        bool comp_true2(res2==(res/std::vector<size_t>({1,2}))) ;
        std::cerr << "Test res2 == res/{1,2}: " << comp_true2 << std::endl ;
        assert(comp_true2) ;
    }

    {
        std::cerr << "----> Read/write chains" << std::endl;
        Column_chain cc1;
        Row_chain rr1;

        // Column
        CGAL::OSM::write_chain(c1, "tmp/column_chain.osm");
        CGAL::OSM::read_chain(cc1, "tmp/column_chain.osm");
        bool comp_true(c1==cc1) ;
        std::cerr << "Test write/read for column chain: " << comp_true << std::endl ;
        assert(comp_true) ;
        // Row
        CGAL::OSM::write_chain(r1, "tmp/row_chain.osm");
        CGAL::OSM::read_chain(rr1, "tmp/row_chain.osm");
        bool comp_true2(r1==rr1) ;
        std::cerr << "Test write/read for row chain: " << comp_true2 << std::endl ;
        assert(comp_true2) ;

    }



    return 0;
}

