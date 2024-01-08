#include "uniform_radix_sort.h"
#include <random>
#include <iostream>
#include <chrono>

template <class Executable>
auto timeit(Executable exe) -> size_t{

    using namespace std::chrono;
    auto s  = high_resolution_clock::now();
    exe();
    auto l  = duration_cast<milliseconds>(high_resolution_clock::now() - s).count();

    return l;
}

int main(){

    constexpr auto SZ   = size_t{1} << 15;
    auto sz_rand_dev    = std::bind(std::uniform_int_distribution<size_t>(0u, SZ), std::mt19937{});

    while (true){
        auto data_rand_dev  = std::bind(std::uniform_int_distribution<size_t>{}, std::mt19937{std::chrono::high_resolution_clock::now().time_since_epoch().count()});
        auto sz     = sz_rand_dev();        
        auto test   = std::vector<size_t>(sz);
        std::generate(test.begin(), test.end(), data_rand_dev);
        auto test2  = test;

        std::sort(test2.begin(), test2.end());
        dg::uniform_radix_sort::radix_sort(test.data(), test.data() + test.size());

        if (test != test2){
            std::cout << "mayday" << std::endl;
        }
    }
    

}