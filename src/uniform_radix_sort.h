#ifndef __UNIFORM_RADIX_SORT_H__
#define __UNIFORM_RADIX_SORT_H__

#include "limits.h"
#include <cstring>
#include <functional>
#include <numeric>
#include <utility>
#include <algorithm>
#include <memory>

namespace dg::uniform_radix_sort::constants{

    static inline constexpr size_t ALPHABET_BIT_SIZE   = CHAR_BIT;
    static inline constexpr size_t RADIX_SIZE          = size_t{1} << ALPHABET_BIT_SIZE;
    static inline constexpr size_t DECAY_RATE          = ALPHABET_BIT_SIZE - 1;
    static inline constexpr size_t PADDING             = 4u; //to-be-adjusted 
    static inline constexpr size_t BLOCK_QS_THRESHOLD  = 0b1111;
}

namespace dg::uniform_radix_sort::types{

    using count_type    = uint32_t;    
}

namespace dg::uniform_radix_sort::utility{

    static constexpr auto temp_buffer_size(size_t n, size_t decay_rate, size_t padding) -> size_t{

        if (n == 0){
            return 0;
        }

        return n + padding + temp_buffer_size(n >> decay_rate, decay_rate, padding);
    }

    template <class T, size_t IDX, std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
    static constexpr auto extract_radix(T val, const std::integral_constant<size_t, IDX>) -> T{
    
        static_assert(constants::RADIX_SIZE - 1 <= std::numeric_limits<T>::max());

        constexpr size_t OFFS   = (sizeof(T) - IDX - 1) * constants::ALPHABET_BIT_SIZE;
        constexpr T BITMASK     = constants::RADIX_SIZE - 1;

        return (val >> OFFS) & BITMASK;
    }
}

namespace dg::uniform_radix_sort::template_sort{

    template <class CallBack, class First, class Second, class ...Args>
    static void insertion_sort(const CallBack& callback, First first, Second second, Args ...args){
        
        if constexpr(sizeof...(Args) == 0){
            callback(std::min(first, second), std::max(first, second));
        } else{
            auto cb_lambda = [=]<class ...AArgs>(AArgs ...aargs){
                callback(std::min(first, second), aargs...);
            };

            insertion_sort(cb_lambda, std::max(first, second), args...);
        }
    } 

    template <class CallBack, class First, class ...Args>
    static void template_sort(const CallBack& callback, First first, Args ...args){

        if constexpr(sizeof...(Args) == 0){
            callback(first);
        } else{
            auto cb_lambda  = [=]<class ...AArgs>(AArgs ...aargs){
                insertion_sort(callback, first, aargs...);
            };

            template_sort(cb_lambda, args...);
        }
    }

    template <class _Ty, size_t SZ_Arg>
    static void template_sort_arr(_Ty * first, const std::integral_constant<size_t, SZ_Arg>&){

        auto sort_cb    = [=]<class ...Args>(Args ...args){
            
            auto fwd_tup        = std::make_tuple(args...);
            const auto idx_seq  = std::make_index_sequence<sizeof...(Args)>{};

            [=]<class Tup, size_t ...IDX>(Tup&& tup, const std::index_sequence<IDX...>&){
                ((first[IDX]  = std::get<IDX>(tup)), ...);
            }(fwd_tup, idx_seq);

        };

        const auto idx_seq    = std::make_index_sequence<SZ_Arg>{};

        [=]<size_t ...IDX>(const std::index_sequence<IDX...>&){
            template_sort(sort_cb, first[IDX]...);
        }(idx_seq);
    }
}

namespace dg::uniform_radix_sort::block_quicksort{

    template <class T, size_t MAX_SZ, size_t PADDING>
    static void block_quicksort(T * first, T * last, 
                                const std::integral_constant<size_t, MAX_SZ>, 
                                const std::integral_constant<size_t, PADDING>){

        static_assert(PADDING >= 1);
        
        if (std::distance(first, last) <= PADDING){
            template_sort::template_sort_arr(first, std::integral_constant<size_t, PADDING>{});
            return;
        }

        T * less_array[MAX_SZ]; 
        size_t less_sz  = 0u;
        T pivot = last[-1]; 

        for (T * iter = first; iter != last; ++iter){
            less_array[less_sz] = iter; 
            less_sz += *iter < pivot;
        }

        for (size_t i = 0; i < less_sz; ++i){
            std::swap(first[i], *less_array[i]);
        }

        std::swap(first[less_sz], last[-1]);
        block_quicksort(first, first + less_sz, std::integral_constant<size_t, MAX_SZ>{}, std::integral_constant<size_t, PADDING>{});
        block_quicksort(first + less_sz + 1, last, std::integral_constant<size_t, MAX_SZ>{}, std::integral_constant<size_t, PADDING>{});
    }
}

namespace dg::uniform_radix_sort{

    using namespace uniform_radix_sort::types;

    template <class T, size_t FIRST, size_t LAST>
    static void radix_sort(T * first, T * last, 
                           count_type * counter,
                           T * tmp, size_t tmp_sz, size_t decay_rate,
                           const std::integral_constant<size_t, FIRST> f_idx,
                           const std::integral_constant<size_t, LAST> l_idx){

        if constexpr(FIRST != LAST){

            size_t sz   = std::distance(first, last);
            
            if (sz <= constants::BLOCK_QS_THRESHOLD){
                block_quicksort::block_quicksort(first, last, std::integral_constant<size_t, constants::BLOCK_QS_THRESHOLD>{}, std::integral_constant<size_t, constants::PADDING>{});
                return;
            } 
            
            if (sz < constants::RADIX_SIZE || sz > tmp_sz){
                std::sort(first, last);
                return;
            }

            auto nxt_counter    = counter + constants::RADIX_SIZE;
            auto nxt_tmp        = tmp + (tmp_sz + constants::PADDING);  
            auto last_tmp       = tmp + sz; 
            auto nxt_tmp_sz     = tmp_sz >> decay_rate;

            std::memset(counter, 0u, constants::RADIX_SIZE * sizeof(count_type));
            std::fill(last_tmp, last_tmp + constants::PADDING, std::numeric_limits<T>::max());
            
            for (auto i = first; i != last; ++i){
                auto radix = utility::extract_radix(*i, f_idx);
                counter[radix] += 1;
            }
            
            for (size_t i = 1; i < constants::RADIX_SIZE; ++i){
                counter[i] += counter[i - 1];
            }
        
            for (auto i = first; i != last; ++i){
                auto radix = utility::extract_radix(*i, f_idx);
                tmp[--counter[radix]] = *i;
            }

            for (size_t i = constants::RADIX_SIZE; i > 0; --i){                
                auto ffirst  = tmp + counter[i - 1];
                auto llast   = (i == constants::RADIX_SIZE) ? last_tmp: tmp + counter[i]; 
                radix_sort(ffirst, llast, nxt_counter, nxt_tmp, nxt_tmp_sz, decay_rate, std::integral_constant<size_t, FIRST + 1>{}, l_idx);
            }

            std::copy(tmp, last_tmp, first);
        }
    } 
    
    template <class T, size_t LAST = sizeof(T)>
    void radix_sort(T * first, T * last, const std::integral_constant<size_t, LAST> = std::integral_constant<size_t, LAST>{}){
    
        auto sz = std::distance(first, last);
        
        if (sz <= constants::BLOCK_QS_THRESHOLD || sz > std::numeric_limits<count_type>::max()){
            std::sort(first, last);
            return;
        }
        
        auto tmp_sz     = utility::temp_buffer_size(static_cast<size_t>(sz), constants::DECAY_RATE, constants::PADDING);
        auto tmp_buf    = std::make_unique<T[]>(tmp_sz); 
        auto count_buf  = std::make_unique<count_type[]>(constants::RADIX_SIZE * LAST);

        radix_sort(first, last, count_buf.get(), tmp_buf.get(), sz, constants::DECAY_RATE, std::integral_constant<size_t, 0u>{}, std::integral_constant<size_t, LAST>{});
    }   
}

#endif