#include <utility>
#include <type_traits>
#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>

template<class T, class U>
class HeapScheduler {
    using TU = std::pair<T,U>;
    private:
    static constexpr size_t default_size{100'000}; // initial heap size (why?)
    static constexpr T bigval{std::numeric_limits<T>::max()};
    std::vector<TU> elems; // each item's timestamp is the sorting criterion so it can be processed LIFO?

    // TODO: concept. Without the silly 'bigval' behavior we can have any type implementing operator<
    static_assert(std::is_arithmetic_v<T>);
    public:
    HeapScheduler() {
        elems.reserve(default_size);
    }
    HeapScheduler(size_t n) {
        elems.reserve(n);
    }
    size_t size() const {
        return elems.size();
    }
    size_t capacity() const {
        return elems.capacity();
    }
    void push(T time, U cargo = U{}) { 
        elems.emplace_back(std::make_pair(time,cargo));
        std::push_heap(std::begin(elems),std::end(elems),[](TU& A,TU& B){
            return A.first < B.first;
        });
    }
    TU pop() { // return top of heap, move last to top, shorten list, sift down
    // pops the next (in order) time and its cargo from the heap
    // returns bigval time, and U() cargo, when heap is empty
        if (!elems.empty()) {
            std::pop_heap(std::begin(elems),std::end(elems));
            auto ans = elems.back();
            elems.pop_back();
            return ans;
        } else {
            // this is an 'error code'. May want to just 'forward' the std::vector 
            // undefined behavior when pop_back is called on empty container?
            return std::make_pair(bigval,U{});
        }
    }
    void rewind(size_t n = default_size) { // """zero out""" the heap
        elems.assign(n, std::pair(bigval,U{}));
        elems.shrink_to_fit();
    }

    void print() const { 
        std::cout << "heap is:\n";
        for (std::size_t i=0;i < std::min(100ul,elems.size()); ++i) {
            std::cout << i << ", " << elems[i].first << " " << elems[i].second << '\n';
        }
        std::cout << "end heap\n";
    }
    T default_value() constexpr noexcept {
        return bigval;
    }
};
