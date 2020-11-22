#include <utility>
#include <type_traits>
#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>

template<class T, class U>
class HeapScheduler {
    using TU = std::pair<T,U>;
public:
	constexpr size_t default_size{1'100'000}; // initial heap size (why?)
    constexpr T bigval{std::numeric_limits<T>::max()};
	std::vector<TU> elems; // each item's timestamp is the sorting criterion so it can be processed LIFO?

	HeapScheduler() {
        // TODO: concept
        static_assert<std::is_arithmetic_v<T>>;
        TU.reserve(default_size);
    }
	void push(T time, U cargo = U{}) { 
	    // pushes a time and cargo onto the heap
        elems.push_back(std::make_pair(time,cargo));
        std::push_heap(std::begin(elems),std::end(elems),[](TU& A,TU& B){
            return A.first < B.first;
        });
	}
	TU pop() { // return top of heap, move last to top, shorten list, sift down
	// pops the next (in order) time and its cargo from the heap
	// returns bigval time, and U() cargo, when heap is empty
        if (!elems.empty()) {
            std::pop_heap(std::begin(elems),std::end(elems));
            return elems.pop_back();
        } else {
            // this is an 'error code'
            return std::make_pair(bigval,U{});
        }
	}
	void rewind() { // zero out the heap w/o changing its size in memory
		elems.assign(default_size, std::pair(bigval,U{}));
        elems.shrink_to_fit();
	}
	void reinit() { // zero out the heap and give back memory
		elems.clear();
	}

	void print() const { 
		std::cout << "heap is:\n";
		for (std::size_t i=0;i < elems.size();++i) {
            std::cout << i << ", " << elems[i].first << " " << elems[i].second << '\n';
        }
		std::cout << "end heap\n";
	}
};
