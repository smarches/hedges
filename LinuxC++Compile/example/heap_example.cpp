#include <iostream>
#include "../heapscheduler.h"

int main() {
    HeapScheduler<double,int> heep(100);
    std::cout << "heap size and capacity are " << heep.size() << ", " << heep.capacity() << '\n';
    heep.push(0.1,400);
    heep.push(22.1,30);
    heep.push(1900,234);
    for(int i=0;i<10;++i) heep.push(8.3*i,i);
    std::cout << "heap size and capacity are " << heep.size() << ", " << heep.capacity() << '\n';
    heep.print();
    // items coming off in descending order based on the first of the value pairs
    for(int i=0;i<4;++i) {
        auto item = heep.pop();
        std::cout << item.first << ", " << item.second << '\n';
    }
    return 0;
}