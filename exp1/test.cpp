#include <vector>
#include <iostream>

int main(int argc, char * argb[]){
    std::vector<uint64_t> dataset(10);
    
    for (uint64_t idx; idx < 10; idx++){
        dataset.emplace_back(idx);
    }
    
    std::cout << dataset.data() << std::endl;
}