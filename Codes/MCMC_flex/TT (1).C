
#include <algorithm>
#include <iostream>
#include <vector>
 
using namespace std;

/*
int main()
{
vector<int> v{1,2,3,4,5}; 
return 1;
}
*/

int main()
{
    std::vector<int> v = { 3, 9, 1, 4, 2, 5, 9 };
 
    auto result = std::minmax_element(v.begin(), v.end());
    std::cout << "min element at: " << (result.first - v.begin()) << '\n';
    std::cout << "max element at: " << (result.second - v.begin()) << '\n';
    std::cout << *result.first << '\n';
    std::cout << *result.second << '\n';
    std::vector<int> vv(100,0);
    cout<<vv[0]<<endl;
}

