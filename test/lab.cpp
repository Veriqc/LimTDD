#include <iostream>
#include <string>
#include <ctime>
int test(){
    int res;
    res = 5;
    std::cout << "int res fun: " << res <<" in: " << &res << std::endl;
    return res;
}
int* test_() {
    int* res = new int;
    *res = 5;
    std::cout << "int* res fun: " << *res << " in: " << res << std::endl;
    return res;
}
int& test__() {
    static int res ;
    res = 5;
    std::cout << "int* res fun: " << res << " in: " << &res << std::endl;
    return res;
}
int main(){
    int a{};
    std::cout << "a: " << a <<" in: " << &a << std::endl;
    int b;
    std::cout << "b: " << b <<" in: " << &b << std::endl;
    auto c = test();
    std::cout << "c: " << c <<" in: " << &c << std::endl;
    auto d = test_();
    auto& dd = *test_();
    std::cout << "d: " << *d <<" in: " << d << std::endl;
    std::cout << "dd: " << dd <<" in: " << &dd << std::endl;
    int* e = test_();
    std::cout << "e: " << *e <<" in: " << e << std::endl;
    auto& f = test__();
    std::cout << "f: " << f <<" in: " << &f << std::endl;
    int t;
    std::cin >> t;
    return 0;
}