#include <iostream>

int main(int argc, char** argv)
{
#if defined(__AVX2__)
    std::cout << "have avx2" << std::endl;
#else
    std::cout << "error: no avx2" << std::endl;
#endif

    return 0;
}

