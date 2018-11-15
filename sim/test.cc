#include <iostream>
#include <fstream>
#include "G4MaterialPropertyVector.hh"

int main(int argc, char *argv[])
{
    if (argc<2) return 1;
    
    G4MaterialPropertyVector data(argv[1]);
    std::cout << data << std::endl;
    
    return 0;
}
