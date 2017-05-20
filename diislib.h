/*OOParts DIIS library*/
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
extern"C"
{
#include "mkl.h"
}

using namespace std;

bool Cint_DIIS(double diis_conv);
