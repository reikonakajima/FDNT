#include <iostream>
#include <iomanip>
using namespace std;
#include "SCAMPMap.h"

main() {
  int a = 1<<1;
  int b = 1<<8;
  int flags;
  cout << flags << endl;
  flags &= ~(a && b);
  cout << hex << a << " " << b << endl;
  cout << hex << (a && b) << endl;
  cout << hex << ~(a && b) << endl;
  cout << hex << flags << endl;
  cout << "hello world\n";
  exit(0);
}
