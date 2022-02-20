#include<iostream>
#include"fem.h"

int main(int argc, char* argv[]) {

  FEM fem;
  for(int n=0; n<40; n++) {
    for (int i=0; i<200; i++) {
      fem.update_U();
      fem.advance(2.5e-4);
    }
    fem.save_state(n,"./states/");
  }

  return 0;
}
