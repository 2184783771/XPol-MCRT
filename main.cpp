// Entry point for XPol-MCRT. The transport loop lives in run-simulation.cpp.
// The code of STOKES was written
// by Rene Goosmann and Frédéric Marin. A manual can be downloaded from the web page
// www.stokes-program.info. If you have feedback or questions please
// contact admin@stokes-program.info.
//
// STOKES, version 1.2, March 2014


#include "run-state.h"

int main()
{
  RunState run;
  return run.Run();
}
