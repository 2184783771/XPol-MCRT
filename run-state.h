#ifndef XPOL_RUN_STATE_H
#define XPOL_RUN_STATE_H

#include <ctime>
#include "model-stokes-v1.2.h"
#include "photon-stokes-v1.2.h"
#include "observe-stokes-v1.2.h"

class RunState
{
public:
  int Run();

private:
  void ResetStatistics();
  void SimulatePhoton();
  void RecordPhoton();
  void PrintSummary();
  void SaveMetadata(int elapsedSeconds);

  Model M;
  Monitor MT;
  Photon p;

  // Run and per-photon counters.
  long int i,i_lost,i_lost_corona,i_lost_cloud,i_stuck;
  long int photon_steps,progress_step,max_photon_steps;
  long int reset_calls,shift_calls,check_calls,proceed_calls;
  long int corona_electron_calls,cloud_electron_calls;
  long int register_calls,save_calls;
  long int photon_corona_scatters,max_corona_scatters;
  long int corona_scatter_events,corona_entered_photons;
  long int corona_scattered_photons,corona_unscattered_photons;
  long int detected_after_corona_scatter,detected_without_corona_scatter;
  long int lost_after_corona_scatter,stuck_after_corona_scatter;
  long int xrb_seed_emitted,xrb_seed_entered_corona;
  long int xrb_seed_corona_scattered,xrb_seed_escaped_without_corona;
  long int xrb_seed_detected_after_corona;
  long int disk_hit_photons,disk_hit_after_corona,disk_hit_without_corona;
  long int disk_hit_from_corona_source;
  long int disk_reflection_events,disk_reflected_photons;
  long int detected_after_disk_reflection,max_disk_reflections;
  long int photon_disk_reflections;
  long int compton_shift_events,compton_upshift_events,compton_downshift_events;

  // Aggregate physical diagnostics.
  double total_photon_steps;
  double total_corona_scatters;
  double disk_hit_radius,total_disk_hit_radius;
  double min_disk_hit_radius,max_disk_hit_radius;
  double lambda_before,lambda_after,amp,total_compton_amp,total_compton_ratio;
  double max_compton_amp,max_photon_energy;

  // Timings and current photon state.
  double t_reset,t_shift,t_check,t_proceed,t_corona_electron;
  double t_cloud_electron,t_register,t_save;
  bool found_region,proceeded,photon_entered_corona;
  clock_t tick;
  time_t start,end;
};

#endif
