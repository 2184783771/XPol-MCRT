// Main program of STOKES containing the declaration of all objects
// and the main loop for a model run. The code of STOKES was written
// by Rene Goosmann and Frédéric Marin. A manual can be downloaded from the web page
// www.stokes-program.info. If you have feedback or questions please
// contact admin@stokes-program.info.
//
// STOKES, version 1.2, March 2014


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>  
#ifndef _WIN32
#include <sys/resource.h>
#endif
#include "run-state.h"
#include "photon-stokes-v1.2.h"
#include "observe-stokes-v1.2.h"
#include "randomc-stokes-v1.2.h"                   

using namespace std;

extern int seed;
extern CRandomMersenne RanGen;

#ifndef _WIN32
static void ConfigureStackLimit()
{
  int result;
  const rlim_t kStackSize = 64L * 1024L * 1024L;   // min stack size = 64 Mb
  struct rlimit rl;

  // increase the stack limit if STOKES requires it
  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
  {
      if (rl.rlim_cur < kStackSize)
      {
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0)
          {
              fprintf(stderr, "setrlimit returned result = %d\n", result);
          }
      }
  }
}
#else
static void ConfigureStackLimit()
{
  // On Windows the stack reserve is configured by the linker /STACK option.
}
#endif

static void PrintProgress(long int done,long int total,time_t start,
                          long int lost,long int stuck)
{
  int width,filled,k,percent,elapsed;

  width = 40;
  if (total<=0)
  {
    return;
  }

  percent = int(100.0*done/total);
  filled = int(double(width)*done/total);
  elapsed = int(difftime(time(NULL),start));

  cout << "\r[";
  for (k=0; k<width; k++)
  {
    if (k<filled)
    {
      cout << '#';
    }
    else
    {
      cout << '-';
    }
  }
  cout << "] " << percent << "% "
       << done << "/" << total
       << " elapsed " << elapsed << "s"
       << " lost " << lost
       << " stuck " << stuck << flush;
}

static double ClockSeconds(clock_t ticks)
{
  return double(ticks)/CLOCKS_PER_SEC;
}

static void PrintProfileLine(const char name[],double seconds,double total,
                             long int calls)
{
  double percent,average;

  percent = 0;
  average = 0;
  if (total>0)
  {
    percent = 100.0*seconds/total;
  }
  if (calls>0)
  {
    average = seconds/calls;
  }
  cout << "  " << name << ": " << seconds << " s"
       << " (" << percent << "%)"
       << ", calls " << calls
       << ", avg " << average << " s\n";
}


static const long int MaxPhotonSteps = 5000;

void RunState::ResetStatistics()
{
    i_lost = 0;
    i_lost_corona = 0;
    i_lost_cloud = 0;
    i_stuck = 0;
    total_photon_steps = 0;
    total_corona_scatters = 0;
    max_photon_steps = 0;
    max_corona_scatters = 0;
    corona_scatter_events = 0;
    corona_entered_photons = 0;
    corona_scattered_photons = 0;
    corona_unscattered_photons = 0;
    detected_after_corona_scatter = 0;
    detected_without_corona_scatter = 0;
    lost_after_corona_scatter = 0;
    stuck_after_corona_scatter = 0;
    xrb_seed_emitted = 0;
    xrb_seed_entered_corona = 0;
    xrb_seed_corona_scattered = 0;
    xrb_seed_escaped_without_corona = 0;
    xrb_seed_detected_after_corona = 0;
    disk_hit_photons = 0;
    disk_hit_after_corona = 0;
    disk_hit_without_corona = 0;
    disk_hit_from_corona_source = 0;
    disk_reflection_events = 0;
    disk_reflected_photons = 0;
    detected_after_disk_reflection = 0;
    max_disk_reflections = 0;
    total_disk_hit_radius = 0;
    min_disk_hit_radius = 0;
    max_disk_hit_radius = 0;
    compton_shift_events = 0;
    compton_upshift_events = 0;
    compton_downshift_events = 0;
    total_compton_amp = 0;
    total_compton_ratio = 0;
    max_compton_amp = 0;
    max_photon_energy = 0;
    reset_calls = 0;
    shift_calls = 0;
    check_calls = 0;
    proceed_calls = 0;
    corona_electron_calls = 0;
    cloud_electron_calls = 0;
    register_calls = 0;
    save_calls = 0;
    t_reset = 0;
    t_shift = 0;
    t_check = 0;
    t_proceed = 0;
    t_corona_electron = 0;
    t_cloud_electron = 0;
    t_register = 0;
    t_save = 0;
}

void RunState::SimulatePhoton()
{
      photon_steps = 0;
      photon_corona_scatters = 0;
      photon_disk_reflections = 0;
      photon_entered_corona = false;

      // create a new photon
      tick = clock();
      p.Reset(M);
      t_reset = t_reset+ClockSeconds(clock()-tick);
      reset_calls++;
      if (M.HasXRBSeedSource())
      {
        xrb_seed_emitted++;
      }


      // search for scattering clouds
      while(true)
      {
        tick = clock();
        found_region = p.ShiftToNextScatt(M);
        t_shift = t_shift+ClockSeconds(clock()-tick);
        shift_calls++;
        if (!found_region)
        {
          if (p.ReflectEscapedXRBDisk(M))
          {
            if (p.GetXRBDiskReflectionNum()>photon_disk_reflections)
            {
              photon_disk_reflections = p.GetXRBDiskReflectionNum();
            }
            if (p.Absorbed())
            {
              break;
            }
            continue;
          }
          else
          {
            break;
          }
        }

        photon_steps++;
        if (photon_steps%1000==0)
        {
          PrintProgress(i,M.PhotonNum(),start,i_lost,i_stuck);
          cout << " photon " << i+1 << " steps " << photon_steps << flush;
        }
        if (photon_steps>MaxPhotonSteps)
        {
          p.SetAbsorbed();
          i_stuck++;
          break;
        }

        // if for numerical reasons the photon is not inside the cloud
        // then destroy it
        tick = clock();
	if (p.InCorona() && !p.InsideCorona(M.GetCorona(),M.GetGap()))
	{  
          t_check = t_check+ClockSeconds(clock()-tick);
          check_calls++;
          p.SetAbsorbed();
          i_lost++;
          i_lost_corona++;
          break;
     	}
	if (!p.InCorona() && !p.InsideCloud(M.GetCloud(p.GetFlag()),M.GetGap()))
	{  
          t_check = t_check+ClockSeconds(clock()-tick);
          check_calls++;
          p.SetAbsorbed();
          i_lost++;
          i_lost_cloud++;
          break;
     	}
        t_check = t_check+ClockSeconds(clock()-tick);
        check_calls++;
        if (p.InCorona())
        {
          if (!photon_entered_corona)
          {
            corona_entered_photons++;
            photon_entered_corona = true;
          }
        }

	// photon meets scattering region, proceeds
	// and eventually becomes scattered
        while(true)
        {
          tick = clock();
          proceeded = p.Proceed(M);
          t_proceed = t_proceed+ClockSeconds(clock()-tick);
          proceed_calls++;
          if (!proceeded)
          {
            break;
          }

          photon_steps++;
          if (photon_steps%1000==0)
          {
            PrintProgress(i,M.PhotonNum(),start,i_lost,i_stuck);
            cout << " photon " << i+1 << " steps " << photon_steps << flush;
          }
          if (photon_steps>MaxPhotonSteps)
          {
            p.SetAbsorbed();
            i_stuck++;
            break;
          }
          if (p.GetXRBDiskReflectionNum()>photon_disk_reflections)
          {
            photon_disk_reflections = p.GetXRBDiskReflectionNum();
          }
          if (p.InCorona())
          {
            lambda_before = p.GetLambda();
            tick = clock();
            p.ElectronScattering(M);
            t_corona_electron = t_corona_electron+ClockSeconds(clock()-tick);
            lambda_after = p.GetLambda();
            if (M.GetCorona().compton && lambda_before>0 && lambda_after>0)
            {
              amp = lambda_before/lambda_after;
              if (amp!=1)
              {
                compton_shift_events++;
                total_compton_ratio = total_compton_ratio+amp;
              }
              if (amp>1)
              {
                compton_upshift_events++;
                total_compton_amp = total_compton_amp+amp;
                if (amp>max_compton_amp)
                {
                  max_compton_amp = amp;
                }
              }
              if (amp<1)
              {
                compton_downshift_events++;
              }
              if (12.3984193/lambda_after>max_photon_energy)
              {
                max_photon_energy = 12.3984193/lambda_after;
              }
            }
            corona_electron_calls++;
            corona_scatter_events++;
            photon_corona_scatters++;
          }
          else
          {
            tick = clock();
            p.ElectronScattering(M);
            t_cloud_electron = t_cloud_electron+ClockSeconds(clock()-tick);
            cloud_electron_calls++;
          }
        }
        if (p.Absorbed()) break;

      }

  RecordPhoton();
}

void RunState::RecordPhoton()
{
      // photon escaped or absorbed
      if (p.GetXRBDiskReflectionNum()>photon_disk_reflections)
      {
        photon_disk_reflections = p.GetXRBDiskReflectionNum();
      }
      if (!p.Absorbed())
      {
           if (p.HasHitXRBDisk(disk_hit_radius) ||
               p.WouldHitXRBDisk(M,disk_hit_radius))
           {
             disk_hit_photons++;
             total_disk_hit_radius = total_disk_hit_radius+disk_hit_radius;
             if (disk_hit_photons==1 || disk_hit_radius<min_disk_hit_radius)
             {
               min_disk_hit_radius = disk_hit_radius;
             }
             if (disk_hit_radius>max_disk_hit_radius)
             {
               max_disk_hit_radius = disk_hit_radius;
             }
             if (photon_corona_scatters>0)
             {
               disk_hit_after_corona++;
             }
             else
             {
               disk_hit_without_corona++;
             }
             if (M.HasCoronaSource())
             {
               disk_hit_from_corona_source++;
             }
           }
           tick = clock();
           MT.RegisterPhoton(p,M);
           t_register = t_register+ClockSeconds(clock()-tick);
          register_calls++;
          if (photon_corona_scatters>0)
          {
            detected_after_corona_scatter++;
            if (M.HasXRBSeedSource())
            {
              xrb_seed_detected_after_corona++;
            }
          }
          else
          {
            detected_without_corona_scatter++;
            if (M.HasXRBSeedSource() && !photon_entered_corona)
            {
              xrb_seed_escaped_without_corona++;
            }
          }
          if (photon_disk_reflections>0)
          {
            detected_after_disk_reflection++;
          }
      }
      else if (photon_corona_scatters>0)
      {
        if (photon_steps>MaxPhotonSteps)
        {
          stuck_after_corona_scatter++;
        }
        else
        {
          lost_after_corona_scatter++;
        }
      }
      total_photon_steps = total_photon_steps+photon_steps;
      total_corona_scatters = total_corona_scatters+photon_corona_scatters;
      if (photon_steps>max_photon_steps)
      {
        max_photon_steps = photon_steps;
      }
      if (photon_corona_scatters>max_corona_scatters)
      {
        max_corona_scatters = photon_corona_scatters;
      }
      if (p.GetXRBDiskReflectionNum()>photon_disk_reflections)
      {
        photon_disk_reflections = p.GetXRBDiskReflectionNum();
      }
      disk_reflection_events = disk_reflection_events+photon_disk_reflections;
      if (photon_disk_reflections>0)
      {
        disk_reflected_photons++;
      }
      if (photon_disk_reflections>max_disk_reflections)
      {
        max_disk_reflections = photon_disk_reflections;
      }
      if (photon_corona_scatters>0)
      {
        corona_scattered_photons++;
        if (M.HasXRBSeedSource())
        {
          xrb_seed_corona_scattered++;
        }
      }
      else
      {
        corona_unscattered_photons++;
      }
}

int RunState::Run()
{
  ConfigureStackLimit();

  // initialize model parameters
  if (M.Initialize())
  {
    if (M.HasFixedRandomSeed())
    {
      seed = M.GetRandomSeed();
      RanGen.RandomInit(seed);
    }
    srand(seed);
    cout << "Random seed set to " << seed;
    if (!M.HasFixedRandomSeed())
    {
      cout << " (system time)";
    }
    cout << '\n';

    // STOKES running time
    start = time(NULL);

    cout << "XPol-MCRT: starting photon transport...\n";

    ResetStatistics();
    progress_step = M.PhotonNum()/100;
    if (progress_step<1)
    {
      progress_step = 1;
    }
    PrintProgress(0,M.PhotonNum(),start,i_lost,i_stuck);

    // main loop to sample photons
    for (i=0; i<M.PhotonNum(); i++)
    {
      SimulatePhoton();
      // regularly show photon number and save results
      if ((i+1)%10000000==0)
      {
        cout << "\n... " << (i+1) << " photons sampled\n";
        tick = clock();
        MT.SaveResults(i+1,M);
        t_save = t_save+ClockSeconds(clock()-tick);
        save_calls++;
        PrintProgress(i+1,M.PhotonNum(),start,i_lost,i_stuck);
      }

      if ((i+1)%progress_step==0 || i+1==M.PhotonNum())
      {
        PrintProgress(i+1,M.PhotonNum(),start,i_lost,i_stuck);
      }

    }
    cout << "\n";

    // STOKES running time
    end = time(NULL);
    int time = int(difftime(end, start)); 
    int days = time/86400;
    int hours = (time-days*86400)/3600;
    int minutes = (time-(days*86400+hours*3600))/60;
    int secondes = time-(days*86400+hours*3600+minutes*60);
    cout << "XPol-MCRT simulation completed in " << days << " days " << hours << " hours " << minutes << " minutes " << secondes << " seconds\n";

    // save final results
    tick = clock();
    MT.SaveResults(i,M);
    t_save = t_save+ClockSeconds(clock()-tick);
    save_calls++;
    PrintSummary();
    SaveMetadata(time);
    cout << "XPol-MCRT finished.\n\n\n";
  }
  return 0;
}
