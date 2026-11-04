#include "run-state.h"
#include <fstream>
#include <iostream>
#include <cstdio>

using namespace std;

extern int seed;

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


template <class T>
static void WriteMetadataValue(ofstream& stream,const char name[],const T& value)
{
  stream << name << "," << value << "\n";
}


static const char* YesNo(bool value)
{
  return value ? "yes" : "no";
}


static const char* CoronaTypeName(int type)
{
  if (type==0) return "slab";
  if (type==1) return "toroidal";
  if (type==2) return "lamp-post";
  return "unknown";
}


static const char* SeedSpectrumName(int spectrum)
{
  if (spectrum==0) return "powerlaw";
  if (spectrum==1) return "blackbody";
  return "unknown";
}


static const char* AngularName(int angular)
{
  if (angular==0) return "isotropic";
  if (angular==1) return "lambert";
  return "unknown";
}


static const char* SeedPolarizationName(int polarization)
{
  if (polarization==0) return "none";
  if (polarization==1) return "chandrasekhar";
  return "unknown";
}


void RunState::PrintSummary()
{
    cout << "XPol-MCRT results saved\n";
    cout << "Number of lost photons: " << i_lost << '\n';
    cout << "  lost in corona: " << i_lost_corona << '\n';
    cout << "  lost in cloud: " << i_lost_cloud << '\n';
    cout << "Photons stopped by step limit: " << i_stuck << '\n';
    if (M.PhotonNum()>0)
    {
      cout << "Average photon steps: " << total_photon_steps/M.PhotonNum() << '\n';
    }
    cout << "Maximum photon steps: " << max_photon_steps << '\n';
    cout << "Corona scattering diagnostics:\n";
    cout << "  photons entering corona: " << corona_entered_photons << '\n';
    cout << "  corona scattering events: " << corona_scatter_events << '\n';
    cout << "  photons with corona scattering: " << corona_scattered_photons << '\n';
    cout << "  photons without corona scattering: " << corona_unscattered_photons << '\n';
    if (M.PhotonNum()>0)
    {
      cout << "  average corona scatterings per photon: "
           << total_corona_scatters/M.PhotonNum() << '\n';
    }
    cout << "  maximum corona scatterings in one photon: "
         << max_corona_scatters << '\n';
    cout << "  detected after corona scattering: "
         << detected_after_corona_scatter << '\n';
    cout << "  detected without corona scattering: "
         << detected_without_corona_scatter << '\n';
    cout << "  lost after corona scattering: "
         << lost_after_corona_scatter << '\n';
    cout << "  stopped after corona scattering: "
         << stuck_after_corona_scatter << '\n';
    if (M.HasXRBSeedSource())
    {
      xrb_seed_entered_corona = corona_entered_photons;
      cout << "XRB disk seed diagnostics:\n";
      cout << "  disk seed photons emitted: "
           << xrb_seed_emitted << '\n';
      cout << "  disk seed photons entering corona: "
           << xrb_seed_entered_corona << '\n';
      cout << "  disk seed photons corona-scattered: "
           << xrb_seed_corona_scattered << '\n';
      cout << "  disk seed photons escaped without corona: "
           << xrb_seed_escaped_without_corona << '\n';
      cout << "  disk seed photons detected after corona: "
           << xrb_seed_detected_after_corona << '\n';
    }
    cout << "Disk hit diagnostics:\n";
    cout << "  escaped photons that would hit disk before detector: "
         << disk_hit_photons << '\n';
    cout << "  disk-hit fraction of all photons: ";
    if (M.PhotonNum()>0)
    {
      cout << double(disk_hit_photons)/M.PhotonNum();
    }
    else
    {
      cout << 0;
    }
    cout << '\n';
    cout << "  disk hits after corona scattering: "
         << disk_hit_after_corona << '\n';
    cout << "  disk hits without corona scattering: "
         << disk_hit_without_corona << '\n';
    cout << "  disk hits from corona-source photons: "
         << disk_hit_from_corona_source << '\n';
    cout << "  disk reflection events: "
         << disk_reflection_events << '\n';
    cout << "  photons with disk reflection: "
         << disk_reflected_photons << '\n';
    cout << "  detected after disk reflection: "
         << detected_after_disk_reflection << '\n';
    cout << "  maximum disk reflections in one photon: "
         << max_disk_reflections << '\n';
    if (disk_hit_photons>0)
    {
      cout << "  average disk hit radius: "
           << total_disk_hit_radius/disk_hit_photons << '\n';
      cout << "  minimum disk hit radius: "
           << min_disk_hit_radius << '\n';
      cout << "  maximum disk hit radius: "
           << max_disk_hit_radius << '\n';
    }
    if (M.GetXRBSeedDisk().reflection)
    {
      cout << "  disk reflection handling: single Thomson/";
      if (M.GetXRBSeedDisk().reflectAngular==1)
      {
        cout << "Lambert reflection\n";
      }
      else
      {
        cout << "specular reflection\n";
      }
    }
    else
    {
      cout << "  disk reflection handling: disabled\n";
    }
    if (M.GetCorona().compton)
    {
      cout << "Compton diagnostics:\n";
      cout << "  corona Compton energy changes: "
           << compton_shift_events << '\n';
      cout << "  up-scattering events: "
           << compton_upshift_events << '\n';
      cout << "  down-scattering events: "
           << compton_downshift_events << '\n';
      if (compton_shift_events>0)
      {
        cout << "  average energy ratio after/before: "
             << total_compton_ratio/compton_shift_events << '\n';
      }
      if (compton_upshift_events>0)
      {
        cout << "  average up-scattering amplification: "
             << total_compton_amp/compton_upshift_events << '\n';
      }
      cout << "  max energy amplification: "
           << max_compton_amp << '\n';
      cout << "  max photon energy keV: "
           << max_photon_energy << '\n';
    }
    cout << "Runtime profile:\n";
    double prof_total = t_reset+t_shift+t_check+t_proceed+
                        t_corona_electron+t_cloud_electron+
                        t_register+t_save;
    PrintProfileLine("Reset photons",t_reset,prof_total,reset_calls);
    PrintProfileLine("Find next region",t_shift,prof_total,shift_calls);
    PrintProfileLine("Inside checks",t_check,prof_total,check_calls);
    PrintProfileLine("Proceed in region",t_proceed,prof_total,proceed_calls);
    PrintProfileLine("Corona electron scattering",t_corona_electron,prof_total,corona_electron_calls);
    PrintProfileLine("Cloud electron scattering",t_cloud_electron,prof_total,cloud_electron_calls);
    PrintProfileLine("Register detector",t_register,prof_total,register_calls);
    PrintProfileLine("Save results",t_save,prof_total,save_calls);
}

void RunState::SaveMetadata(int elapsedSeconds)
{
    {
      char outfile[300],metafile[360];
      CoronaRegion Corona = M.GetCorona();
      XRBSeedDisk SeedDisk = M.GetXRBSeedDisk();

      M.GetOutputFile(outfile);
      sprintf(metafile,"%s_metadata.csv",outfile);
      ofstream meta(metafile,ios::out);
      if (meta.fail()==0)
      {
        meta << "key,value\n";
        WriteMetadataValue(meta,"output_file",outfile);
        WriteMetadataValue(meta,"random_seed",seed);
        WriteMetadataValue(meta,"random_seed_fixed",YesNo(M.HasFixedRandomSeed()));
        WriteMetadataValue(meta,"spin_case_enabled",YesNo(M.HasSpinCase()));
        WriteMetadataValue(meta,"spin_case",M.GetSpinCase());
        WriteMetadataValue(meta,"spin_a_star",M.HasSpinCase() ? (M.GetSpinCase()==1 ? 0.998 : 0.0) : 0.0);
        WriteMetadataValue(meta,"spin_isco_rg",M.GetSpinISCO());
        WriteMetadataValue(meta,"sky_axis_pa_deg_east_of_north",M.GetSkyAxisPA());
        WriteMetadataValue(meta,"inclination_enabled",YesNo(M.HasInclination()));
        WriteMetadataValue(meta,"inclination_input_deg",M.GetInclinationDeg());
        WriteMetadataValue(meta,"single_inclination_band_summary",YesNo(M.HasInclination()));
        if (M.HasInclination())
        {
          WriteMetadataValue(meta,"summary_energy_bands_kev","2-3;3-5;5-8;2-8 total");
        }
        WriteMetadataValue(meta,"source_ra_deg",M.GetSourceRA());
        WriteMetadataValue(meta,"source_dec_deg",M.GetSourceDec());
        WriteMetadataValue(meta,"source_distance_kpc",M.GetSourceDistanceKpc());
        WriteMetadataValue(meta,"black_hole_mass_solar",M.GetBlackHoleMassSolar());
        WriteMetadataValue(meta,"photon_num",M.PhotonNum());
        WriteMetadataValue(meta,"lambda_min_angstrom",M.LambdaMin());
        WriteMetadataValue(meta,"lambda_max_angstrom",M.LambdaMax());
        WriteMetadataValue(meta,"spectral_bins",M.res());
        WriteMetadataValue(meta,"theta_view_bins",M.ScreenNumTheta());
        WriteMetadataValue(meta,"phi_view_bins",M.ScreenNumPhi());
        WriteMetadataValue(meta,"plane_symmetry",YesNo(M.GetPlaneSymmetry()));
        WriteMetadataValue(meta,"axis_symmetry",YesNo(M.GetAxisSymmetry()));
        WriteMetadataValue(meta,"res_x",M.GetResX());
        WriteMetadataValue(meta,"res_y",M.GetResY());
        WriteMetadataValue(meta,"save_v",YesNo(M.ShouldSaveV()));
        WriteMetadataValue(meta,"save_components",YesNo(M.ShouldSaveComponents()));
        WriteMetadataValue(meta,"save_stokes_files",YesNo(M.ShouldSaveStokesFiles()));
        WriteMetadataValue(meta,"save_map",YesNo(M.ShouldSaveMap()));
        WriteMetadataValue(meta,"map_theta_index",M.GetMapThetaIndex());
        WriteMetadataValue(meta,"map_energy_min_kev",M.GetMapEnergyMin());
        WriteMetadataValue(meta,"map_energy_max_kev",M.GetMapEnergyMax());
        WriteMetadataValue(meta,"has_corona",YesNo(M.HasCorona()));
        WriteMetadataValue(meta,"corona_type",CoronaTypeName(Corona.type));
        WriteMetadataValue(meta,"corona_source",YesNo(Corona.source));
        WriteMetadataValue(meta,"corona_compton",YesNo(Corona.compton));
        WriteMetadataValue(meta,"corona_pos_x",Corona.pos[0]);
        WriteMetadataValue(meta,"corona_pos_y",Corona.pos[1]);
        WriteMetadataValue(meta,"corona_pos_z",Corona.pos[2]);
        WriteMetadataValue(meta,"corona_radius",Corona.radius);
        WriteMetadataValue(meta,"corona_height",Corona.height);
        WriteMetadataValue(meta,"corona_tau",Corona.tau);
        WriteMetadataValue(meta,"corona_density",Corona.density);
        WriteMetadataValue(meta,"corona_temp_kev",Corona.kTe);
        WriteMetadataValue(meta,"corona_gamma",Corona.gamma);
        WriteMetadataValue(meta,"corona_energy_min_kev",Corona.emin);
        WriteMetadataValue(meta,"corona_energy_max_kev",Corona.emax);
        WriteMetadataValue(meta,"xrb_seed_source",YesNo(SeedDisk.source));
        WriteMetadataValue(meta,"xrb_seed_spectrum",SeedSpectrumName(SeedDisk.spectrum));
        WriteMetadataValue(meta,"xrb_seed_angular",AngularName(SeedDisk.angular));
        WriteMetadataValue(meta,"xrb_seed_polarization",
                           SeedPolarizationName(SeedDisk.polarization));
        WriteMetadataValue(meta,"xrb_seed_rin",SeedDisk.rin);
        WriteMetadataValue(meta,"xrb_seed_rout",SeedDisk.rout);
        WriteMetadataValue(meta,"xrb_seed_gamma",SeedDisk.gamma);
        WriteMetadataValue(meta,"xrb_seed_tin_kev",SeedDisk.tin);
        WriteMetadataValue(meta,"xrb_seed_energy_min_kev",SeedDisk.emin);
        WriteMetadataValue(meta,"xrb_seed_energy_max_kev",SeedDisk.emax);
        WriteMetadataValue(meta,"disk_reflection",YesNo(SeedDisk.reflection));
        WriteMetadataValue(meta,"disk_reflect_angular",AngularName(SeedDisk.reflectAngular));
        WriteMetadataValue(meta,"disk_reflect_rin",SeedDisk.reflectRin);
        WriteMetadataValue(meta,"disk_reflect_rout",SeedDisk.reflectRout);
        WriteMetadataValue(meta,"disk_reflect_albedo",SeedDisk.albedo);
        WriteMetadataValue(meta,"elapsed_seconds",elapsedSeconds);
        WriteMetadataValue(meta,"lost_photons",i_lost);
        WriteMetadataValue(meta,"lost_in_corona",i_lost_corona);
        WriteMetadataValue(meta,"lost_in_cloud",i_lost_cloud);
        WriteMetadataValue(meta,"stopped_by_step_limit",i_stuck);
        WriteMetadataValue(meta,"average_photon_steps",M.PhotonNum()>0 ? total_photon_steps/M.PhotonNum() : 0);
        WriteMetadataValue(meta,"maximum_photon_steps",max_photon_steps);
        WriteMetadataValue(meta,"photons_entering_corona",corona_entered_photons);
        WriteMetadataValue(meta,"corona_scattering_events",corona_scatter_events);
        WriteMetadataValue(meta,"photons_with_corona_scattering",corona_scattered_photons);
        WriteMetadataValue(meta,"photons_without_corona_scattering",corona_unscattered_photons);
        WriteMetadataValue(meta,"average_corona_scatterings_per_photon",M.PhotonNum()>0 ? total_corona_scatters/M.PhotonNum() : 0);
        WriteMetadataValue(meta,"maximum_corona_scatterings_one_photon",max_corona_scatters);
        WriteMetadataValue(meta,"detected_after_corona_scattering",detected_after_corona_scatter);
        WriteMetadataValue(meta,"detected_without_corona_scattering",detected_without_corona_scatter);
        WriteMetadataValue(meta,"lost_after_corona_scattering",lost_after_corona_scatter);
        WriteMetadataValue(meta,"stopped_after_corona_scattering",stuck_after_corona_scatter);
        WriteMetadataValue(meta,"xrb_seed_photons_emitted",xrb_seed_emitted);
        WriteMetadataValue(meta,"xrb_seed_photons_entering_corona",xrb_seed_entered_corona);
        WriteMetadataValue(meta,"xrb_seed_photons_corona_scattered",xrb_seed_corona_scattered);
        WriteMetadataValue(meta,"xrb_seed_photons_escaped_without_corona",xrb_seed_escaped_without_corona);
        WriteMetadataValue(meta,"xrb_seed_photons_detected_after_corona",xrb_seed_detected_after_corona);
        WriteMetadataValue(meta,"disk_hit_photons",disk_hit_photons);
        WriteMetadataValue(meta,"disk_hit_fraction_all_photons",M.PhotonNum()>0 ? double(disk_hit_photons)/M.PhotonNum() : 0);
        WriteMetadataValue(meta,"disk_hits_after_corona_scattering",disk_hit_after_corona);
        WriteMetadataValue(meta,"disk_hits_without_corona_scattering",disk_hit_without_corona);
        WriteMetadataValue(meta,"disk_hits_from_corona_source",disk_hit_from_corona_source);
        WriteMetadataValue(meta,"disk_reflection_events",disk_reflection_events);
        WriteMetadataValue(meta,"photons_with_disk_reflection",disk_reflected_photons);
        WriteMetadataValue(meta,"detected_after_disk_reflection",detected_after_disk_reflection);
        WriteMetadataValue(meta,"maximum_disk_reflections_one_photon",max_disk_reflections);
        WriteMetadataValue(meta,"average_disk_hit_radius",disk_hit_photons>0 ? total_disk_hit_radius/disk_hit_photons : 0);
        WriteMetadataValue(meta,"minimum_disk_hit_radius",disk_hit_photons>0 ? min_disk_hit_radius : 0);
        WriteMetadataValue(meta,"maximum_disk_hit_radius",disk_hit_photons>0 ? max_disk_hit_radius : 0);
        WriteMetadataValue(meta,"compton_energy_change_events",compton_shift_events);
        WriteMetadataValue(meta,"compton_upscattering_events",compton_upshift_events);
        WriteMetadataValue(meta,"compton_downscattering_events",compton_downshift_events);
        WriteMetadataValue(meta,"average_compton_energy_ratio",compton_shift_events>0 ? total_compton_ratio/compton_shift_events : 0);
        WriteMetadataValue(meta,"average_upscattering_amplification",compton_upshift_events>0 ? total_compton_amp/compton_upshift_events : 0);
        WriteMetadataValue(meta,"maximum_compton_amplification",max_compton_amp);
        WriteMetadataValue(meta,"maximum_photon_energy_kev",max_photon_energy);
        WriteMetadataValue(meta,"profile_reset_seconds",t_reset);
        WriteMetadataValue(meta,"profile_find_next_region_seconds",t_shift);
        WriteMetadataValue(meta,"profile_inside_checks_seconds",t_check);
        WriteMetadataValue(meta,"profile_proceed_region_seconds",t_proceed);
        WriteMetadataValue(meta,"profile_corona_electron_seconds",t_corona_electron);
        WriteMetadataValue(meta,"profile_cloud_electron_seconds",t_cloud_electron);
        WriteMetadataValue(meta,"profile_register_detector_seconds",t_register);
        WriteMetadataValue(meta,"profile_save_results_seconds",t_save);
        meta.close();
        cout << "Run metadata saved to " << metafile << '\n';
      }
      else
      {
        cout << "Warning: unable to save run metadata to " << metafile << '\n';
      }
    }
}
