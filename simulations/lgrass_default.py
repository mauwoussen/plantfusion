import os
import time

from plantfusion.lgrass_wrapper import Lgrass_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter

def simulation(in_folder, out_folder, simulation_length=500, write_geo=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    plants_name = "lgrass"
    index_log = Indexer(global_order=[plants_name], lgrass_names=[plants_name])
    planter = Planter(generation_type="default", indexer=index_log)
    
    lgrass = Lgrass_wrapper(in_folder=in_folder, out_folder=out_folder)

    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky=sky,
        planter=planter, 
        indexer=index_log,
        writegeo=write_geo
    )

    current_time_of_the_system = time.time()
    for t in range(lgrass.start_time, simulation_length, lgrass.SENESClgrass_TIMESTEP):
        if (t % light_timestep == 0) and (lgrass.PARi_next_hours(t) > 0):
            lgrass_input, stems = lgrass.light_inputs(planter)
            lighting.run(scenes=[lgrass_input], day=lgrass.doy(t), hour=lgrass.hour(t), parunit="micromol.m-2.s-1", stems=stems)
            lgrass.light_results(energy=lgrass.energy(t), lighting=lighting)

        lgrass.run(t)

    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


    lgrass.end(run_postprocessing=run_postprocessing)


if __name__ == "__main__":
    in_folder = "inputs_fspmlgrass"
    out_folder = "outputs/cnlgrass_default"
    simulation_length = 500
    write_geo = True

    simulation(in_folder, out_folder, simulation_length, write_geo=write_geo)
