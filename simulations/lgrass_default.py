import os
import time
import datetime

from lgrass import param_reproduction_functions as prf

from plantfusion.lgrass_wrapper import Lgrass_wrapper, lgrass_soil_domain
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter

def simulation(in_folder, genetic_model_folder, out_folder, id_scenario=0, write_geo=False, graphs=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    plants_name = "lgrass"
    index_log = Indexer(global_order=[plants_name], lgrass_names=[plants_name])
    planter = Planter(generation_type="default", indexer=index_log, xy_plane=lgrass_soil_domain())

    lgrass = Lgrass_wrapper(name=plants_name,
                            indexer=index_log,
                            in_folder=in_folder, 
                            out_folder=out_folder, 
                            id_scenario=id_scenario,
                            activate_genetic_model=True, 
                            genetic_model_folder=genetic_model_folder,
                            outputs_graphs=graphs)


    lighting = Light_wrapper(
        lightmodel="caribu", 
        out_folder=out_folder, 
        sky="turtle46",
        planter=planter, 
        indexer=index_log,
        writegeo=write_geo
    )

    current_time_of_the_system = time.time()
    prf.rungenet(lgrass.genet_src, lgrass.genet_dst, lgrass.genet_exe, None, 0)
    
    # genetic generation loop
    for i in range(1, lgrass.setup['num_gener'] + 1):
        lgrass = Lgrass_wrapper(name=plants_name,
                                indexer=index_log,
                                in_folder=in_folder, 
                                out_folder=out_folder, 
                                id_scenario=id_scenario, 
                                activate_genetic_model=True, 
                                genetic_model_folder=genetic_model_folder,
                                generation_index=i,
                                outputs_graphs=graphs)
        
        # daily loop
        for t in range(0, lgrass.lsystem.derivationLength):
            lgrass.derive(t)

            if lgrass.setup["option_morphogenetic_regulation_by_carbone"] :
                scene_lgrass = lgrass.light_inputs()
                lighting.run(energy=lgrass.energy(), scenes=[scene_lgrass])
                lgrass.light_results(lighting=lighting)
            
            # empty run
            # lgrass.run()

        lgrass.end()

        if lgrass.setup["option_reproduction"] != "False":
            prf.rungenet(lgrass.genet_src, lgrass.genet_dst, lgrass.genet_exe, lgrass.genet_mat, 1)


    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


if __name__ == "__main__":
    in_folder = "inputs_lgrass"
    genetic_model_folder = os.path.join(os.getcwd(), in_folder, "modelgenet")
    out_folder = "outputs/lgrass_default"
    write_geo = False
    graphs = True
    scenario_index = 0

    simulation(in_folder, genetic_model_folder, out_folder, scenario_index, write_geo, graphs)
