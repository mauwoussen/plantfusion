from lgrass import param_reproduction_functions as prf

from plantfusion.lgrass_wrapper import Lgrass_wrapper, lgrass_soil_domain

from plantfusion.l_egume_wrapper import L_egume_wrapper
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.soil_wrapper import Soil_wrapper
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer

import os
import time
import datetime


def simulation(in_folder_legume, in_folder_lgrass, out_folder, scenario_legume, scenario_lgrass, write_geo=False, lgrass_graph=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    # planter and indexer
    legume_name = "legume"
    lgrass_name = "lgrass"
    index_log = Indexer(global_order=[legume_name, lgrass_name], legume_names=[legume_name], lgrass_names=[lgrass_name])
    planter = Planter(indexer=index_log, legume_cote={legume_name : 40.}, legume_number_of_plants={legume_name : 64})

    # l-egume
    legume = L_egume_wrapper(
        name=legume_name, indexer=index_log, in_folder=in_folder_legume, out_folder=out_folder, IDusm=scenario_legume, caribu_scene=True, planter=planter
    )

    # l-grass
    lgrass = Lgrass_wrapper(
        name=lgrass_name,
        indexer=index_log,
        in_folder=in_folder_lgrass,
        out_folder=out_folder,
        id_scenario=scenario_lgrass,
        activate_genetic_model=False,
        genetic_model_folder=os.path.join(in_folder_lgrass, "modelgenet"),
        outputs_graphs=lgrass_graph,
    )

    # lighting
    lighting = Light_wrapper(
        lightmodel="caribu",
        indexer=index_log, 
        planter=planter, 
        legume_wrapper=legume,
        sky="turtle46",
        out_folder=out_folder,
        writegeo=write_geo,
    )

    # soil
    soil = Soil_wrapper(out_folder=out_folder, legume_wrapper=legume,  legume_pattern=True, planter=planter)

    try:
        current_time_of_the_system = time.time()
        for t in range(legume.lsystem.derivationLength):
            legume.derive(t)

            thermal_day = lgrass.lsystem.current_day
            lgrass.derive(t)

            scene_legume = legume.light_inputs(elements="triangles")
            scene_lgrass = lgrass.light_inputs()
            lighting.run(
                scenes=[scene_legume, scene_lgrass], day=legume.doy(), parunit="RG"
            )
        
            if thermal_day < lgrass.lsystem.current_day :
                if lgrass.setup["option_morphogenetic_regulation_by_carbone"]:
                    lgrass.light_results(lighting=lighting)

            legume.light_results(legume.energy(), lighting)    
            (
                N_content_roots_per_plant,
                roots_length_per_plant_per_soil_layer,
                plants_soil_parameters,
                plants_light_interception,
            ) = legume.soil_inputs()
            soil.run(
                legume.doy(),
                [N_content_roots_per_plant],
                [roots_length_per_plant_per_soil_layer],
                [plants_soil_parameters],
                [plants_light_interception],
            )
            legume.soil_results(soil.results, planter)

            legume.run()

        execution_time = int(time.time() - current_time_of_the_system)
        print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    finally:
        legume.end()
        lgrass.end()


if __name__ == "__main__":
    in_folder_legume = "inputs_soil_legume"
    in_folder_lgrass = "inputs_lgrass"
    out_folder = "outputs/lgrass_legume"
    scenario_legume = 1122
    scenario_lgrass = 5
    write_geo = True
    lgrass_graph = False

    simulation(in_folder_legume, in_folder_lgrass, out_folder, scenario_legume, scenario_lgrass, write_geo=write_geo, lgrass_graph=lgrass_graph)
