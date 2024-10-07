import os
import time
import datetime

from lgrass import param_reproduction_functions as prf

from plantfusion.lgrass_wrapper import Lgrass_wrapper, lgrass_soil_domain
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter


def simulation(in_folder, genetic_model_folder, out_folder, id_scenario=0, write_geo=False, graphs=False):
    """

    Parameters
    ----------
    in_folder : str
        input folder path
    genetic_model_folder : str
        folder path with genetic model parameters and executables
    out_folder : str
        output folder path
    id_scenario : int
        scenario index in plan_simulation.csv file
    write_geo : bool
        write 3D VTK file of the scene
    graphs : bool
        write graphs from output results


    """
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")

    plants_name = "lgrass"
    index_log = Indexer(global_order=[plants_name], lgrass_names=[plants_name])
    planter_9plants_plantfusion = Planter(generation_type="default", indexer=index_log, xy_plane=lgrass_soil_domain(espacement=1., rows=3, columns=3), save_plant_positions=True)
    planter_9plants_default = Planter(generation_type="default", indexer=index_log, xy_plane=lgrass_soil_domain(espacement=1., rows=3, columns=3), save_plant_positions=True)
    planter_1plant = Planter(generation_type="default", indexer=index_log, xy_plane=lgrass_soil_domain(espacement=1., rows=1, columns=1), save_plant_positions=True)

    lgrass_9plants_plantfusion = Lgrass_wrapper(
        name=plants_name,
        indexer=index_log,
        planter=planter_9plants_plantfusion,
        in_folder=in_folder,
        out_folder=out_folder,
        id_scenario=id_scenario,
        lai="plantfusion",
        number_of_plants=9,
        activate_genetic_model=True,
        genetic_model_folder=genetic_model_folder,
        outputs_graphs=graphs,
    )

    lgrass_9plants_default = Lgrass_wrapper(
        name=plants_name,
        indexer=index_log,
        planter=planter_9plants_default,
        in_folder=in_folder,
        out_folder=out_folder,
        id_scenario=id_scenario,
        lai="lgrass",
        number_of_plants=9,
        activate_genetic_model=True,
        genetic_model_folder=genetic_model_folder,
        outputs_graphs=graphs,
    )

    lgrass_1plant = Lgrass_wrapper(
        name=plants_name,
        indexer=index_log,
        planter=planter_1plant,
        in_folder=in_folder,
        out_folder=out_folder,
        id_scenario=id_scenario,
        lai="plantfusion",
        number_of_plants=1,
        activate_genetic_model=True,
        genetic_model_folder=genetic_model_folder,
        outputs_graphs=graphs,
    )

    lighting_9plants = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        sky="turtle46",
        planter=planter_9plants_plantfusion,
        indexer=index_log,
        writegeo=write_geo,
    )

    lighting_1plant = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        sky="turtle46",
        planter=planter_1plant,
        indexer=index_log,
        writegeo=write_geo,
    )

    scanning_ray = 0.015
    planter_9plants_plantfusion.scan_nearest_plants_neighours(scanning_ray)
    planter_9plants_default.scan_nearest_plants_neighours(scanning_ray)
    planter_1plant.scan_nearest_plants_neighours(scanning_ray)

    current_time_of_the_system = time.time()

    # daily loop
    for t in range(0, lgrass_1plant.lsystem.derivationLength):
        thermal_day = lgrass_1plant.lsystem.current_day
        lgrass_1plant.derive(t)
        lgrass_9plants_default.derive(t)
        lgrass_9plants_plantfusion.derive(t)

        # execute CARIBU on each new thermal day
        if thermal_day < lgrass_9plants_default.lsystem.current_day or t == 0:
            scene_lgrass = lgrass_1plant.light_inputs()
            lighting_1plant.run(energy=lgrass_1plant.energy(), scenes=[scene_lgrass])
            lgrass_1plant.light_results(lighting=lighting_1plant)

            scene_lgrass = lgrass_9plants_default.light_inputs()
            lighting_9plants.run(energy=lgrass_9plants_default.energy(), scenes=[scene_lgrass])
            lgrass_9plants_default.light_results(lighting=lighting_9plants)

            scene_lgrass = lgrass_9plants_plantfusion.light_inputs()
            lighting_9plants.run(energy=lgrass_9plants_plantfusion.energy(), scenes=[scene_lgrass])
            lgrass_9plants_plantfusion.light_results(lighting=lighting_9plants)
        
        # compute local LAI
        lgrass_1plant.run(planter_1plant)
        lgrass_9plants_default.run(planter_9plants_default)
        lgrass_9plants_plantfusion.run(planter_9plants_plantfusion)

        # print
        print("--- Local LAI comparaison ---")
        print("plant ID \t one plant plantfusion \t 9 plants lgrass \t 9 plants plantfusion")
        for k,v in lgrass_1plant.lsystem.rapportS9_SSol_dict.items() :
            for k2, v2 in lgrass_9plants_default.lsystem.rapportS9_SSol_dict.items():
                print("%i \t %f \t %f \t %f"%(k2, v, v2, lgrass_9plants_plantfusion.lsystem.rapportS9_SSol_dict[k2]))
        print("\n")

    lgrass_1plant.end()
    lgrass_9plants_default.end()
    lgrass_9plants_plantfusion.end()


    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


if __name__ == "__main__":
    in_folder = "inputs_lgrass"
    genetic_model_folder = os.path.join(os.getcwd(), in_folder, "modelgenet")
    out_folder = "outputs/lgrass_default"
    write_geo = True
    graphs = True
    scenario_index = 0

    simulation(in_folder, genetic_model_folder, out_folder, scenario_index, write_geo, graphs)
