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
    
    plant_density = {plants_name : 400}
    xy_square_length = 0.4 # m
    planter = Planter(generation_type="random", indexer=index_log, plant_density=plant_density, xy_square_length=xy_square_length, save_plant_positions=True)

    lgrass = Lgrass_wrapper(
        name=plants_name,
        indexer=index_log,
        planter=planter,
        in_folder=in_folder,
        out_folder=out_folder,
        id_scenario=id_scenario,
        number_of_plants=planter.number_of_plants[0],
        activate_genetic_model=True,
        genetic_model_folder=genetic_model_folder,
        outputs_graphs=graphs,
    )

    lighting = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        sky="turtle46",
        planter=planter,
        indexer=index_log,
        writegeo=write_geo,
    )

    current_time_of_the_system = time.time()
    prf.rungenet(lgrass.genet_src, lgrass.genet_dst, lgrass.genet_exe, None, 0)

    # genetic generation loop
    for i in range(1, lgrass.setup["num_gener"] + 1):
        lgrass = Lgrass_wrapper(
            name=plants_name,
            indexer=index_log,
            planter=planter,
            in_folder=in_folder,
            out_folder=out_folder,
            id_scenario=id_scenario,
            number_of_plants=planter.number_of_plants[0],
            activate_genetic_model=True,
            genetic_model_folder=genetic_model_folder,
            outputs_graphs=graphs,
        )

        scanning_ray = 0.0142
        planter.scan_nearest_plants_neighours(scanning_ray)

        # daily loop
        for t in range(0, lgrass.lsystem.derivationLength):
            thermal_day = lgrass.lsystem.current_day
            lgrass.derive(t)

            # execute CARIBU on each new thermal day
            if thermal_day < lgrass.lsystem.current_day or t == 0:
                if lgrass.setup["option_morphogenetic_regulation_by_carbone"]:
                    scene_lgrass = lgrass.light_inputs()
                    lighting.run(energy=lgrass.energy(), scenes=[scene_lgrass])
                    lgrass.light_results(lighting=lighting)

            # compute local LAI
            lgrass.run(planter)

        lgrass.end()

        if lgrass.setup["option_reproduction"] != "False":
            prf.rungenet(lgrass.genet_src, lgrass.genet_dst, lgrass.genet_exe, lgrass.genet_mat, 1)

    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))


if __name__ == "__main__":
    in_folder = "inputs_lgrass"
    genetic_model_folder = os.path.join(os.getcwd(), in_folder, "modelgenet")
    out_folder = "outputs/lgrass_random"
    write_geo = True
    graphs = True
    scenario_index = 0

    simulation(in_folder, genetic_model_folder, out_folder, scenario_index, write_geo, graphs)
