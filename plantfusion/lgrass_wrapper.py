import os
import pandas as pd


class Lgrass_wrapper:
    def __init__(self, name="lgrass", indexer=Indexer(), in_folder="inputs", simulation_name="", out_folder=None):

        self.name = name

        # paramètres plantes
        path_param = os.path.join(in_folder, "Parametre_plante_Lgrass.xls")

        onglet1 = "FL"
        onglet2 = "FC"
        TableParamP1 = pd.read_excel(path_param, sheet_name=onglet1)
        TableParamP2 = pd.read_excel(path_param, sheet_name=onglet2)
        self.paramP1 = dict(zip(TableParamP1["name"], TableParamP1["value"]))
        self.paramP2 = dict(zip(TableParamP2["name"], TableParamP2["value"]))

        # plan de simulation

        # données carbone

        # other data

        # init flowering_param

        # init lsystem

        # trnasmet les options au lsystem

    def light_inputs(self):
        return self.scene

    def light_results(self, lighting):
        self.energy = lighting.results_organs().mean()["par Eabs"]

    def soil_inputs(self, soil_dimensions):
        nb_plants = 1

        N_content_roots_per_plant = [0.5] * nb_plants
        plants_light_interception = [0.4] * nb_plants

        roots_length = 6.0  # m
        roots_length_per_plant_per_soil_layer = []
        for i in range(nb_plants):
            # on répartit de manière homogène les racines à travers les couches du sol
            # convertit m en cm # --> peut etre en metre finalement
            rootLen_i = numpy.ones(soil_dimensions) * roots_length / numpy.prod(soil_dimensions)
            roots_length_per_plant_per_soil_layer.append(rootLen_i)

        return (
            N_content_roots_per_plant,
            roots_length_per_plant_per_soil_layer,
            [self.soil_parameters] * nb_plants,
            plants_light_interception,
        )

    def soil_results(self, uptakeN):
        pass

    def run(self):
        pass

    def end(self):
        self.lsystem.clear()
        print("--- END ---")
