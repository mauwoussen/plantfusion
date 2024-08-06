import os
import pandas as pd
import numpy

import lgrass
import openalea.lpy as lpy

from plantfusion.utils import create_child_folder
from plantfusion.indexer import Indexer


class Lgrass_wrapper:
    def __init__(self, name="lgrass", 
                 indexer=Indexer(), 
                 in_folder="inputs", 
                 out_folder=None,
                 plan_sim_file="plan_simulation.csv",
                 param_plant_file = "Parametre_plante_Lgrass.csv",# un fichier de remplacement du modèle génétique qui génère une population de C et détermine le nombre de plantes du couvert        
                 id_scenario = 0,
                 ):
        if out_folder is not None:
            try:
                os.mkdir(os.path.normpath(out_folder))
                print("Directory ", out_folder, " Created ")
            except FileExistsError:
                pass

            # output folder for l-egume
            out_folder = os.path.normpath(out_folder)
            self.out_folder = os.path.join(out_folder, name)
            try:
                os.mkdir(os.path.normpath(self.out_folder))
                print("Directory ", self.out_folder, " Created ")
            except FileExistsError:
                pass
            create_child_folder(self.out_folder, "brut")
            create_child_folder(self.out_folder, "graphs")

        else:
            self.out_folder = ""

        self.name = name
        self.indexer = indexer
        self.global_index = indexer.global_order.index(name)
        self.lgrass_index = indexer.lgrass_names.index(name)

        # paramètres plantes
        path_param = os.path.join(in_folder, "Parametre_plante_Lgrass.xls")

        onglet1 = "FL"
        onglet2 = "FC"
        TableParamP1 = pd.read_excel(path_param, sheet_name=onglet1)
        TableParamP2 = pd.read_excel(path_param, sheet_name=onglet2)
        self.paramP1 = dict(zip(TableParamP1["name"], TableParamP1["value"]))
        self.paramP2 = dict(zip(TableParamP2["name"], TableParamP2["value"]))

        # plan de simulation
        self.plan_sim = pd.read_csv(os.path.join(in_folder, plan_sim_file), sep=',')
        self.row = self.plan_sim.iloc[id_scenario]
        # données carbone

        # other data

        # init flowering_param

        # init lsystem
        self.lpy_filename = os.path.join(lgrass.__path__[0], 'lgrass.lpy')
        self.lsystem = lpy.Lsystem(self.lpy_filename)
        self.lsystem.name_sim = name
        self.lstring = self.lsystem.axiom

        # transmet les options au lsystem

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
