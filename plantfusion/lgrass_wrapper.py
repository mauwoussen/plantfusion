import os
import shutil
import pandas as pd
import numpy

import openalea.lpy as lpy

from lgrass import param_reproduction_functions as prf
from lgrass import meteo_ephem as meteo_ephem
from lgrass import cuts as cuts
from lgrass import gen_lstring as gen_lstring
from lgrass.output_data import CsvGenerator

from plantfusion.utils import create_child_folder
from plantfusion.indexer import Indexer


class Lgrass_wrapper:
    def __init__(self, name="lgrass", 
                 indexer=Indexer(), 
                 in_folder="inputs", 
                 out_folder=None,
                 configuration_file="plan_simulation.csv",
                 id_scenario = 0,
                 activate_genetic_model = False,
                 genetic_model_folder = "modelgenet",
                 generation_index = 1,
                 outputs_graphs = False
                 ):
        # create the output folders
        if out_folder is not None:
            self.out_folder = out_folder
            try:
                os.mkdir(os.path.normpath(out_folder))
                print("Directory ", out_folder, " Created ")
            except FileExistsError:
                pass

            create_child_folder(self.out_folder, "data")
            create_child_folder(self.out_folder, "graphs")

        else:
            self.out_folder = ""

        self.name = name
        self.indexer = indexer
        self.global_index = indexer.global_order.index(name)
        self.lgrass_index = indexer.lgrass_names.index(name)

        # plan de simulation
        self.setup_list = pd.read_csv(os.path.join(in_folder, configuration_file), sep=',')
        self.setup = self.setup_list.iloc[id_scenario]
        self.simulation_name = self.setup["name"]

        # initialisation du modèle génétique
        if activate_genetic_model:
            self.simulation_name =  self.setup["name"] + "_G" + str(generation_index) 

            self.genet_src = os.path.join(in_folder, 'insim.txt')
            self.genet_dst = genetic_model_folder
            if os.name == "posix" :
                self.genet_exe = 'simpraise'    
            else:
                self.genet_exe = 'simpraise.exe'

        # plants carbon data
        plants_carbon_data = os.path.join(in_folder, 'liste_plantes.csv')

        # option reproduction des plantes
        opt_repro = self.setup["option_reproduction"]
        in_genet_file = os.path.join(genetic_model_folder, 'ped.r') if opt_repro in ('spikelets', 'SPPR_2012') else None

        # init lsystem    
        lpy_filename = os.path.join(os.path.dirname(prf.__file__), 'lgrass.lpy')
        self.lsystem = lpy.Lsystem(lpy_filename)
        self.lsystem.name_sim = self.simulation_name

        ## Lsystem parameters
        # initialisation parameters
        self.lsystem.ParamP,\
        self.lsystem.nb_plantes,\
        self.lsystem.NBlignes,\
        self.lsystem.NBcolonnes,\
        self.lsystem.posPlante,\
        self.lsystem.Plantes, \
        self.lsystem.Genotypes,\
        self.lsystem.flowering_model = prf.define_param(in_param_file=plants_carbon_data, 
                                                        in_genet_file=in_genet_file,
                                                        out_param_file=os.path.join(out_folder, "data", self.simulation_name + '.csv'), 
                                                        id_gener=generation_index, 
                                                        opt_repro=opt_repro)
        
        # options and more parameters
        self.lsystem.option_tallage = self.setup["option_tallage"]
        self.lsystem.option_senescence = self.setup["option_senescence"]
        self.lsystem.option_floraison = self.setup["option_floraison"]
        self.lsystem.option_tiller_regression = self.setup["option_tiller_regression"]
        self.lsystem.option_morphogenetic_regulation_by_carbone = self.setup["option_morphogenetic_regulation_by_carbone"]
        self.lsystem.derivationLength = int(self.setup["derivationLength"])
        self.lsystem.sowing_date = self.setup["sowing_date"]
        self.lsystem.site = self.setup["site"]
        self.lsystem.meteo = meteo_ephem.import_meteo_data(self.setup["meteo_path"], self.setup['sowing_date'], self.setup['site'])
        self.lsystem.OUTPUTS_DIRPATH = os.path.join(out_folder, "data")
        self.lsystem.GRAPHS_DIRPATH = os.path.join(out_folder, "graphs")
        self.lsystem.output_induction_file_name = self.simulation_name + '_' + 'induction'
        self.lsystem.output_organ_lengths_file_name = self.simulation_name + '_' + 'organ_lengths'
        self.lsystem.option_graphic_outputs = outputs_graphs
        
        # Gestion des tontes
        if self.setup["option_tontes"]:
            self.lsystem.cutting_dates, self.lsystem.derivationLength = cuts.define_cutting_dates(self.lsystem.meteo,
                                                                                        int(self.setup["derivationLength"]),
                                                                                        self.setup["cutting_freq"],
                                                                                        self.setup["cutting_start"])
        else:
            self.lsystem.cutting_dates = []

        # initalise lstring
        self.lstring = self.lsystem.axiom

        # Gestion caribu
        # opt_caribu = self.setup["option_caribu"]
        # if opt_caribu:
        #     dico_caribu = run_caribu_lgrass.init(meteo=self.lsystem.meteo, nb_plantes=self.lsystem.nb_plantes, scenario=self.setup)
        #     self.lsystem.BiomProd = [0.] * self.lsystem.nb_plantes
        #     # Rédaction d'un fichier de sortie
        #     path_out = os.path.join(out_folder, name + '_caribu.csv')
        #     output = open(path_out, 'w')
        #     output.write("GDD;Date;Day;nb_talles;biomasse_aerienne;surface_foliaire;lstring" + "\n")

    def derive(self, t):
        self.lstring = self.lsystem.derive(self.lstring, t, 1)

    def light_inputs(self):
        return self.lsystem.sceneInterpretation(self.lstring)

    def light_results(self, lighting):
        self.energy = lighting.results_organs().mean()["par Eabs"]
        pass

    def run(self):
        pass

    def end(self):
        # Matrice de croisement des plantes
        if self.setup["option_reproduction"] != "False":
            self.genet_mat = prf.create_seeds(self.lstring, 
                                                self.lsystem.nb_plantes, 
                                                self.lsystem.nb_talle, 
                                                self.setup["option_reproduction"], 
                                                self.setup["cutting_freq"], 
                                                self.lsystem.ParamP)
            numpy.savetxt(os.path.join(self.out_folder, "data", str(self.simulation_name)  + "_mat.csv"), self.genet_mat)
        else:
            mat = 0

        # Sauvegarder la lstring dans un répertoire pour pouvoir la charger dans une prochaine simulation
        if self.setup['option_sauvegarde']:
            gen_lstring.save_lstring(self.lstring, self.lsystem)

        csv_generator = CsvGenerator(self.lstring, self.simulation_name, os.path.join(self.out_folder, "data"))
        csv_generator.metadata_to_csv(self.lsystem)
        csv_generator.leaves_to_csv()
        csv_generator.internodes_to_csv()
        csv_generator.apex_to_csv()

        # move graph picuters created by lsystem in graphs folder
        if self.lsystem.option_graphic_outputs :
            for file in os.listdir(os.path.join(self.out_folder, "data")):
                if file.endswith(".png") or file.endswith(".PNG") or file.endswith(".pdf") :
                    shutil.move(os.path.join(self.out_folder, "data", file), 
                                os.path.join(self.out_folder, "graphs", file))

        # Vider le lsystem
        self.lsystem.clear()
        print(''.join((self.simulation_name, " - done")))

    def energy(self):
        return self.lsystem.meteo[self.lsystem.meteo.experimental_day == self.lsystem.current_day].PAR_incident.iloc[0]
