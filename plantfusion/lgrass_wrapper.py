import os
import shutil
import pandas
import numpy
import math
import itertools

import openalea.lpy as lpy

from lgrass import param_reproduction_functions as prf
from lgrass import meteo_ephem as meteo_ephem
from lgrass import cuts as cuts
from lgrass import gen_lstring as gen_lstring
from lgrass.output_data import CsvGenerator

from plantfusion.utils import create_child_folder
from plantfusion.indexer import Indexer
from plantfusion.planter import Planter
from plantfusion.light_wrapper import Light_wrapper


class Lgrass_wrapper:
    """Wrapper for l-grass model

    Note
    ---- 
    Only one lsystem per wrapper instance

    Parameters
    ----------
    name : str, optional
        name of the fspm instance, by default "lgrass"
    indexer : Indexer, optional
        indexer for listing FSPM in the simulation, by default Indexer()
    in_folder : str, optional
        input folder path, by default "inputs"
    out_folder : str, optional
        output folder path, by default None
    configuration_file : str, optional
        csv file name which contains parameters for one simulation, by default "plan_simulation.csv"
    id_scenario : int, optional
        which simulation to run in the configuration_file, by default 0
    caribu_parameters_file : str, optional
        csv file with parameters for the lighting post-process, by default "param_caribu"
    activate_genetic_model : bool, optional
        if genetic generation is activated, by default False
    genetic_model_folder : str, optional
        folder path which contains parameters and executable og the genetic model, by default "modelgenet"
    generation_index : int, optional
        which generation during the simulation, by default 1
    outputs_graphs : bool, optional
        if you wish to plot graphs from output results, by default False
    """  
    def __init__(
        self,
        name="lgrass",
        indexer=Indexer(),
        planter=Planter(),
        in_folder="inputs",
        out_folder=None,
        configuration_file="plan_simulation.csv",
        id_scenario=0,
        caribu_parameters_file="param_caribu",
        lai="plantfusion",
        number_of_plants=None,
        activate_genetic_model=False,
        genetic_model_folder="modelgenet",
        generation_index=1,
        outputs_graphs=False,
    ):
        """Constructor
        """        
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
        self.lai = lai

        # plan de simulation
        self.setup_list = pandas.read_csv(os.path.join(in_folder, configuration_file), sep=",")
        self.setup = self.setup_list.iloc[id_scenario]
        self.simulation_name = self.setup["name"]

        # initialisation du modèle génétique
        if activate_genetic_model:
            self.simulation_name = self.setup["name"] + "_G" + str(generation_index)

            self.genet_src = os.path.join(in_folder, "insim.txt")

            self.genet_dst = genetic_model_folder + "_" + self.simulation_name
            os.system(f"mkdir -p {self.genet_dst}")
            for f in os.listdir(genetic_model_folder):
                if os.name == "posix":
                    os.system(f"cp {os.path.join(genetic_model_folder, f)} {os.path.join(self.genet_dst, f)}")
                else:
                    os.system(f"copy {os.path.join(genetic_model_folder, f)} {os.path.join(self.genet_dst, f)}")

            if os.name == "posix":
                self.genet_exe = "simpraise"
            else:
                self.genet_exe = "simpraise.exe"

        # plants carbon data
        plants_carbon_data = os.path.join(in_folder, "liste_plantes.csv")

        # CARIBU additionnal parameters
        self.lighting_parameters = pandas.read_csv(
            os.path.join(in_folder, caribu_parameters_file + ".csv"), sep=";", header=0
        )
        self.lighting_parameters = dict(zip(self.lighting_parameters, self.lighting_parameters.iloc[0, :]))


        # option reproduction des plantes
        opt_repro = self.setup["option_reproduction"]
        in_genet_file = os.path.join(genetic_model_folder, "ped.r") if opt_repro in ("spikelets", "SPPR_2012") else None

        # init lsystem
        lpy_filename = os.path.join(os.path.dirname(prf.__file__), "lgrass.lpy")
        self.lsystem = lpy.Lsystem(lpy_filename)
        self.lsystem.name_sim = self.simulation_name

        ## Lsystem parameters
        # initialisation parameters
        (
            self.lsystem.ParamP,
            self.lsystem.nb_plantes,
            self.lsystem.NBlignes,
            self.lsystem.NBcolonnes,
            posPlante,
            self.lsystem.Plantes,
            self.lsystem.Genotypes,
            self.lsystem.flowering_model,
        ) = prf.define_param(
            in_param_file=plants_carbon_data,
            in_genet_file=in_genet_file,
            out_param_file=os.path.join(out_folder, "data", self.simulation_name + ".csv"),
            id_gener=generation_index,
            opt_repro=opt_repro,
            number_of_plants=number_of_plants,
        )
        
        # plant positions
        self.generation_type = planter.generation_type
        if planter.generation_type == "default":
            self.lsystem.posPlante = posPlante
            if len(posPlante) > 1:
                self.lsystem.Espacement = float(posPlante[1][1] - posPlante[0][1])
            else:
                self.lsystem.Espacement = 1.
        
        elif planter.generation_type == "random":
            self.lsystem.posPlante = planter.generate_random_lgrass(indice_instance=self.lgrass_index)

        if planter.save_plant_positions:
            for i, p in enumerate(posPlante):
                planter.plants_information.loc[len(planter.plants_information)] = [i, 
                                                                                    self.name, 
                                                                                    self.global_index, 
                                                                                    [float(p[0])*0.01, float(p[1])*0.01, 0.], 
                                                                                    []]

        # options and more parameters
        self.lsystem.option_tallage = self.setup["option_tallage"]
        self.lsystem.option_senescence = self.setup["option_senescence"]
        self.lsystem.option_floraison = self.setup["option_floraison"]
        self.lsystem.option_tiller_regression = self.setup["option_tiller_regression"]
        self.lsystem.option_morphogenetic_regulation_by_carbone = self.setup[
            "option_morphogenetic_regulation_by_carbone"
        ]
        self.lsystem.derivationLength = int(self.setup["derivationLength"])
        self.lsystem.sowing_date = self.setup["sowing_date"]
        self.lsystem.site = self.setup["site"]
        self.lsystem.meteo = meteo_ephem.import_meteo_data(
            self.setup["meteo_path"], self.setup["sowing_date"], self.setup["site"]
        )
        self.lsystem.OUTPUTS_DIRPATH = os.path.join(out_folder, "data")
        self.lsystem.GRAPHS_DIRPATH = os.path.join(out_folder, "graphs")
        self.lsystem.output_induction_file_name = self.simulation_name + "_" + "induction"
        self.lsystem.output_organ_lengths_file_name = self.simulation_name + "_" + "organ_lengths"
        self.lsystem.option_graphic_outputs = outputs_graphs

        # set display lsystem
        self.lsystem.display_only_leaves = True

        # external coupling mode
        self.lsystem.option_external_coupling = True

        # Gestion des tontes
        if self.setup["option_tontes"]:
            self.lsystem.cutting_dates, self.lsystem.derivationLength = cuts.define_cutting_dates(
                self.lsystem.meteo,
                int(self.setup["derivationLength"]),
                self.setup["cutting_freq"],
                self.setup["cutting_start"],
            )
        else:
            self.lsystem.cutting_dates = []

        # initalise lstring
        self.lstring = self.lsystem.axiom

        self.lsystem.current_day = 1

    def derive(self, t):
        """Derive the lsystem

        Parameters
        ----------
        t : int
            timestep
        """        
        self.lstring = self.lsystem.derive(self.lstring, t, 1)

    def light_inputs(self):
        """Return a geometric scene of plant leaves

        Returns
        -------
        plantgl.Scene 
            plantgl.Scene of the entire canopy, aerial organs
        """        
        return self.lsystem.sceneInterpretation(self.lstring)

    def light_results(self, lighting: Light_wrapper):
        """Interpret the lighting results, return the biomass production

        Note
        ----
        tiller regression can be computed but is not integrated in lgrass yet

        Parameters
        ----------
        lighting : Light_wrapper
            lighting results and parameters
        """        

        # filtering by vegetation type in the global simulation
        results = pandas.DataFrame({})
        out = lighting.results_organs()
        if isinstance(self.global_index, list):
            filter = out["VegetationType"].isin(self.global_index)
        else:
            filter = out.VegetationType == self.global_index
        results = out[filter]
        
        
        self.leaves_area = {}

        # biomass computation
        self.biomass_production = [0.] * self.lsystem.nb_plantes
        self.leaves_area = [0.] * self.lsystem.nb_plantes
        if not results.empty:
            for plant_id in range(self.lsystem.nb_plantes):
                # filter to only organ within the current plant
                list_organs = [i for i in results["Organ"] if self.lstring[i][0].id_plante == plant_id]
                df_plant = results[results["Organ"].isin(list_organs)]

                # radiation integered with organ surface
                energy_per_plant = (df_plant["par Ei"] * df_plant["Area"]).sum()
                self.biomass_production[plant_id] = energy_per_plant * self.lighting_parameters["RUE"] # Ray: MJ PAR ; RUE : g MJ-1

                self.leaves_area[plant_id] = df_plant["Area"].sum()

        # # tiller regression option (not returned yet)
        # if self.lsystem.option_tiller_regression:

        #     # time condition (always true ?)
        #     if self.lsystem.current_day - self.lighting_parameters['period_considered_tiller_regression'] <= self.lsystem.current_day:
        #         if len(self.lsystem.tiller_appearance) > 0:
        #             # plant loop
        #             for id_plante in numpy.unique(self.lsystem.tiller_appearance.id_plante):
        #                 plant_tillers = self.lsystem.tiller_appearance[self.lsystem.tiller_appearance.id_plante == id_plante]
        #                 youngest_tillers = plant_tillers[
        #                     plant_tillers.appearance_date == max(plant_tillers.appearance_date)]
        #                 youngest_tillers_radiations = pd.DataFrame()

        #                 # talle loop
        #                 for id_talle in youngest_tillers.id_talle:
        #                     # filter organ in the actual talle and plant
        #                     list_organs = [ i for i in results["Organ"] if self.lstring[i][0].id_plante == id_plante and self.lstring[i][0].id_talle == id_talle]
        #                     df_plant_talle = results[results["Organ"].isin(list_organs)]
        #                     tiller_raditation = (df_plant_talle["par Ei"] * df_plant_talle["Area"]).sum() / df_plant_talle["Area"].sum()

        #                     youngest_tillers_radiations = pd.concat([youngest_tillers_radiations,
        #                         pd.DataFrame([{'id_talle': [id_talle], 'Ei_tiller': [tiller_raditation]}])], ignore_index=True)

        #                 potential_tiller_to_remove = youngest_tillers_radiations[
        #                     youngest_tillers_radiations.Ei_tiller == min(youngest_tillers_radiations.Ei_tiller)]
        #                 if potential_tiller_to_remove.Ei_tiller.item() <= self.lighting_parameters['radiation_threshold']:
        #                     tiller_to_remove = pd.concat([tiller_to_remove, pd.DataFrame([
        #                         {'id_plante': [id_plante], 'id_talle': [potential_tiller_to_remove.id_talle.item()]}])], ignore_index=True)

    def run(self, planter):
        
        # Calcul LAI proximite
        if self.lai == "plantfusion":
            rapportS9_SSol_dict = self._proximity_LAI(planter)
        elif self.lai == "lgrass":
            rapportS9_SSol_dict = self._lgrass_default_lai(planter)
        
        # dependance avec biomprod (couplage lumière)
        Biomasse_racinaire, PourcentageFeuilGrowthRealized, PourcentageRootGrowthRealized = self._potential_plant_growth()

        # copie dans le lsystem
        self.lsystem.rapportS9_SSol_dict = rapportS9_SSol_dict
        self.lsystem.Biomasse_racinaire = Biomasse_racinaire
        self.lsystem.PourcentageFeuilGrowthRealized = PourcentageFeuilGrowthRealized
        self.lsystem.PourcentageRootGrowthRealized = PourcentageRootGrowthRealized
                

    def end(self):
        """Writes outputs and close the instance in the simulation
        It return a plants crossover matrix from genetic data
        """        
        # Matrice de croisement des plantes
        if self.setup["option_reproduction"] != "False":
            self.genet_mat = prf.create_seeds(
                self.lstring,
                self.lsystem.nb_plantes,
                self.lsystem.nb_talle,
                self.setup["option_reproduction"],
                self.setup["cutting_freq"],
                self.lsystem.ParamP,
            )
            numpy.savetxt(os.path.join(self.out_folder, "data", str(self.simulation_name) + "_mat.csv"), self.genet_mat)
        else:
            mat = 0

        # Sauvegarder la lstring dans un répertoire pour pouvoir la charger dans une prochaine simulation
        if self.setup["option_sauvegarde"]:
            gen_lstring.save_lstring(self.lstring, self.lsystem)

        csv_generator = CsvGenerator(self.lstring, self.simulation_name, os.path.join(self.out_folder, "data"))
        csv_generator.metadata_to_csv(self.lsystem)
        csv_generator.leaves_to_csv()
        csv_generator.internodes_to_csv()
        csv_generator.apex_to_csv()

        # move graph picuters created by lsystem in graphs folder
        if self.lsystem.option_graphic_outputs:
            for file in os.listdir(os.path.join(self.out_folder, "data")):
                if file.endswith(".png") or file.endswith(".PNG") or file.endswith(".pdf"):
                    shutil.move(
                        os.path.join(self.out_folder, "data", file), os.path.join(self.out_folder, "graphs", file)
                    )

        # Vider le lsystem
        self.lsystem.clear()
        print("".join((self.simulation_name, " - done")))

    def _proximity_LAI(self, planter):
        local_LAI_by_plant = {}

        for row in planter.plants_information.itertuples():

            if self.generation_type == "default":
                soil_surface = (len(row._5)+1) * (planter.scanning_ray/math.sqrt(2))**2
            else:
                soil_surface = math.pi * planter.scanning_ray**2

            local_LAI_by_plant[row.plant] = sum([self.leaves_area[i[0]] for i in row._5] + [self.leaves_area[row.plant]]) / soil_surface
        
        return local_LAI_by_plant
    
    def _lgrass_default_lai(self, planter): 
        def LAI_proximite(id_plante):
            def posTOnum(pos): #Determination de l'id de la plante a partir de la position
                id_plante = -1
                for i in range((self.lsystem.NBlignes)*(self.lsystem.NBcolonnes)):
                    if (self.lsystem.posPlante[i] == pos):
                        id_plante = i
                    return id_plante
            
            def numTOpos(id_plante): #Determination de la position de la plante a partir de l'id
                return self.lsystem.posPlante[id_plante]
            
            def deplacPos(pos, depl): #Selectionne la plante adjacente donnee
                ret=list(pos)  #Cree une copie de pos
                if ((depl[0]==-1 or depl[0]==0 or depl[0]==1) and (depl[1]==-1 or depl[1]==0 or depl[1]==1)): #Valeurs entre -1 et 1 : plantes adjacentes et non a 2 espacements
                    ret[0]=pos[0] + depl[0]
                    ret[1]=pos[1] + depl[1]
                    #Simulation d'un couvert de taille infini:
                    if (ret[0] == -1): ret[0] = self.lsystem.NBlignes-1
                    if (ret[0] == self.lsystem.NBlignes): ret[0] = 0
                    if (ret[1] == -1): ret[1] = self.lsystem.NBcolonnes-1
                    if (ret[1] == self.lsystem.NBcolonnes): ret[1] = 0
                return ret
            
            def adj9Plantes(pos): #Fait la liste des plantes adjacentes + plante cible
                position_plantes_adjacentes = []
                for p in itertools.product([-1,0,1], [-1,0,1]):
                    position_plantes_adjacentes.append(deplacPos(pos, p))
                return position_plantes_adjacentes
            
            def surfol9Plantes(id_plante): #Calcule la surface foliaires des 9 plantes adjacentes
                position_plantes_adjacentes = adj9Plantes(numTOpos(id_plante)) # liste des positions des plantes adjacentes
                surface_9_plantes = 0
                for p in position_plantes_adjacentes:
                    if posTOnum(p)>-1:
                        surface_9_plantes += self.lsystem.surface_foliaire_emergee[posTOnum(p)]
                return surface_9_plantes
            
            nb_plantes = 9
            LAI_proximite = surfol9Plantes(id_plante) / (50**2 * nb_plantes)
            return LAI_proximite
        
        local_LAI_by_plant = {}

        for row in planter.plants_information.itertuples():
            local_LAI_by_plant[row.plant] = LAI_proximite(row.plant)

        return local_LAI_by_plant
        
    def _potential_plant_growth(self):
        Biomasse_racinaire = self.lsystem.Biomasse_racinaire
        PourcentageFeuilGrowthRealized = self.lsystem.PourcentageFeuilGrowthRealized
        PourcentageRootGrowthRealized = self.lsystem.PourcentageRootGrowthRealized

        for id_plante in range(self.lsystem.nb_plantes):
            # Mise a jour biomasse et pourcentage echelle plante
            if self.lsystem.Biomasse_aerienne[id_plante]>0 and self.lsystem.t1[id_plante]*math.exp(self.lsystem.TPS*self.lsystem.t2[id_plante])-self.lsystem.t1[id_plante]>0 and self.lsystem.surface_foliaire_emergee[id_plante] > 0 and self.lsystem.Demande_feuille[id_plante] != 0: #Verifie qu'il y a une partie aerienne et que TPS, t1 et t2 sont positifs          # surface_foliaire_emergee[ID] > 0 and Demande_feuille[ID] > 0 : ajout SIMON
                # Calcul biomasse creee
                if self.lsystem.option_morphogenetic_regulation_by_carbone == True:
                    # Approche RUE
                    BiomassCree = self.biomass_production[id_plante]
                else:
                    # Approche croissance biomasse exponentielle (cf th?se Vincent Migault)
                    BiomassCree=self.lsystem.t2[id_plante]*(self.lsystem.Biomasse_aerienne[id_plante]+self.lsystem.t1[id_plante])*(1+self.lsystem.Beta[id_plante]*self.lsystem.Alpha[id_plante]*(self.lsystem.Biomasse_aerienne[id_plante])**(self.lsystem.Alpha[id_plante]-1))

                # Allocation aerien/racinaire + ratio offre/demande
                if BiomassCree>0 and 0.99*BiomassCree > self.lsystem.Demande_feuille[id_plante]:  # demande satisfaite
                    RootBiomassCree = BiomassCree - self.lsystem.Demande_feuille[id_plante]  # Allocation surplus non utilise par parties aeriennes aux racines
                    Biomasse_racinaire[id_plante] += RootBiomassCree
                    PourcentageFeuilGrowthRealized[id_plante] = 1  # Pas de limitation croissance foliaire pour prochain pas de temps
                else:  # Demande non satisafite
                    RootBiomassCree = 0.01*BiomassCree  # Racines ne prennent qu'1% de l'offre
                    Biomasse_racinaire[id_plante] += RootBiomassCree
                    BiomassFeuilCree = 0.99*BiomassCree
                    PourcentageFeuilGrowthRealized[id_plante] = max(0, BiomassFeuilCree / self.lsystem.Demande_feuille[id_plante]) # Reduction croissance foliaire potentielle au prorata offre/demande pour prochain pas de temps
            else:
                RootBiomassCree=0
                BiomassCree=0

            
            if self.lsystem.RootPotentialNewBiomass[id_plante]==0:
                PourcentageRootGrowthRealized[id_plante]=0
                
            else:
                PourcentageRootGrowthRealized[id_plante]=RootBiomassCree/self.lsystem.RootPotentialNewBiomass[id_plante]
                if PourcentageRootGrowthRealized[id_plante]>1:
                    PourcentageRootGrowthRealized[id_plante]=1
                else:
                    if self.lsystem.RootPotentialNewBiomass[id_plante]-RootBiomassCree < self.lsystem.Reserve[id_plante]:
                        PourcentageRootGrowthRealized[id_plante]=1
                
                    else:
                        PourcentageRootGrowthRealized[id_plante]=(RootBiomassCree+self.lsystem.Reserve[id_plante])/self.lsystem.RootPotentialNewBiomass[id_plante]

        return Biomasse_racinaire, PourcentageFeuilGrowthRealized, PourcentageRootGrowthRealized

    def doy(self):
        """Current thermal day

        Returns
        -------
        int
            current day in the simulation
        """        
        return self.lsystem.current_day

    def energy(self):
        """Light energy from meteo input during current day

        Returns
        -------
        float
            light energy
        """        
        return self.lsystem.meteo[self.lsystem.meteo.experimental_day == self.lsystem.current_day].PAR_incident.iloc[0]


def lgrass_soil_domain(espacement=50, rows=1, columns=1):
    """Soil domain for computing infinite pattern

    Parameters
    ----------
    espacement : int, optional
        length of one side in cm, by default 50
    rows : int, optional
        number of plant rows, by default 1
    columns : int, optional
        number of plant columns, by default 1

    Returns
    -------
    tuple of tuple
        ((xmin, ymin), (xmax, ymax))
    """   
    # s = [[-espacement / 2, -espacement / 2],
    #     [espacement * (rows - 1) + espacement / 2, espacement * (columns - 1) + espacement / 2]]

    return ((0., 0.), (0.01 * espacement * (columns), 0.01 * espacement * (rows)))
