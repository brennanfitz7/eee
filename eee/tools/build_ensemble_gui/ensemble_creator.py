import eee

from eee.tools.build_ensemble_gui.base import VariableWidgetStack
from eee.tools.build_ensemble_gui.species import SpeciesWidget
from eee.tools.build_ensemble_gui.titration import TitrationWidget
from eee.tools.build_ensemble_gui.fitness import FitnessWidget
from eee.tools.build_ensemble_gui.basic_info import BasicInfoWidget

from eee.tools import plots

import ipywidgets as widgets
from matplotlib import pyplot as plt

import copy

class EnsembleCreator:
    """
    Create an ensemble.
    """
    
    def __init__(self,update_callback):
        
        self._update_callback = update_callback

        w = []
        w.append(widgets.HTML(value="<h1>Create an ensemble</h1>"))
        w.append(widgets.HTML(value="<b><hr/></b>"))

        # Basic ensemble information
        self._basic_info = BasicInfoWidget(update_callback=self._update_follower)
        w.append(self._basic_info.widget)
        w.append(widgets.HTML(value="<hr/>"))
        
        # Species builder
        w.append(widgets.HTML(value="<h3>Define species</h3>"))
        self._species = VariableWidgetStack(update_callback=self._update_follower,
                                            widget_to_stack=SpeciesWidget,
                                            button_description="Add species")
        w.append(self._species.widget)
        w.append(widgets.HTML(value="<hr/>"))
        
        # Titration definition
        w.append(widgets.HTML(value="<h3>Define ligand titration</h3>"))
        self._titration = TitrationWidget(update_callback=self._update_follower)
        w.append(self._titration.widget)
        w.append(widgets.HTML(value="<hr/>"))

        # # Fitness definition
        w.append(widgets.HTML(value="<h3>Define selection conditions</h3>"))
        self._fitness = VariableWidgetStack(update_callback=self._update_follower,
                                            widget_to_stack=FitnessWidget,
                                            button_description="Add condition")
        w.append(self._fitness.widget)
        w.append(widgets.HTML(value="<hr/>"))

        # Display panel
        self._out = widgets.widget_output.Output()
        w.append(self._out)
        
        self._w = w
        self._main = widgets.VBox(self._w)

        self._load_defaults()

    def _update_follower(self,*args,**kwargs):
                
        values = self.get_values()
        self._update_callback(self.get_values())
        
        # Start creating an ensemble
        ens = eee.core.Ensemble(gas_constant=values["gas_constant"])
        
        # Try to add species to the ensemble. If this throws an error, the
        # species info is not complete. Set good_ensemble to False
        good_ensemble = True
        for species in values["species"]:
            tmp_species = copy.deepcopy(species)
            ligand_stoich = {}
            if "ligands" in tmp_species:
                for lig in tmp_species["ligands"]:
                    ligand_stoich[lig[0]] = lig[1]
                tmp_species.pop("ligands")
            
            tmp_species["ligand_stoich"] = ligand_stoich
            try:
                ens.add_species(**tmp_species)
            except ValueError:
                good_ensemble = False
                break
            
        # Make sure we can calculate an observable from this ensemble
        if good_ensemble:
            try:
                ens.get_obs()
            except ValueError:
                good_ensemble = False
                
        # Only draw titration of ensemble if we have titration X values
        ligand_dict = {}
        if good_ensemble:
            if "titration" not in values:
                good_ensemble = False
            else:
                ligand_dict = values["titration"]
                
        with self._out:

            self._out.clear_output()
            if len(values["species"]) == 0:
                return
            
            if good_ensemble:

                fig, axes = plt.subplots(1,2,figsize=(12,6))
                fig, axes[0] = plots.energy_diagram(values,
                                                    fig=fig,
                                                    ax=axes[0],
                                                    context_manager=self._out)
                fig, axes[1] = plots.titration(ens=ens,
                                               ligand_dict=ligand_dict,
                                               fig=fig,
                                               ax=axes[1],
                                               context_manager=self._out)
            else:
                fig, ax = plt.subplots(1,figsize=(6,6))
                fig, ax = plots.energy_diagram(values,
                                               fig=fig,
                                               ax=ax,
                                               context_manager=self._out)
                
            fig.tight_layout()
            plt.show()
        
        if good_ensemble:
            with self._out:
                print(ens.to_dict())

    def get_values(self,*args,**kwargs):
        
        out = self._basic_info.get_values()
        species = self._species.get_values()
        
        # Get all available ligands. 
        available_ligands = []
        species = self._species.get_values()
        for s in species:
            if "ligands" in s:
                for lig in s["ligands"]:
                    available_ligands.append(lig[0])
        available_ligands = list(set(available_ligands))
        available_ligands.sort()
        
        self._titration.update(current_ligands=available_ligands)
        for f in self._fitness.widgets:
            f.update(current_ligands=available_ligands)
    
        ligand_titrations = self._titration.get_values()
        
        out["species"] = species
        out["titration"]= ligand_titrations
               
        return out

    def _load_defaults(self):
        
        self._species.add_widget(species_name="A",
                                  dG0=0,
                                  observable=True,
                                  folded=True)
        self._species.add_widget(species_name="B",
                                  dG0=10,
                                  observable=False,
                                  folded=True,
                                  ligands=[("X",1)])
        self._species.add_widget(species_name="U",
                                  dG0=10,
                                  observable=False,
                                  folded=False)

        
    def main(self):
        return self._main