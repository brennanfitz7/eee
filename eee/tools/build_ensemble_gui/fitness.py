import eee
from eee.tools.build_ensemble_gui.base import MetaWidget

import ipywidgets as widgets

class FitnessWidget(MetaWidget):

    def __init__(self,
                 update_callback,
                 ligand_dict=None,
                 fitness_fcn=None,
                 select_on="fx_obs",
                 select_on_folded=True):
            
        super().__init__(update_callback=update_callback)

        self.build_widget(ligand_dict=ligand_dict,
                          fitness_fcn=fitness_fcn,
                          select_on=select_on,
                          select_on_folded=select_on_folded)
    
    def build_widget(self,
                     ligand_dict=None,
                     fitness_fcn=None,
                     select_on="fx_obs",
                     select_on_folded=True):
        
        if ligand_dict is None:
            ligand_dict = {}

        self._conditions_w = {}
        for k in ligand_dict:
            self._conditions_w[k] = widgets.FloatText(value=ligand_dict[k],
                                                      description=f"{k}:",
                                                      continuous_update=False)
            self._conditions_w[k].observe(self._watcher)

        # select_on hard-coded in as fx_obs. all we have at this point. 
        self._select_on_w = widgets.Select(options=["fx_obs"],
                                           value="fx_obs",
                                           description="select on:")
        self._select_on_w.observe(self._watcher)
        
        # Select on folded at this condition
        self._on_folded_w = widgets.Checkbox(value=select_on_folded,
                                             description="select on folded")
        self._on_folded_w.observe(self._watcher)

        # Get list of possible fitness functions
        ff_options = list(eee.core.FF_AVAILABLE.keys())
        if fitness_fcn is None:
            fitness_fcn = ff_options[0]
        self._ff_w = widgets.Select(options=ff_options,
                                    value=fitness_fcn,
                                    description="select for:")
        self._ff_w.observe(self._watcher)
        
        # Get threshold (for fitness function, if needed)
        self._threshold_w = widgets.FloatText(value=0,
                                              description="threshold",
                                              continuous_update=False)
        self._threshold_w.observe(self._watcher)
        
        cond_box = widgets.VBox(list(self._conditions_w.values()))
        select_box = widgets.HBox([self._select_on_w,
                                   self._on_folded_w,
                                   self._ff_w,
                                   self._threshold_w])
        
        self._widget = widgets.VBox([cond_box,
                                     select_box])

    def update(self,current_ligands):
        """
        Update the widget programatically. Key purpose is to synchronize the
        widget with the currently defined ligands in the species. 
        """
    
        interface_ligands = set(self._conditions_w.keys())
        current_ligands = set(current_ligands)

        # interface and current ligands match; do not update. 
        if interface_ligands == current_ligands:
            return None
        
        condition_box = self._widget.children[0]
        list_of_widgets = list(condition_box.children)

        # If a ligand is in current ligands, keep it. If a ligand is in the 
        # interface but not in current ligands, record that it should be removed.
        conditions_w = {}
        to_remove = []
        for i, k in enumerate(self._conditions_w.keys()):

            if k in current_ligands:
                conditions_w[k] = self._conditions_w[k]
            else:
                to_remove.append(i)
        
        # Go from highest to lowest in to_remove, popping extra ligands from 
        # the interface
        to_remove = to_remove[::-1]
        for r in to_remove:
            list_of_widgets.pop(r)
            
        # Get list of ligands missing from the interface
        missing_from_interface = list(current_ligands - interface_ligands)

        # Create new widgets for ligands missing from the interface
        for k in missing_from_interface:
            conditions_w[k] = widgets.FloatText(value=0,
                                                description=f"{k}:",
                                                continuous_update=False)
            conditions_w[k].observe(self._watcher)
            list_of_widgets.append(conditions_w[k])

        # Record new fitness conditions dictionary into self
        self._conditions_w = conditions_w

        # Update ligand widgets
        self._widget.children[0].children[0].children = tuple(list_of_widgets)

    def get_values(self):

        cond_box = self._widget.children[0]
        select_box = self._widget.children[1]

        ligand_dict = {}
        for cond in cond_box.children:
            k = cond.description.strip()[:-1]
            v = cond.get_interact_value()
            ligand_dict[k] = v
        
        out = {}
        out["select_on"] = select_box.children[0].get_interact_value()
        out["select_on_folded"] = select_box.children[1].get_interact_value()
        out["fitness_fcn"] = select_box.children[2].get_interact_value()
        out["fitness_kwargs"] = {"threshold":select_box.children[3].get_interact_value()}
        out["ligand_dict"] = ligand_dict

        return
            
            

