
from eee.tools.build_ensemble_gui.base import MetaWidget

import numpy as np
import ipywidgets as widgets


class TitrationWidget(MetaWidget):
    """
    Widget holding subwidgets to define titration behavior.
    """
    
    def __init__(self,
                 update_callback=None,
                 min_value=0,
                 max_value=20,
                 num_steps=100):
        """
        Create the widget.

        Parameters
        ----------
        update_callback : callable, optional
            function called when a widget element changes. should take the 
            output of get_values
        min_value : float, default=0
            default minimum value
        max_value : float, default=20
            default maximum value
        num_steps : int, default=100
            default number of steps in titration
        """
        
        super().__init__(update_callback=update_callback)
        self.build_widget(min_value=min_value,
                          max_value=max_value,
                          num_steps=num_steps)
        
    def build_widget(self,
                     min_value=0,
                     max_value=20,
                     num_steps=100):
        
        dd = widgets.Select(options=[],
                            value=None,
                            description="Ligand:",
                            disabled=True)

        min_value = widgets.FloatText(min_value,
                                      description="min",
                                      continuous_update=False)
        max_value = widgets.FloatText(max_value,
                                      description="max",
                                      continuous_update=False)

        num_steps = widgets.BoundedIntText(num_steps,
                                           description="num steps",
                                           continuous_update=False,
                                           min=2,
                                           max=np.iinfo("i").max)

        dd.observe(self._watcher)
        min_value.observe(self._watcher)
        max_value.observe(self._watcher)
        num_steps.observe(self._watcher)

        self._widget = widgets.HBox([dd,min_value,max_value,num_steps])
        
    
    def update(self,current_ligands):
        """
        Update the widget programatically. Key purpose is to synchronize the
        widget with the currently defined ligands in the species. 
        """
        
        # select box
        select_box = self._widget.children[0]
        
        # Change available options
        self._silent_update(select_box,"options",current_ligands)
        
        # Synchronize value with ligands available. 
        if select_box.value not in current_ligands:
            if len(current_ligands) == 0:
                self._silent_update(select_box,"value",None)
            else:
                self._silent_update(select_box,"value",select_box.options[0])
                    
        # Decide what to do with current ligands
        if len(current_ligands) == 0:
            self._silent_update(select_box,"disabled",True)
        else:
            self._silent_update(select_box,"disabled",False)
    
    def get_values(self):
        """
        Get the titration of the ligand as a dictionary (key is ligand name,
        value is numpy array of titration values)
        """
        
        if self._widget is None:
            return None

        w = self._widget
        
        ligand = w.children[0].get_interact_value()
        if ligand is None:
            return {}
        
        min_value = w.children[1].get_interact_value()
        max_value = w.children[2].get_interact_value()
        num_steps = w.children[3].get_interact_value()
        
        ligand_dict = {ligand:np.linspace(start=min_value,
                                          stop=max_value,
                                          num=num_steps)}
        return ligand_dict
    
    @property
    def widget(self):
        return self._widget