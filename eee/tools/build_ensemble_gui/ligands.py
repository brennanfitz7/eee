
from eee.tools.build_ensemble_gui.base import MetaWidget

import ipywidgets as widgets

class LigandWidget(MetaWidget):

    def __init__(self,
                 update_callback,
                 ligand_name="ligand",
                 stoichiometry=1.0):
        
        super().__init__(update_callback=update_callback)
        self.build_widget(ligand_name=ligand_name,
                          stoichiometry=stoichiometry)

    def build_widget(self,
                     ligand_name="ligand",
                     stoichiometry=1.0):

        # Create set of fields
        left_field = widgets.Text(value=ligand_name,
                                  placeholder="ligand name",
                                  continuous_update=False)
        right_field = widgets.BoundedFloatText(stoichiometry,
                                               description="stoichiometry",
                                               continuous_update=False,
                                               min=0)
                
        # Make sure we are monitoring both fields
        left_field.observe(self._watcher)
        right_field.observe(self._watcher)
        
        # Add with a remove button
        self._widget = widgets.HBox([left_field,
                                     right_field])
        
    def get_values(self):
        
        if self._widget is None:
            return None

        lig = self._widget
        ligand_name = lig.children[0].get_interact_value()
        ligand_stoich = lig.children[1].get_interact_value()

        return (ligand_name,ligand_stoich)
    