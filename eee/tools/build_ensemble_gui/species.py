from eee.tools.build_ensemble_gui.base import MetaWidget
from eee.tools.build_ensemble_gui.base import VariableWidgetStack
from eee.tools.build_ensemble_gui.ligands import LigandWidget

import ipywidgets as widgets

class SpeciesWidget(MetaWidget):

    def __init__(self,
                 update_callback,
                 species_name="",
                 dG0=0,
                 observable=True,
                 folded=True,
                 ligands=None):
    
        super().__init__(update_callback=update_callback)

        self.build_widget(species_name=species_name,
                          dG0=dG0,
                          observable=observable,
                          folded=folded,
                          ligands=ligands)
        
    def build_widget(self,
                     species_name="",
                     dG0=0,
                     observable=True,
                     folded=True,
                     ligands=None):
        
        left_items = []
        left_items.append(widgets.Text(value=species_name,
                                       placeholder="species name",
                                       continuous_update=False))
        left_items.append(widgets.FloatText(dG0,
                                            description="dG0:",
                                            continuous_update=False))

        right_items = []
        right_items.append(widgets.Checkbox(observable,description="observable"))
        right_items.append(widgets.Checkbox(folded,description="folded"))

        # Watch these widgets
        left_items[0].observe(self._watcher)
        left_items[1].observe(self._watcher)
        right_items[0].observe(self._watcher)
        right_items[1].observe(self._watcher)
        
        top_box = widgets.HBox([widgets.VBox(left_items),
                                widgets.VBox(right_items)])
        
        # Ligand for building ligands
        self._ligands = VariableWidgetStack(update_callback=self._watcher,
                                            widget_to_stack=LigandWidget,
                                            button_description="Add ligand")

        if ligands is not None:
            for lig in ligands:
                self._ligands.add_widget(ligand_name=lig[0],
                                         stoichiometry=lig[1])
        
        self._widget = widgets.VBox([top_box,
                                     self._ligands.widget])
        
    def get_values(self):

        if self._widget is None:
            return None

        fields = self._widget.children[0]
    
        v = {}
        v["name"] = fields.children[0].children[0].get_interact_value()
        v["dG0"] = fields.children[0].children[1].get_interact_value()
        v["observable"] = fields.children[1].children[0].get_interact_value()
        v["folded"] = fields.children[1].children[1].get_interact_value()
        
        v["ligands"] = self._ligands.get_values()

        return v
   
        


