import eee
from eee.tools.build_ensemble_gui.base import MetaWidget
from eee.core.data import GAS_CONSTANT

import ipywidgets as widgets
import numpy as np

import inspect

class BasicInfoWidget(MetaWidget):

    def __init__(self,
                 update_callback=None):
        
        super().__init__(update_callback=update_callback)
        self.build_widget()

    def build_widget(self):

        # temperature
        ens = eee.core.Ensemble()
        p = inspect.signature(ens.get_obs).parameters
        default_temperature = p["temperature"].default
        
        temperature = widgets.BoundedFloatText(default_temperature,
                                               description="T:",
                                               continuous_update=False,
                                               min=np.nextafter(0, 1),
                                               max=np.finfo('d').max)
        temperature.observe(self._watcher)
        
        # Gas constant
        
        gas_constant = widgets.BoundedFloatText(GAS_CONSTANT,
                                                description="R:",
                                                continuous_update=False,
                                                min=np.nextafter(0, 1),
                                                max=np.finfo('d').max)
        gas_constant.observe(self._watcher)
        
        self._widget = widgets.HBox([temperature,gas_constant])

    def get_values(self):
        
        temperature = self._widget.children[0].get_interact_value()
        gas_constant = self._widget.children[1].get_interact_value()

        return {"temperature":temperature,
                "gas_constant":gas_constant}
