
import ipywidgets as widgets

class MetaWidget:
    """
    Class holding multiple widgets that observes them and returns their values 
    in a sane, modular way.
    """

    def __init__(self,update_callback=None):

        self._widget = None
        self._update_callback = update_callback

    def build_widget(self,*args,**kwargs):
        """
        Redefine this in a subclass. Should build the block of widgets and store
        it as self._widget.
        """

        pass

    def update(self,*args,**kwargs):
        """
        Redefine this in a subclass. Update the widgets in response to changes
        elsewhere in the interface. 
        """

        pass

    def get_values(self):
        """
        Redefine this in a subclass. Should return widget values as a list or
        dictionary. 
        """
        
        return []

    def _watcher(self,*args,**kwargs):
        """
        Attach this to a particular subwidget (say, a field)  using 'observe'
        and it will pass the values from entire collection of widgets to 
        the update callback.
        
        some_file.observe(self._watcher)
        """
        
        if self._update_callback is not None:
            self._update_callback(self.get_values())

    def _silent_update(self,w,key,value):
        """
        Update specific widget without propagating to observable.
        """

        w.unobserve(self._watcher)
        w.set_trait(key,value)
        w.observe(self._watcher)

    @property
    def widget(self):
        """
        Widget that can be added to the user interface. 
        """
        
        return self._widget
    
class VariableWidgetStack(MetaWidget):
    """
    Super MetaWidget that allows creation of new sub-widgets with a "remove"
    button. Should be subclassed before use. Minimally, the developer will want 
    to add some sort of "add_blah" method and add a "get_values" method. 
    """
    
    def __init__(self,
                 update_callback=None,
                 widget_to_stack=None,
                 button_description="Add",
                 hr_between_rows=False):
        """
        Create new set of widgets, each with a remove button.
        
        Parameters
        ----------
        update_callback : callable, optional
            if defined, call this function without arguments whenever a 
            widget is added or subtracted. This function should take 
            the widget values as its only argument.
        widget_to_stack : MetaWidget, optional
            widget to stack. Should be the class, not an instance. If not
            specified, the add button will not do anything
        button_description : str, optional
            words to put on the add button
        hr_between_rows : bool, default=False
            add <hr/> between elements added
        """

        super().__init__(update_callback=update_callback)

        self._button_to_widget = {}
        self._widget_to_metawidget = {}

        self._widget_to_stack = widget_to_stack
    
        self._add_button = widgets.Button(description=button_description,
                                          disabled=False,
                                          icon='fa-plus')
        self._add_button.on_click(self.add_widget)

        self._hr_between_rows = hr_between_rows

        self._widget = widgets.VBox([self._add_button])

    def add_widget(self,
                   button=None,
                   **kwargs):
        """
        Add a new object to the stack. Can be called as a button click callback
        or directly via the api. NOTE: No error checking is done in this
        function. We assume this is done via the widgets in the gui. 
         
        Parameters
        ----------
        button : None or ipywidgets.Button
            The "button" argument is here so this can be called as a button
            callback, which passes the button instance as the first argument. 
        **kwargs : 
            passed directly to the __init__ function of self._widget_to_stack.
        """

        if self._widget_to_stack is None:
            return
        
        # Pass in update callback for this class if none specified. 
        if "update_callback" not in kwargs:
            kwargs["update_callback"] = self._update_callback

        meta_widget = self._widget_to_stack(**kwargs)
                
        self._add_with_remove_button(some_meta_widget=meta_widget)

    def _append_widget(self,some_widget):
        """
        Append some_widget to the children of self._widget.
        """

        # Actually append the widget to _base_box. 
        current_widgets = list(self._widget.children)
        current_widgets.append(some_widget)
        self._widget.children = tuple(current_widgets)
    
    def _insert_widget(self,index,some_widget):
        """
        Insert some_widget to the children of self._widget at position index.
        """

        # Actually append the widget to _base_box. 
        current_widgets = list(self._widget.children)
        current_widgets.insert(index,some_widget)
        self._widget.children = tuple(current_widgets)
        
    def _pop_widget(self,index):
        """
        Remove a widget (by position) from self._widget.
        """

        current_widgets = list(self._widget.children)
        
        try:
            w = current_widgets[index]
        except IndexError:
            return None

        self._remove_widget(w)
        
    def _remove_widget(self,some_widget):
        """
        Pop a widget (by identity) from _widget.
        """
        
        current_widgets = list(self._widget.children)
        for i, c in enumerate(current_widgets):
            if c is some_widget:
                current_widgets.pop(i)
                break
        self._widget.children = tuple(current_widgets)

        
    def _add_with_remove_button(self,
                                some_meta_widget,
                                widget_container=widgets.HBox,
                                index=-1):
        """
        Add a widget to the ui with an attached remove button.
        """
        
        # Create a button to remove the widget
        remove_button = widgets.Button(description='',
                                       disabled=False,
                                       tooltip='Remove entry',
                                       icon='fa-minus',
                                       layout=widgets.Layout(width='35px'))
        
        # Remove call back will remove the widget_container we
        # create, and thus some_widget and the remove button
        remove_button.on_click(self._remove_button_callback)
        
        # We actually add some_widget within a widget_container 
        # that has both the widget and the remove button
        to_add = widget_container([remove_button,
                                   some_meta_widget.widget])
        if self._hr_between_rows:
            to_add = widgets.VBox([to_add,
                                   widgets.HTML("<hr/>")])
                    
        # Add the widget to the gui
        if index is None:
            self._append_widget(to_add)
        else:
            self._insert_widget(index,to_add)
            
        if self._update_callback is not None:
            self._update_callback(self.get_values())
    
        # Record that we added the widget
        self._button_to_widget[remove_button] = to_add
        self._widget_to_metawidget[to_add] = some_meta_widget
                
    def _remove_button_callback(self,button):
        """
        Callback for a remove button.
        """

        # Figure out the widget to remove
        to_remove = self._button_to_widget[button]
                
        # Remove wiget
        self._remove_widget(to_remove)
        
        # Remove it from the control dictionaries
        widget = self._button_to_widget.pop(button)
        self._widget_to_metawidget.pop(widget)
        
        # Call the removal callback for this object
        if self._update_callback is not None:
            self._update_callback(self.get_values())
        
    def get_values(self):
        """
        Return a list of values from all subwidgets. 
        """
        
        values = []
        for k in self._widget_to_metawidget:
            values.append(self._widget_to_metawidget[k].get_values())
        
        return values
    
    @property
    def widgets(self):
        """
        Return a list of metawidgets that are in the stack. 
        """

        return list(self._widget_to_metawidget.values())
        