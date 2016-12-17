import sys
sys.path.append('Dependencies')
import Pmw
import string

class functional_pmw_counter(Pmw.Counter):
    def __init__(self, a, b, c, d, e, f, g):
        Pmw.Counter.__init__(self, a, labelpos = b, label_text = c, entryfield_value = d, entryfield_validate = e, datatype = f, increment = g)
        self.parent_frame        = a
        self.labelpos            = b
        self.label_text          = c
        self.entryfield_value    = d
        self.entryfield_validate = e
        self.datatype            = f
        self.increment           = g
        self.function = self.empty_function
        self.enabled = 1
        self.active = 0
    def disable(self):
        self.component('entryfield').component('entry').config(state='disabled')
        self.enabled = 0
    def enable(self):
        self.component('entryfield').component('entry').config(state='normal')
        self.enabled = 1
    def empty_function(self):
        pass
    def update_function(self, new_function):
        self.function = new_function
    def _countUp(self, event):
        # override for disable control and attaching the increment/decrement to the filename
        self.active=1
        if self.enabled:
            Pmw.Counter._countUp(self, event)
    def _countDown(self, event):
        self.active=1
        if self.enabled:
            Pmw.Counter._countDown(self, event)
    def _stopCounting(self, event):
        Pmw.Counter._stopCounting(self, event)
        if self.active==1:
            t = self.component('entryfield').getvalue()
            if self.entryfield_validate['validator'] == 'time':
                t_tokens = string.split(t, ':')
                a = string.atoi(t_tokens[0])*60*60
                b = string.atoi(t_tokens[1])*60
                c = string.atoi(t_tokens[2])
                seconds = float(a + b + c)
                self.function(seconds)
            elif self.entryfield_validate['validator'] == 'integer':
                self.function()
            self.active=0
