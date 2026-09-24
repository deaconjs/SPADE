"""Tk widget that hosts a VTK render window.

create_render_widget() uses VTK's own vtkTkRenderWindowInteractor when the
compiled vtkRenderingTk module is available (e.g. some conda-forge builds).
The pip wheels don't ship it, so otherwise it falls back to TkVTKWidget, which
embeds the render window in a Tk frame by its native window handle and
forwards Tk mouse/keyboard events to a vtkGenericRenderWindowInteractor.
The fallback works on Windows and X11, but not on macOS.

Both widgets are ordinary Tk widgets (pack/grid/forget) and provide
GetRenderWindow(), GetInteractor() and Render().
"""
import tkinter
import vtk


class TkVTKWidget(tkinter.Frame):
    _press_events = {1: 'LeftButtonPressEvent', 2: 'MiddleButtonPressEvent', 3: 'RightButtonPressEvent'}
    _release_events = {1: 'LeftButtonReleaseEvent', 2: 'MiddleButtonReleaseEvent', 3: 'RightButtonReleaseEvent'}

    def __init__(self, master, width=400, height=400, rw=None, **kw):
        # an empty background keeps Tk from painting over the OpenGL surface
        kw.setdefault('background', '')
        tkinter.Frame.__init__(self, master, width=width, height=height, takefocus=1, **kw)
        self._render_window = rw if rw else vtk.vtkRenderWindow()
        self._interactor = vtk.vtkGenericRenderWindowInteractor()
        self._interactor.SetRenderWindow(self._render_window)
        self._interactor.AddObserver('CreateTimerEvent', self._create_timer)
        self._interactor.AddObserver('DestroyTimerEvent', lambda obj, event: 1)
        # winfo_id() forces Tk to create the native window, so VTK can render
        # into it even before the frame is mapped
        self._render_window.SetWindowInfo(str(self.winfo_id()))
        self._render_window.SetSize(width, height)
        self._interactor.Initialize()
        self._bind_events()

    def _bind_events(self):
        self.bind('<Configure>', self._on_configure)
        self.bind('<Expose>', lambda event: self.Render())
        self.bind('<Destroy>', self._on_destroy)
        for button in (1, 2, 3):
            self.bind('<ButtonPress-%d>' % button, lambda event, b=button: self._on_press(event, b))
            self.bind('<ButtonRelease-%d>' % button, lambda event, b=button: self._on_release(event, b))
        self.bind('<Motion>', self._on_motion)
        self.bind('<MouseWheel>', self._on_wheel)                                   # Windows
        self.bind('<Button-4>', lambda event: self._on_wheel(event, forward=True))  # X11
        self.bind('<Button-5>', lambda event: self._on_wheel(event, forward=False))
        self.bind('<Enter>', self._on_enter)
        self.bind('<Leave>', self._on_leave)
        self.bind('<KeyPress>', self._on_key_press)
        self.bind('<KeyRelease>', self._on_key_release)

    def GetRenderWindow(self):
        return self._render_window

    def GetInteractor(self):
        return self._interactor

    def Render(self):
        self._render_window.Render()

    def _create_timer(self, obj, event):
        self.after(10, self._interactor.TimerEvent)

    def _set_event_info(self, event, keycode=0, keysym=None):
        ctrl = int((event.state & 0x4) != 0)
        shift = int((event.state & 0x1) != 0)
        self._interactor.SetEventInformationFlipY(event.x, event.y, ctrl, shift, chr(keycode), 0, keysym)

    def _on_configure(self, event):
        self._interactor.UpdateSize(event.width, event.height)
        self._interactor.ConfigureEvent()
        self.Render()

    def _on_destroy(self, event):
        if event.widget is self:
            self._render_window.Finalize()

    def _on_press(self, event, button):
        self.focus_set()
        self._set_event_info(event)
        self._interactor.InvokeEvent(self._press_events[button])

    def _on_release(self, event, button):
        self._set_event_info(event)
        self._interactor.InvokeEvent(self._release_events[button])

    def _on_motion(self, event):
        self._set_event_info(event)
        self._interactor.MouseMoveEvent()

    def _on_wheel(self, event, forward=None):
        if forward is None:
            forward = event.delta > 0
        self._set_event_info(event)
        if forward:
            self._interactor.MouseWheelForwardEvent()
        else:
            self._interactor.MouseWheelBackwardEvent()

    def _on_enter(self, event):
        self.focus_set()
        self._set_event_info(event)
        self._interactor.EnterEvent()

    def _on_leave(self, event):
        self._set_event_info(event)
        self._interactor.LeaveEvent()

    def _on_key_press(self, event):
        keycode = ord(event.char) if len(event.char) == 1 else 0
        self._set_event_info(event, keycode, event.keysym)
        self._interactor.KeyPressEvent()
        self._interactor.CharEvent()

    def _on_key_release(self, event):
        keycode = ord(event.char) if len(event.char) == 1 else 0
        self._set_event_info(event, keycode, event.keysym)
        self._interactor.KeyReleaseEvent()


def native_tk_available(master):
    """True if this VTK build includes the compiled vtkRenderingTk module."""
    try:
        from vtkmodules.tk.vtkLoadPythonTkWidgets import vtkLoadPythonTkWidgets
        vtkLoadPythonTkWidgets(master.tk)
        return True
    except (ImportError, tkinter.TclError):
        return False


def create_render_widget(master, width=400, height=400):
    if native_tk_available(master):
        from vtkmodules.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
        widget = vtkTkRenderWindowInteractor(master, width=width, height=height)
        widget.Initialize()
        widget.Start()
    else:
        widget = TkVTKWidget(master, width=width, height=height)
    # drag-to-rotate, like the old vtkTkRenderWidget; 'j'/'t' keys still switch styles
    widget.GetRenderWindow().GetInteractor().GetInteractorStyle().SetCurrentStyleToTrackballCamera()
    return widget
