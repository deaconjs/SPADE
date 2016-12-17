import bwidget, Tkinter, sys, os

app = Tkinter.Tk(); app.wm_title("pybwidget demo")
notebook = bwidget.NoteBook(app, arcradius=2); notebook.pack()

files = sys.argv[0], Tkinter.__file__, bwidget.__file__
for i, f in enumerate(files):
    if f.endswith(".pyc"): f = f[:-1]
    page = notebook.insert(Tkinter.END, i, text=os.path.basename(f))
    sw = bwidget.ScrolledWindow(page)
    text = Tkinter.Text(sw)
    text.insert(Tkinter.END, open(f).read())
    sw.setwidget(text)
    sw.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=True)
notebook.compute_size()

app.mainloop()
