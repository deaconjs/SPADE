from Tkinter import *
from bwidget import *
t = Tk()
b = Button(t, relief=LINK, text="bwidget button", command=t.destroy)
b.pack()
t.mainloop()
