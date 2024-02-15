#!/usr/bin/python

"""
This is a demonstration of pyVO's global dataset discovery
functionality.

It is a Tk GUI that lets you define a region in space, time, and spectrum and
then runs the pyvo query, reporting on its progress.  You can save the current
result into a VOTable at any time.

Distributed under CC-0 (a.k.a public domain).
"""

import re
import tkinter
from tkinter import messagebox
from tkinter import ttk

from astropy import units as u
from astropy import time

from pyvo import discover


DEFAULT_SPACE = "274.6880, -13.7920 0.1"
DEFAULT_SPECTRUM = "2"
DEFAULT_TIME = "1995-01-01 1995-12-31"

FLOAT_RE = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
SEPARATOR_RE = r"(?:\s*,\s*|\s+)"
DATE_RE = r"\d\d\d\d-\d\d-\d\d"


class InputSyntaxError(Exception):
    """raised when we cannot parse a constraint.
    """


class TkDiscoverer(tkinter.Tk):
    """The main UI.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.space_tk = tkinter.StringVar(self, DEFAULT_SPACE)
        self.spectrum_tk = tkinter.StringVar(self, DEFAULT_SPECTRUM)
        self.time_tk = tkinter.StringVar(self, DEFAULT_TIME)
        self.destname_tk = tkinter.StringVar(self, "discovered.vot")
        self.inclusive_tk = tkinter.IntVar(self)
        self.style = ttk.Style()
        self.discoverer = None
        self._make_widgets()

    def _make_widgets(self):
        stickall = tkinter.NE+tkinter.SW

        # Constraint input

        ttk.Label(self, text="Space (RA, Dec, [SR]) [deg]"
            ).grid(column=1, row=1, sticky=tkinter.W)
        space_entry = ttk.Entry(self, textvariable=self.space_tk)
        space_entry.grid(column=2, row=1, sticky=stickall)

        ttk.Label(self, text="Spectrum (λ [nm])"
            ).grid(column=1, row=2, sticky=tkinter.W)
        spect_entry = ttk.Entry(self, textvariable=self.spectrum_tk)
        spect_entry.grid(column=2, row=2, sticky=stickall)

        ttk.Label(self, text="Time (one or two YYYY-MM-DD)"
            ).grid(column=1, row=3, sticky=tkinter.W)
        date_entry = ttk.Entry(self, textvariable=self.time_tk)
        date_entry.grid(column=2, row=3, sticky=stickall)

        # Results configuration

        ttk.Checkbutton(
            self,
            text="Search inclusive",
            variable=self.inclusive_tk
            ).grid(column=1, columnspan=2, row=4, sticky=tkinter.E)

        ttk.Label(self, text="Save results to"
            ).grid(column=1, row=9, sticky=tkinter.W)
        filename_entry = ttk.Entry(self, textvariable=self.destname_tk)
        filename_entry.grid(column=2, row=9, sticky=stickall)

        # Button bar

        buttonbar = ttk.Frame(self)
        buttonbar.grid(column=1, row=10, columnspan=2, sticky=tkinter.S)

        ttk.Button(buttonbar,
            text="Run", command=self._run).pack(side=tkinter.LEFT)
        ttk.Button(buttonbar,
            text="Save", command=self._save).pack(side=tkinter.LEFT)
        ttk.Button(buttonbar,
            text="Quit", command=self.quit).pack(side=tkinter.LEFT)

        # Statusbar

        status_frame = ttk.Frame(self)
        status_frame.grid(column=1, row=11, columnspan=2, sticky=stickall)
        self.style.map("Label.statusbar", relief=tkinter.SUNKEN)

        self.status_label = ttk.Label(status_frame, relief=tkinter.SUNKEN)
        self.counter_label = ttk.Label(status_frame, text="0/0")
        self.status_label.grid(row=0, column=1, sticky=stickall)
        self.counter_label.grid(row=0, column=2, sticky=stickall)
        status_frame.columnconfigure(1, weight=1)

        self.rowconfigure(10, weight=1)
        self.columnconfigure(2, weight=1)

    def _update_log(self, log_message):
        self.status_label.configure(text=log_message)

    ################## Parsing and UI

    def _parse_space_constraint(self):
        literal = self.space_tk.get().strip()
        if not literal:
            return None

        mat = re.match("({}){}({}){}({})$".format(
            FLOAT_RE, SEPARATOR_RE, FLOAT_RE, SEPARATOR_RE, FLOAT_RE),
            literal)
        if not mat:
            raise InputSyntaxError(f"Cannot parse as a cone:\n`{literal}'"
                +"\nWrite RA, Dec, Radius, all in decimal degrees.")

        return ("space", tuple(float(s) for s in mat.groups()))

    def _parse_spectral_constraint(self):
        literal = self.spectrum_tk.get().strip()
        if not literal:
            return None

        mat = re.match(f"{FLOAT_RE}", literal)
        if not mat:
            raise InputSyntaxError(
                f"Cannot parse as a spectral coordinate:\n`{literal}'"
                +"\nThis is just a float in nanometers now")

        return ("spectrum", float(mat.group())*u.nm)

    def _parse_temporal_constraint(self):
        literal = self.time_tk.get().strip()
        if not literal:
            return None

        mat = re.match(f"({DATE_RE}){SEPARATOR_RE}({DATE_RE})?", literal)

        if not mat:
            raise InputSyntaxError(
                f"Cannot parse as a time constraint:\n`{literal}'"
                +"\nUse one or two YYYY-MM-DD dates.")

        time_min = time.Time(mat.group(1))
        if mat.group(2) is not None:
            return ("time", (time_min, time.Time(mat.group(2))))
        else:
            return ("time", (time_min, time_max))

    def _make_constraints(self):
        return [c for c in (
            self._parse_space_constraint(),
            self._parse_spectral_constraint(),
            self._parse_temporal_constraint())
            if c]

    def _run(self):
        """UI callback: start the discovery.
        """
        try:
            constraints = self._make_constraints()
        except InputSyntaxError as exc:
            messagebox.showerror("Bad Constraint Input", str(exc))
            return

        self.disco = discover.ImageDiscoverer(
            **dict(constraints),
            watcher=self._update_log)

        self._update_log("Obtaining services...")
        self.update_idletasks()

        self.disco.discover_services()

        self.counter_label.configure(
            text="{}/{}".format(*self.disco.get_query_stats()))

    def _save(self):
        """UI callback: save the current results.
        """

def main():
    ui = TkDiscoverer()
    ui.title("Global VO Image Search")
    ui.mainloop()


if __name__=="__main__":
    main()
