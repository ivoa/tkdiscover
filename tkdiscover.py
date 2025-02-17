#!/usr/bin/python

"""
This is a demonstration of pyVO's global dataset discovery
functionality.

It is a Tk GUI that lets you define a region in space, time, and spectrum and
then runs the pyvo query, reporting on its progress.  You can save the current
result into a VOTable at any time.

Distributed under CC-0 (a.k.a public domain).
"""

import datetime
import queue
import io
import re
import os
import pickle
import sys
import tempfile
import threading
import tkinter
from tkinter import messagebox
from tkinter import ttk

from astropy.io import votable
from astropy.io.votable import tree as V
from astropy import time
from astropy import units as u

from pyvo import discover
from pyvo import samp


DEFAULT_SPACE = "274.6880, -13.7920 1"
DEFAULT_SPECTRUM = "2"
DEFAULT_TIME = "1995-01-01 1995-12-31"

FLOAT_RE = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
SEPARATOR_RE = r"(?:\s*,\s*|\s+)"
DATE_RE = r"\d\d\d\d-\d\d-\d\d"

DEBUG_PICKLES = ("last_results.pickle", "last_logs.pickle")

# this a map from VOTable datatypes to ignored standin values
# in case something is None.  Since we're writing binary2 with its
# safe NULLs, this shouldn't matter.  We should really fix astropy
# to translate None-s in their arrays to None automatically.
VOT_NULL_STANDIN = {
    "long": -9999,
    "int": -9999,
    "short": -9999,
    "double": float("NaN"),
    "float": float("NaN"),
    "char": " "}

# The following metadata is generated from DaCHS //obscore RD
# using the extract_obscore_metadata script in the distribution.
# To regenerate it, have DaCHS installed, go to the line with the
# assignment and on an ex commandline say:
OBSCORE_METADATA = \
[{'arraysize': '*',
  'datatype': 'char',
  'description': 'High level scientific classification of the data product, '
                 'taken from an enumeration',
  'name': 'dataproduct_type',
  'ucd': 'meta.code.class',
  'utype': 'obscore:obsdataset.dataproducttype',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Data product specific type',
  'name': 'dataproduct_subtype',
  'ucd': 'meta.code.class',
  'utype': 'obscore:obsdataset.dataproductsubtype',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'short',
  'description': 'Amount of data processing that has been applied to the data',
  'name': 'calib_level',
  'null': -9999,
  'ucd': 'meta.code;obs.calib',
  'utype': 'obscore:obsdataset.caliblevel',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Name of a data collection (e.g., project name) this data '
                 'belongs to',
  'name': 'obs_collection',
  'ucd': 'meta.id',
  'utype': 'obscore:dataid.collection',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Unique identifier for an observation',
  'name': 'obs_id',
  'ucd': 'meta.id',
  'utype': 'obscore:DataID.observationID',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Free-from title of the data set',
  'name': 'obs_title',
  'ucd': 'meta.title;obs',
  'utype': 'obscore:dataid.title',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Dataset identifier assigned by the publisher.',
  'name': 'obs_publisher_did',
  'ucd': 'meta.ref.ivoid',
  'utype': 'obscore:curation.publisherdid',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Dataset identifier assigned by the creator.',
  'name': 'obs_creator_did',
  'ucd': 'meta.id',
  'utype': 'obscore:dataid.creatordid',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'The URL at which to obtain the data set.',
  'name': 'access_url',
  'ucd': 'meta.ref.url',
  'utype': 'obscore:access.reference',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'MIME type of the resource at access_url',
  'name': 'access_format',
  'ucd': 'meta.code.mime',
  'utype': 'obscore:access.format',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'long',
  'description': 'Estimated size of data product',
  'name': 'access_estsize',
  'null': -9999,
  'ucd': 'phys.size;meta.file',
  'utype': 'obscore:access.size',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Object a targeted observation targeted',
  'name': 'target_name',
  'ucd': 'meta.id;src',
  'utype': 'obscore:Target.Name',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Class of the target object (star, QSO, ...)',
  'name': 'target_class',
  'ucd': 'src.class',
  'utype': 'obscore:target.class',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'RA of (center of) observation, ICRS',
  'name': 's_ra',
  'ucd': 'pos.eq.ra',
  'utype': 'obscore:char.spatialaxis.coverage.location.coord.position2d.value2.c1',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Dec of (center of) observation, ICRS',
  'name': 's_dec',
  'ucd': 'pos.eq.dec',
  'utype': 'obscore:char.spatialaxis.coverage.location.coord.position2d.value2.c2',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Approximate spatial extent for the region covered by the '
                 'observation',
  'name': 's_fov',
  'ucd': 'phys.angSize;instr.fov',
  'utype': 'obscore:char.spatialaxis.coverage.bounds.extent.diameter',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Region covered by the observation, as a polygon',
  'name': 's_region',
  'ucd': 'pos.outline;obs.field',
  'utype': 'obscore:char.spatialaxis.coverage.support.area',
  'xtype': 'adql:REGION'},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Best spatial resolution within the data set',
  'name': 's_resolution',
  'ucd': 'pos.angResolution',
  'utype': 'obscore:Char.SpatialAxis.Resolution.refval.value',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Lower bound of times represented in the data set',
  'name': 't_min',
  'ucd': 'time.start;obs.exposure',
  'utype': 'obscore:char.timeaxis.coverage.bounds.limits.starttime',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Upper bound of times represented in the data set',
  'name': 't_max',
  'ucd': 'time.end;obs.exposure',
  'utype': 'obscore:char.timeaxis.coverage.bounds.limits.stoptime',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'float',
  'description': 'Total exposure time',
  'name': 't_exptime',
  'ucd': 'time.duration;obs.exposure',
  'utype': 'obscore:char.timeaxis.coverage.support.extent',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'float',
  'description': 'Minimal significant time interval along the time axis',
  'name': 't_resolution',
  'ucd': 'time.resolution',
  'utype': 'obscore:char.timeaxis.resolution.refval.value',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Minimal wavelength represented within the data set',
  'name': 'em_min',
  'ucd': 'em.wl;stat.min',
  'utype': 'obscore:char.spectralaxis.coverage.bounds.limits.lolimit',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Maximal wavelength represented within the data set',
  'name': 'em_max',
  'ucd': 'em.wl;stat.max',
  'utype': 'obscore:char.spectralaxis.coverage.bounds.limits.hilimit',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Spectral resolving power lambda/delta lambda',
  'name': 'em_res_power',
  'ucd': 'spect.resolution',
  'utype': 'obscore:char.spectralaxis.resolution.resolpower.refval',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': "UCD for the product's observable",
  'name': 'o_ucd',
  'ucd': 'meta.ucd',
  'utype': 'obscore:char.observableaxis.ucd',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'List of polarization states in the data set',
  'name': 'pol_states',
  'ucd': 'meta.code;phys.polarization',
  'utype': 'obscore:Char.PolarizationAxis.stateList',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Name of the facility at which data was taken',
  'name': 'facility_name',
  'ucd': 'meta.id;instr.tel',
  'utype': 'obscore:Provenance.ObsConfig.facility.name',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'Name of the instrument that produced the data',
  'name': 'instrument_name',
  'ucd': 'meta.id;instr',
  'utype': 'obscore:Provenance.ObsConfig.instrument.name',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'long',
  'description': 'Number of elements (typically pixels) along the first '
                 'spatial axis.',
  'name': 's_xel1',
  'null': -9999,
  'ucd': 'meta.number',
  'utype': 'obscore:Char.SpatialAxis.numBins1',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'long',
  'description': 'Number of elements (typically pixels) along the second '
                 'spatial axis.',
  'name': 's_xel2',
  'null': -9999,
  'ucd': 'meta.number',
  'utype': 'obscore:Char.SpatialAxis.numBins2',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'long',
  'description': 'Number of elements (typically pixels) along the time axis.',
  'name': 't_xel',
  'null': -9999,
  'ucd': 'meta.number',
  'utype': 'obscore:Char.TimeAxis.numBins',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'long',
  'description': 'Number of elements (typically pixels) along the spectral '
                 'axis.',
  'name': 'em_xel',
  'null': -9999,
  'ucd': 'meta.number',
  'utype': 'obscore:Char.SpectralAxis.numBins',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'long',
  'description': 'Number of elements (typically pixels) along the polarization '
                 'axis.',
  'name': 'pol_xel',
  'null': -9999,
  'ucd': 'meta.number',
  'utype': 'obscore:Char.PolarizationAxis.numBins',
  'xtype': None},
 {'arraysize': None,
  'datatype': 'double',
  'description': 'Sampling period in world coordinate units along the spatial '
                 'axis',
  'name': 's_pixel_scale',
  'ucd': 'phys.angSize;instr.pixel',
  'utype': 'obscore:Char.SpatialAxis.Sampling.RefVal.SamplingPeriod',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': "Nature of the product's spectral axis (typically, em.freq, "
                 'em.wl, or em.energy)',
  'name': 'em_ucd',
  'ucd': 'meta.ucd',
  'utype': 'obscore:Char.SpectralAxis.ucd',
  'xtype': None},
 {'arraysize': '*',
  'datatype': 'char',
  'description': 'IVOID of the originating service',
  'name': 'origin_service',
  'ucd': 'meta.ref.ivoid',
  'utype': None,
  'xtype': None}]


class InputSyntaxError(Exception):
    """raised when we cannot parse a constraint.
    """


class TkDiscoverer(tkinter.Tk):
    """The main UI.
    """
    def __init__(self, **kwargs):
        self.debugging = kwargs.pop("debug", False)
        super().__init__(**kwargs)
        self.minsize(900, -1)

        self.space_tk = tkinter.StringVar(self, DEFAULT_SPACE)
        self.spectrum_tk = tkinter.StringVar(self, DEFAULT_SPECTRUM)
        self.time_tk = tkinter.StringVar(self, DEFAULT_TIME)
        self.destname_tk = tkinter.StringVar(self, "discovered.vot")
        self.inclusive_tk = tkinter.IntVar(self)
        self.style = ttk.Style()
        self.disco = None
        self.thread = None

        # several widgets get dis- and enabled depending on
        # whether a thread is running.  They are added
        # here in _make_widgets and then furnished with
        # attributes thread_[running|stopped]_state, which
        # are in turn evaluated in the _update methods below
        self.thread_sensitive_widgets = []
        self._make_widgets()

    def _make_thread_sensitive(self, widget, enabled_when_running):
        # See the comment on thread_sensitive_widgets
        if enabled_when_running:
            widget.thread_running_state = tkinter.NORMAL
            widget.thread_stopped_state = tkinter.DISABLED
        else:
            widget.thread_running_state = tkinter.DISABLED
            widget.thread_stopped_state = tkinter.NORMAL
        self.thread_sensitive_widgets.append(widget)
        return widget

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
        buttonbar.grid(column=1, row=10, columnspan=2, sticky=tkinter.SE)

        if self.debugging:
            ttk.Button(buttonbar,
                text="Debugger", command=self._enter_pdb).pack(
                    side=tkinter.LEFT)

        self._make_thread_sensitive(
            ttk.Button(buttonbar, text="Run", command=self._run),
            False).pack(side=tkinter.LEFT)
        self._make_thread_sensitive(
            ttk.Button(buttonbar, text="Stop", command=self._stop),
            True).pack(side=tkinter.LEFT)
        ttk.Button(buttonbar,
            text="Broadcast", command=self._broadcast).pack(side=tkinter.LEFT)
        ttk.Button(buttonbar,
            text="Save", command=self._save).pack(side=tkinter.LEFT)
        self._make_thread_sensitive(
            ttk.Button(buttonbar, text="Quit", command=self.quit),
            False).pack(side=tkinter.LEFT)

        # Statusbar

        status_frame = ttk.Frame(self)
        status_frame.grid(column=1, row=11, columnspan=2, sticky=stickall)
        self.style.map("Label.statusbar", relief=tkinter.SUNKEN)

        self.status_label = ttk.Label(status_frame, relief=tkinter.SUNKEN)
        self.counter_label = ttk.Label(status_frame, text="-/-")
        self.status_label.grid(row=0, column=1, sticky=stickall)
        self.counter_label.grid(row=0, column=2, sticky=stickall)
        status_frame.columnconfigure(1, weight=1)

        self.rowconfigure(10, weight=1)
        self.columnconfigure(2, weight=1)

        self._update_for_stopped_thread()

    def _update_status(self, log_message):
        self.status_label.configure(text=log_message)
        self._update_service_counter()

    def _update_service_counter(self):
        self.counter_label.configure(
            text="{}/{}/{}".format(*self.disco.get_query_stats())
            +": {}".format(len(self.disco.results)))

    ################## Parsing and UI

    def _enter_pdb(self):
        """a callback for the debug button.

        Note that pdb will run in the tkinter event loop, so the whole
        thing will freeze until someone leaves the debugger.  Only
        enable this when there's a terminal on which people can do this.
        """
        import pdb;pdb.Pdb(nosigint=True).set_trace()

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

        return ("time", (time_min, time_min))

    def _make_constraints(self):
        return [c for c in (
            self._parse_space_constraint(),
            self._parse_spectral_constraint(),
            self._parse_temporal_constraint())
            if c]

    ################ the thread watcher and related stuff
    # (the thread watcher is a tkinter timer that pulls messages
    # from the thread's message queue and additionally notices
    # when the thread has finished and cleans up after it
    def _watch_thread(self):
        try:
            msg = None
            while True:
                try:
                    msg = self._message_queue.get_nowait()
                except queue.Empty:
                    break
            if msg is not None:
                self._update_status(msg)

            if not self.thread.is_alive():
                self.thread.join()
                self.thread = None
                self._update_for_stopped_thread()
        finally:
            if self.thread:
                self.after(100, self._watch_thread)

    def _update_for_running_thread(self):
        """configure the UI for a running discovery.

        That's disabling some buttons, creating the communication
        queue, and starting the thread watcher.
        """
        self._message_queue = queue.Queue()
        for widget in self.thread_sensitive_widgets:
            widget.configure(state=widget.thread_running_state)

        self.after(100, self._watch_thread)

    def _update_for_stopped_thread(self):
        """configure the UI for a stopped discovery.

        That, in particular, stops the thread watcher, too.
        """
        if hasattr(self, "_message_queue"):
            del self._message_queue
        for widget in self.thread_sensitive_widgets:
            widget.configure(state=widget.thread_stopped_state)

        # we are called at initialisation; don't make any claims then.
        if self.disco is not None:
            self._update_status(
                "Discovery finished. {} records found.".format(
                    len(self.disco.results)))

    def _stop(self):
        self.disco.reset_services()
        if self.thread:
            self._update_status(
                "Please wait for the last service to complete")

    def _queue_message(self, msg):
        """queues a status bar message.

        This is so the subordinate thread can communicate with the
        front end.
        """
        self._message_queue.put(msg)

    ################# discovery and results management

    def _start_thread(self):
        if self.thread is not None:
            raise Exception("We are already processing")

        def thread_target():
            self.disco.discover_services()
            self.disco.query_services()

        self._update_for_running_thread()
        self.thread = threading.Thread(target=thread_target)
        self.thread.start()

    def _run(self):
        """UI callback: start the discovery.
        """
        try:
            constraints = self._make_constraints()
        except InputSyntaxError as exc:
            messagebox.showerror("Tkdiscover Error", str(exc))
            return

        self.disco = discover.ImageDiscoverer(
            **dict(constraints),
            inclusive=self.inclusive_tk.get(),
            watcher=lambda _, msg: self._queue_message(msg))
        if self.debugging and os.path.exists(DEBUG_PICKLES[0]):
            with open(DEBUG_PICKLES[0], "rb") as f:
                self.disco.results = pickle.load(f)
            with open(DEBUG_PICKLES[1], "rb") as f:
                self.disco.log_messages = pickle.load(f)
            return

        self._update_status("Obtaining services...")
        self.update_idletasks()
        self._start_thread()

    def _get_current_votable(self):
        """returns bytes for a VOTable of the current results.

        This returns does some user interaction and returns None when
        there are no results yet.  That's lazy and probably should
        be changed.
        """
        if self.disco is None:
            messagebox.showerror(
                "Tkdiscover Error",
                "No results so far (try Run?)")
            return

        # Copy these things more or less simultaneously because we
        # may still be querying.  We *should* use a lock here, really.
        results = self.disco.results[:]
        log_messages = self.disco.log_messages[:]
        if self.debugging:
            with open(DEBUG_PICKLES[0], "wb") as f:
                pickle.dump(results, f)
            with open(DEBUG_PICKLES[1], "wb") as f:
                pickle.dump(log_messages, f)

        vot = V.VOTableFile(config={'version_1_2_or_later': True})
        resource = V.Resource()
        vot.resources.append(resource)
        resource.infos.append(
            V.Info(
                name="provenance",
                value="tkdiscover.py global query at {}".format(
                    datetime.datetime.utcnow())))
        # TODO: serialise query parameters
        resource.infos.append(
            V.Info(
                name="discovery_log",
                value="\n".join(log_messages)))

        table = V.Table(vot)
        resource.tables.append(table)

        for meta in OBSCORE_METADATA:
            meta = meta.copy()  # so we can pop the description
            description = meta.pop("description")
            null = meta.pop("null", None)

            field = V.Field(vot, **meta)
            field.description = description
            if null is not None:
                field.values = V.Values(votable, field, null=null)

            table.fields.append(field)

        def make_record(imageFound):
            values, mask = [], []
            for meta in OBSCORE_METADATA:
                val = getattr(imageFound, meta["name"], None)
                if val is None:
                    values.append(VOT_NULL_STANDIN[meta["datatype"]])
                    mask.append(True)
                else:
                    values.append(val)
                    mask.append(False)

            return tuple(values), tuple(mask)

        table.create_arrays(len(results))
        for index, row in enumerate(results):
            try:
                serialised, mask = make_record(row)
                table.array[index] = serialised
                table.array.mask[index] = mask
            except Exception as exc:
                import pdb;pdb.set_trace()
                sys.stderr.write(
                    f"Cannot serialise a row ({exc}): {serialised}\n")

        dest = io.BytesIO()
        votable.writeto(vot, dest, "binary2")
        return dest.getvalue()

    def _save(self):
        """UI callback: save the current results.
        """
        vot = self._get_current_votable()
        if vot:
            with open(self.destname_tk.get(), "wb") as f:
                f.write(vot)

    def _send_to_samp_client(self, conn, client_id, message):
        """sends a SAMP message to client_id via the SAMP connection
        conn.

        This is purely a helper for _broadcast_samp.
        """
        try:
            conn.call_and_wait(client_id, message, "10")
        except Exception:
            # we don't want to crash on a misbehaved thing
            # on the SAMP bus; I don't think the user
            # should even learn of any trouble.  But we don't want
            # to entirely silent, either, so, for now:
            import traceback; traceback.print_exc()

    def _broadcast_samp(self, serialised_data):
        """does a samp broadcast of serialised_data on a connection
        of its own.

        In the context of tkdiscover, this is running in a dedicated
        thread that terminates when the transmission is done.
        """
        # TODO: Errors in here are basically ignored (well, you'll
        # see tracebacks on the console).  To change this, we'd probably
        # have to re-use the threadwatcher queue; but that would then
        # have to be examined all the time.
        handle, f_name = tempfile.mkstemp(suffix=".xml")
        with open(handle, "wb") as f:
            f.write(serialised_data)

        try:
            with samp.connection(
                    client_name="tkdiscover",
                    description="PyVO global dataset discovery") as conn:
                message = {
                    "samp.mtype": "table.load.votable",
                    "samp.params": {
                        "url": "file://"+f_name,
                        "name": "tkdiscover current result",
                    },
                }
                for client_id in conn.get_subscribed_clients(
                        "table.load.votable"):
                    self._send_to_samp_client(
                        conn, client_id, message)
        finally:
            os.unlink(f_name)

    def _broadcast(self):
        """UI callback: SAMP-broadcast the current results
        """
        vot = self._get_current_votable()
        if vot:
            samp_thread = threading.Thread(
                target=lambda: self._broadcast_samp(vot))
            samp_thread.start()

            def harvest_thread():
                if samp_thread.is_alive():
                    self.after(1000, harvest_thread)
                else:
                    samp_thread.join()

            self.after(1000, harvest_thread)


def parse_command_line():
    import argparse
    parser = argparse.ArgumentParser(
        description="Global dataset discovery for the Virtual Observatory")
    parser.add_argument("--debug",
        action="store_true",
        help="Add a debug button to drop into pdb, and cache query results"
            " in pickle files.  Don't use unless you understand the"
            " implications.")

    return parser.parse_args()


def main():
    args = parse_command_line()
    ui = TkDiscoverer(debug=args.debug)
    ui.title("Global VO Image Search")
    ui.mainloop()


if __name__=="__main__":
    main()
