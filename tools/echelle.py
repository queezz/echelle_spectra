"""
Tools for Echelle spectrometer images
"""
from sif_reader import np_open
import numpy as np
import os
from os.path import join

DIMW = 1024  # pixel dimension in the wavelength direction
DIMO = 1024  # pixel dimension in the order direction

remove_npnans = lambda a: a[~np.isnan(a)]  # remove nans from numpy array


class EchelleImage:
    """"""

    def __init__(self, fpth, clbr=None):
        """ Sif Image container with tools to convert Echelle image into spectra
        fpth - path to sif file
        """
        self.fpth = fpth
        self.read_image()
        if clbr is not None:
            self.clbr = clbr

    def read_image(self):
        """
        """
        images, self.info = np_open(self.fpth)
        images = images.repeat(self.info["ybin"], axis=1)
        self.images = images.repeat(self.info["xbin"], axis=2)

    def plot_frame(self, frame, **kws):
        """show a frame of the sif image from fpth"""
        import matplotlib.pylab as plt
        import matplotlib as mpl

        # font sizes
        fs = {"cb": 6, "cbname": 8, "pattern": 5, "hint": 3}
        fs = {"cb": 12, "cbname": 16, "pattern": 6, "hint": 3}

        scale = kws.get("scale", 1)

        image = self.images[frame]
        ymax, xmax = image.shape

        names = ["cmap", "under", "over", "pcolor", "ord_color"]
        if kws.get("dark", True):
            colors = ["viridis", "k", "w", "k", "w"]
        else:
            colors = ["gist_yarg", "r", "gold", "royalblue", "k"]
        clrs = {n: c for n, c in zip(names, colors)}

        fig = plt.figure()
        ax = fig.add_subplot(111)
        imin = np.min(image)
        if image.min() > imin:
            impin = image.min()
        im = ax.imshow(
            image, cmap=clrs["cmap"], norm=mpl.colors.LogNorm(imin, image.max() * scale)
        )
        # im = plt.imshow(a, norm = mpl.colors.LogNorm(1e2,5e3))
        im.cmap.set_under(clrs["under"])
        im.cmap.set_over(clrs["over"])

        if kws.get("axis", True):
            # ax.text(.65,-.16,'%s frame %s'%(os.path.basename(fpth), frame),
            #             transform = ax.transAxes, fontsize = fs['shotname'])
            cbaxes = fig.add_axes([0.255, 0.99, 0.51, 0.05])
            cb = plt.colorbar(im, extend="both", cax=cbaxes, orientation="horizontal")
            cb.ax.set_xticklabels(
                cb.ax.get_xticklabels(), rotation=0, fontsize=fs["cb"]
            )
            cb.ax.text(
                0.4,
                0.2,
                "counts",
                fontsize=fs["cbname"],
                rotation=0,
                color="w",
                transform=cb.ax.transAxes,
            )
            if kws.get("axlabel", True):
                ax.text(-0.2 * ymax, ymax * 0.95, "vertical pixel", rotation=90)
                ax.text(xmax * 0.65, -xmax * 0.15, "horizontal pixel")
        else:
            plt.axis("off")
            ax.set_xticks([])
            ax.set_yticks([])

        if kws.get("pattern", False):
            # Need calibration to show the pattern
            clbr = kws.get("clbr", None)
            if clbr is None:
                clbr = self.clbr
            clbr.show_masks(dv=8)
            plt.imshow(
                np.invert(clbr.masks.sum(axis=0)),
                origin="lower",
                alpha=0.4,
                cmap="binary",
            )

        ax.set_xlim(0, xmax)
        ax.set_ylim(ymax, 0)
        ax.invert_yaxis()
        plt.sca(ax)
        return

    def order_image(self, frame, ordind, sm=False):
        """ Cut the image along the order
        clbr - instance of Calibrations, already "started"
        """
        clbr = self.clbr
        img = self.images[frame][clbr.cutting_masks[ordind]]
        try:
            ordr = img.reshape(clbr.dv * 2 + 1, clbr.DIMW)
        except ValueError:
            ordr = img.reshape(clbr.dv * 2 + 1, int(img.shape[0] / (clbr.dv * 2 + 1)))
        if sm:
            return ordr.sum(axis=0)
        else:
            return ordr

    def calculate_order_spectra(self):
        """ """
        clbr = self.clbr
        frames = range(self.info["NumberOfFrames"])
        orders = range(clbr.pattern.shape[1])

        self.order_spectra = np.array(
            [
                np.array([self.order_image(fi, o, sm=True) for o in orders])
                for fi in frames
            ]
        )

    def correct_order_shapes(self):
        """
        Correction of the short orders (28th for CMOS camera)
        fill spectra with np.nan to make a np.array with all orders spectra
        """
        clbr = self.clbr
        # Do nothing if there are no np.nan's in orders spectra
        if not len(clbr.orders_bad_shape):
            return

        frames = range(self.info["NumberOfFrames"])

        for o, i in zip(clbr.orders_bad_shape, clbr.orders_bad_froms):
            # print(o,i)
            lnan = clbr.DIMW - self.order_spectra[0][o].shape[0]
            nans = np.ones(lnan) * np.nan
            for fi in frames:
                # self.order_spectra[fi][o][np.isnan(clbr.order_wavel[o])] = np.nan
                self.order_spectra[fi][o] = np.append(self.order_spectra[fi][o], nans)

        self.order_spectra = np.hstack(np.hstack(self.order_spectra)).reshape(
            self.info["NumberOfFrames"], clbr.pattern.shape[1], clbr.DIMW
        )

    def calculate_spectra(self):
        """ calculate continuous spectra from order_spectra
        using order borders from calibration clbr
        """
        clbr = self.clbr

        self.spectra = np.array([i[clbr.order_borders] for i in self.order_spectra])

        a = remove_npnans(self.spectra)

        self.spectra = a.reshape(
            self.info["NumberOfFrames"], int(len(a) / self.info["NumberOfFrames"])
        )

        self.wavelength = remove_npnans(clbr.wavelength)

    def plot_order_image(self, frame, ordind, aspect=2):
        """ Plot order image"""
        import matplotlib.pylab as plt

        plt.imshow(
            self.order_image(frame, ordind), origin="lower", aspect=aspect,
        )

    def plot_cut_image(self, frame, aspect=2):
        """ Plot images of all orders in one picture
        """
        import matplotlib.pylab as plt

        NORD = self.clbr.pattern.shape[1]

        fig, axs = plt.subplots(NORD, 1)
        plt.subplots_adjust(
            left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=0.0
        )
        [ax.set_xticks([]) for ax in axs[:-1]]
        [ax.set_yticks([]) for ax in axs[:-1]]
        [ax.set_xlim([0, self.clbr.DIMW]) for ax in axs]

        for o, ax in enumerate(reversed(axs)):
            ax.imshow(self.order_image(frame, o), aspect=aspect)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            # ax.axis('off')
            ax.text(
                0,
                0.5,
                o,
                transform=ax.transAxes,
                color="#fff83a",
                va="center",
                ha="left",
            )


class Calibrations:
    """ """

    folder = join("../calibration_files")
    dv = 8  # dv*2 + 1 - width of the order in pixels

    def __init__(self, folder=folder):
        """ """
        self.folder = folder
        self.filenames = {
            "orders": "pattern.txt",
            "wavelength": "Th_wavelength.txt",
            "sphr": "absolute_20170613_b8_0.2_v2.sif",
            "bkgr": "absolute_20170613_b8_0.2_bkg.sif",
            "integral": "integrating_sphere.txt",
        }
        # print(os.listdir(folder))

    def start(self):
        """ Prepare calibrations:
            1. Cutting patterns for all diffraction orders
            2. Wavelength calibration
            3. Absolute calibration
        """
        self.load_pattern()
        self.load_sphere()
        self.make_cutting_masks()

        self.wavelength_calibration()
        self.calculate_order_borders()

        self.absolute_calibration()

    # =============================
    # PATTERN
    # =============================
    def load_sphere(self):
        """
        Load sphere image to
        1. Get image dimentions to make cutting patterns
        2. Use it in absolute calibraiton
        """
        self.sphr = EchelleImage(join(self.folder, self.filenames["sphr"]), clbr=self,)

        self.bkgr = EchelleImage(join(self.folder, self.filenames["bkgr"]), clbr=self,)

        self.DIMW, self.DIMO = self.sphr.info["size"] * np.array(
            [self.sphr.info["xbin"], self.sphr.info["ybin"],]
        )

        bDIMW, bDIMO = self.bkgr.info["size"] * np.array(
            [self.bkgr.info["xbin"], self.bkgr.info["ybin"],]
        )

        if self.DIMW != bDIMW or self.DIMO != bDIMO:
            raise ValueError(
                "Wrong background size: {}x{} vs {}x{}".format(
                    self.DIMW, self.DIMO, bDIMW, bDIMO
                )
            )

    def load_pattern(self):
        """ Read file with coordinats for each difraction order of 
        the Echelle spectrometer
        """
        pth = join(self.folder, self.filenames["orders"])
        self.pattern = np.loadtxt(pth, dtype="int")

    def make_mask(self, ordind, show=False, **kws):
        """
        converts linear coordinates into 2d mask to mask the image
        """
        dv = kws.get("dv", self.dv)

        l = self.pattern[:, ordind]
        cc = np.arange(-dv, dv + 1, 1)
        ii = ((np.zeros([self.DIMW, 1]) + cc).T + l).flatten()

        jj = np.repeat(
            np.arange(self.DIMW)[np.newaxis, ...], dv * 2 + 1, axis=0,
        ).flatten()

        mask = (ii.astype(int, copy=False), jj.astype(int, copy=False))
        if not show:
            return mask
        else:
            pp = np.zeros((self.DIMW, self.DIMO), dtype=bool)
            pp[mask] = True
            return pp

    def make_cutting_masks(self, **kws):
        """ make cutting masks for each diffraction order"""
        self.dv = kws.get("dv", self.dv)
        self.cutting_masks = [
            self.make_mask(i, dv=self.dv) for i in range(self.pattern.shape[1])
        ]

        # check if pixel in the mask is out of top bound (for 28th order, CMOS)
        ordrs = [
            i for i, j in enumerate(self.cutting_masks) if j[0].max() > self.DIMO - 1
        ]

        for o in ordrs:
            ind = self.cutting_masks[o][0] > self.DIMO - 1
            maxlam = self.cutting_masks[o][1][ind].min()
            ind1 = self.cutting_masks[o][1] < maxlam
            self.cutting_masks[o] = (
                self.cutting_masks[o][0][ind1],
                self.cutting_masks[o][1][ind1],
            )
            # dimw = int(cmsk[0].shape[0]/self.dv)

    def show_masks(self, **kws):
        """ show all masks """
        dv = kws.get("dv", self.dv)
        self.masks = np.array(
            [self.make_mask(i, show=True, dv=dv) for i in range(self.pattern.shape[1])]
        )

    # =============================
    # WAVELENGTH
    # =============================
    def wavelength_calibration(self, **kws):
        """Read calibration data, containing list of identified lines from
        lamp spectra in format [order pixfrom pixto center wavelength]
        """
        fpth = join(self.folder, self.filenames["wavelength"])
        # Columns: order, from, to, center, wavelength
        w = np.loadtxt(fpth, skiprows=2, usecols=(0, 1, 2, 3, 4))

        ords = np.unique(w[:, 0]).astype(int)
        msks = [np.where(w[:, 0] == o) for o in ords]
        ab = lambda x: 1 if x < 3 else 2
        pwrs = [ab(len(i[0])) for i in msks]
        x = np.arange(self.DIMW)

        self.order_wavel = np.array(
            [
                np.poly1d(np.polyfit(w[m][:, 3], w[m][:, 4], p))(x)
                for m, p in zip(msks, pwrs)
            ]
        )
        if kws.get("rtrn", False):
            return self.order_wavel
        if kws.get("plot", False):
            import matplotlib.pylab as plt

            _ = [plt.plot(self.order_wavel[i]) for i in ords]
            _ = [plt.plot(w[m][:, 3], w[m][:, 4], "o") for m in msks]

        """ For CMOS camera the 28th order has H-gamma, but only part of
        the order is on the sensor. So the length of this last order is shorter.
        To compensate for this and make all spectra the same length,
        np.nan values should be put in the wavelength vector and
        in each specrum for all (only 28th for now) orders which are out of the
        boundary.
        """
        self.sphr.calculate_order_spectra()

        spectra_shapes = np.array([s.shape[0] for s in self.sphr.order_spectra[0]])
        lambda_shapes = np.array([ow.shape[0] for ow in self.order_wavel])
        wrong_shapes = spectra_shapes - lambda_shapes
        wrong_shapes = np.where(wrong_shapes != 0)[0]

        froms = []
        for o in wrong_shapes:
            newshape = spectra_shapes[o]
            self.order_wavel[o][newshape:] = np.nan
            froms.append(newshape)

        self.orders_bad_shape = wrong_shapes
        self.orders_bad_froms = np.array(froms)

    def print_bad_shapes(self):
        """ Print orders with bad shapes and the indices where the np.nan starts
        """
        for o, i in zip(self.orders_bad_shape, self.orders_bad_froms):
            print(o, i)

    def print_wavelength(self):
        """ Print pandas dataframe to show the nan values in the
        wavelength calibration to confirm the dimentions.

        Usage:

        cb.print_wavelength().style.highlight_null(null_color='red')

        where cb -> Calibrations class
        """
        import pandas as pd

        np.set_printoptions(precision=2)
        d = pd.DataFrame(self.order_wavel)
        i = 2173
        return d.loc[[0, 1, 26, 27, 28], range(i, i + 12)].round(2)

    def calculate_order_borders(self):
        """ Calculate mask for the orders from a Integrating Sphere Image
        sphr - EchelleImage('IntegratingSphere.SIF')
        """

        def isc(x1, y1, x2, y2):
            """ For two neighbor orders with lambda \propto -pix
            find end index for x1 and start index for x2 so they will 
            not intersect
            x1 < ind1 and x2 >= ind2
            """
            from scipy.interpolate import interp1d

            # remove np.nan from vectors
            i1 = ~np.isnan(x1)
            i2 = ~np.isnan(x2)
            x1 = x1[i1]
            y1 = y1[i1]
            x2 = x2[i2]
            y2 = y2[i2]

            f1 = interp1d(x1, y1, 1)
            f2 = interp1d(x2, y2, 1)
            x = np.linspace(*sorted([x1[0], x1[-1], x2[0], x2[-1]])[1:-1], len(x2))
            x0 = x[np.argmin(np.abs(f1(x) - f2(x)))]
            ind2 = np.where(x2 <= x0)
            ind1 = np.where(x1 > x0)

            return ind1[0][-1] + 1, ind2[0][0]

        self.sphr.correct_order_shapes()

        nord = self.order_wavel.shape[0]

        brdr = [
            isc(
                remove_npnans(self.order_wavel[o]),
                remove_npnans(self.sphr.order_spectra[0, o]),
                remove_npnans(self.order_wavel[o + 1]),
                remove_npnans(self.sphr.order_spectra[0, o + 1]),
            )
            for o in range(nord - 1)
        ]

        brdrs = np.append(np.insert(np.array(brdr), 0, 0), [self.DIMW]).reshape(nord, 2)

        r = np.arange(self.DIMW)[:, None]
        b = (brdrs[:, 0] <= r) & (r <= brdrs[:, 1])
        self.order_borders = b.T
        # self.wavelength = remove_npnans(self.order_wavel[self.order_borders])
        self.wavelength = self.order_wavel[self.order_borders]

    def plot_order_borders(self):
        """ Plot order spectra vs wavelength and show how the borders are
        selected
        """
        import matplotlib.pylab as plt

        clrs = ["#498fff", "#a3e21b"]

        ax = plt.gca()

        for o in range(self.order_wavel.shape[0]):
            lam = self.order_wavel[o]
            sig = self.sphr.order_spectra[0][o]

            plt.plot(
                lam, sig, c="gray",
            )

            try:
                ind = self.order_borders[o]
                plt.plot(
                    lam[ind], sig[ind], c=clrs[o % 2],
                )
            except:
                pass

            ax.text(
                self.order_wavel[o][1500],
                self.sphr.order_spectra[0][o][1500] * 1.5,
                o,
                va="center",
                ha="center",
                color=clrs[o % 2],
            )

        ax.set_yscale("log")
        plt.gcf().set_size_inches([20, 5])

    # =============================
    # Absolute
    # =============================
    def absolute_calibration(self):
        """ Read the integrating sphere wavelength characteristic
        """
        from scipy.interpolate import interp1d
        from scipy.constants import speed_of_light, Planck

        data = np.loadtxt(join(self.folder, self.filenames["integral"]))
        self.integral = interp1d(data[:, 0] * 1000, data[:, 1], 3)

        self.sphr.calculate_spectra()

        self.bkgr.calculate_order_spectra()
        self.bkgr.correct_order_shapes()
        self.bkgr.calculate_spectra()

        y = self.sphr.spectra[0] - self.bkgr.spectra[0]
        x = self.sphr.wavelength
        # np.nans Warning here:
        wmsr = (
            self.integral(x) * self.sphr.info["ExposureTime"] / y * 1e-2
        )  # W/(m2 sr nm)
        wm = wmsr * 4 * np.pi  # W/(m2 sr nm)
        phmsr = (
            wmsr * x * 1e-9 / (speed_of_light * Planck)
        )  # convert W (or J/s) to Nph/s

        self.absolute = {"wmsr": wmsr, "wm": wm, "phmsr": phmsr}


class Spectrum:

    """ Echelle spectra converted from EchelleImage
    Contains image info, input - EchelleImage.order_image(....,sm=True)

    First frame, number 0, is normally used for background subtraction.
    To use arbitrary frame (or another image), an update is required
    """

    def __init__(self, image):
        """ """
        try:
            self.wavelength = image.wavelength
        except AttributeError:
            image.calculate_order_spectra()
            image.calculate_spectra()
        finally:
            self.wavelength = image.wavelength

        self.info = image.info
        self.fpth = image.fpth
        if image.spectra.shape[0] > 1:
            self.counts = image.spectra - image.spectra[0]
        else:
            self.counts = image.spectra

        self.absolute = image.clbr.absolute
        image = None

        # Flip wavelength, it is from high to low by default for black Echelle
        self.wavelength = np.flip(self.wavelength)
        self.counts = np.flip(self.counts, axis=1)
        self.absolute = {i: np.flip(j) for i, j in self.absolute.items()}

        self.wm = self.counts * self.absolute["wm"] / self.info["ExposureTime"]
        self.wmsr = self.counts * self.absolute["wmsr"] / self.info["ExposureTime"]
        self.phmsr = self.counts * self.absolute["phmsr"] / self.info["ExposureTime"]

        self.savepath = None
        self.shotnumber = None
        self.trigdelay = "none"
        self.path_output = None

        self.saveunits = "wm"
        self.units_kws = ["counts", "wm", "wmsr", "phmsr"]
        units_names = ["counts", "W/(m2 nm)", "W/(m2 sr nm)", "Nph./(m2 sr)"]
        spectra_to_save = [self.counts, self.wm, self.wmsr, self.phmsr]
        self.units_names = {i: j for i, j in zip(self.units_kws, units_names)}
        self.spectra_to_save = {i: j for i, j in zip(self.units_kws, spectra_to_save)}

    def plot_frame_incounts(self, frame):
        import matplotlib.pylab as plt

        plt.plot(self.wavelength, self.counts[frame])

    def save(self):
        """ save spectra
        """
        import datetime
        import time

        if self.shotnumber is None:
            self.shotnumber = os.path.basename(self.fpth)[:-4]

        frames = np.arange(self.info["NumberOfFrames"])

        vnames = "'" + "','".join(str(i) for i in frames) + "'"
        vunit = "'" + "','".join(self.units_names[self.saveunits] for i in frames) + "'"
        fmt = ["%.6f"] + ["%.6e" for i in frames]
        date = datetime.datetime.fromtimestamp(time.time())
        date = date.strftime("%m/%d/%Y %H:%M")
        header = header_template.format(
            shot=self.shotnumber,
            date=date,
            size=len(frames),
            nval=len(frames),
            vnames=vnames,
            vunit=vunit,
            cycletime=self.info["CycleTime"],
            exposure=self.info["ExposureTime"],
            trigdelay=self.trigdelay,
        )

        pth = join(
            self.path_output, "{}_{}.txt".format(self.shotnumber, "echelle_spec")
        )

        data = np.vstack((self.wavelength, self.spectra_to_save[self.saveunits]))
        try:
            np.savetxt(pth, data.T, fmt=fmt, delimiter=", ", header=header, comments="")
        except Exception as err:
            print("failed to save spectra\n{}".format(err))


header_template = """# [Parameters]
# Name = 'Echelle Spectra'
# ShotNo = {shot}
# Date = '{date}'
#
# DimNo = 1
# DimName = 'wavelength'
# DimSize = {size}
# DimUnits = 'nm'
#
# ValNo = {nval}
# ValName = {vnames}
# ValUnit = {vunit}
# [Comments]
# time = {trigdelay} + frameNo*{cycletime} (s)
# exposure = {exposure} (s)
#
# [data]"""

if __name__ == "__main__":
    pass
