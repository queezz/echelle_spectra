import numpy as np
from lmfit.lineshapes import gaussian


def baseline_als(y, lam, p, niter=10):
    """Asymmetric Least Squares Smoothing" by P. Eilers and H. Boelens
    y - signal, lam - parameter (float or array), p - weight (float, array?)
    """
    from scipy.sparse import csc_matrix, spdiags
    from scipy.sparse.linalg import spsolve

    L = len(y)
    D = csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def banddata(band, frame, s, **kws):
    """Extract data from the spctra at the frame,
    x = band.bounds
    """
    wl = s.wavelength
    dx = 0.1  # nm
    m = (wl >= band.bounds[0]) & (wl <= band.bounds[-1])
    m1 = (wl >= band.bounds[0] - dx) & (wl <= band.bounds[-1] + dx)

    x = s.wavelength[m]
    x1 = s.wavelength[m1]
    p = x.argsort()
    p1 = x1.argsort()
    x = x[p]
    x1 = x1[p1]
    # s.wm - W/(m2 nm), s.wmsr - W/(m2 sr nm) s.phmsr - Nph/(m2 sr)
    y = s.wm[frame][m][p]
    y1 = s.wm[frame][m1][p1]
    doblin = kws.get("doblin", True)
    if doblin:
        bline = baseline_als(y1, 1e7, 1e-5, 5)
        ym = bline.mean()
    else:
        di = max(5, int(len(x1) * 0.01))
        c1 = y1[:di].mean()
        c2 = y1[-di:].mean()
        d1 = x1[:di].mean()
        d2 = x1[-di:].mean()
        a = (c2 - c1) / (d2 - d1)
        b = c1 - a * x1[0]
        # ym = (c1+c2)/2.
        ym = a * x + b

    return x, y - ym


class EB:
    """Mimic EmissionBand object for intensity calculations
    bounds - interval where intesities are wanted
    name - name of this band
    method calculate intensity
      no baseline calculation, no gaussian fitting
      simply subtract mean values from left and right and numerically integrate
    """

    def __init__(self, **kws):
        self.bounds = kws.get("bounds", None)
        self.name = kws.get("name", None)

    def intensity(self, s):
        """s - Spectrum Object from echelle.py
        (spectra for all frames, caiblrated)
        """

    HTML_REPORT_TEMP = (
        "<div class='report'>"
        "<span style='font-size: 20px; color: red;'>{}</span><br>"
        "<span style='font-size: 14px; color: #3b9cd4;'>"
        "Boundaries, nm<br>"
        "</span>"
        "{}<br>"
        "</div>"
    )

    def report_html(self):
        """Generate band report in html friendly format"""
        return self.HTML_REPORT_TEMP.format(
            self.name,
            f"{self.bounds[0]:.2f} - {self.bounds[-1]:.2f}",
        )


class EmissionBand(object):
    "Transition data for CR models"

    def __init__(
        self, name, wavelengths, configs, Js, Ak, energ_low, energ_high, **kws
    ):
        """ """
        self.name = name
        self.cw = wavelengths
        self.configs = configs
        self.Js = Js
        self.terms = kws.get("terms")
        self.Ak = Ak
        self.energ_low = energ_low
        self.energ_high = energ_high
        self.wavelengths = wavelengths
        self.nw = len(self.cw)  # number of lines
        self.stat_weights = self.calc_stat_weights()
        self.bounds = np.array([np.mean(self.cw) - 0.3, np.mean(self.cw) + 0.3])
        self.average_Ak = self.calc_average_Ak()
        self.make_model()

    def copy(self):
        from copy import deepcopy

        return deepcopy(self)

    def stat_weight(self, ind):
        """calculate statweights from J:
        wg = 2*J+1 if J != *
        wg = (2*S+1)*(2*l+1), 2*S+1 - spin, l - angular quantum number
        """
        js = self.Js[ind].split("-")[1]
        if "*" not in js:
            wg = 2 * self.fracfloat(js) + 1
        else:
            lnum = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5, "i": 6, "k": 7}
            l = lnum[self.configs[ind][-1]]
            wg = 2 * (2 * l + 1)
        return wg

    def calc_stat_weights(self):
        return np.array([self.stat_weight(i) for i in range(self.nw)])

    def calc_average_Ak(self):
        """Calculate average Ak for a band"""
        return sum(self.Ak * self.stat_weights) / sum(self.stat_weights)

    def fracfloat(self, frac_str):
        """Convert fraction-like J to float
        1/2 -> 0.5, 1 1/2 -> 1.5"""
        try:
            return float(frac_str)
        except ValueError:
            num, denom = frac_str.split("/")
            try:
                leading, num = num.split(" ")
                whole = float(leading)
            except ValueError:
                whole = 0
            frac = float(num) / float(denom)
            return whole - frac if whole < 0 else whole + frac

    REPORT_TEMP = (
        "{}\nlambdas,nm\n{}\n"
        "------------------------\n"
        "Eterm = {:.2f}     [eV]\n"
        "<Ak>  = {:.2e}\n"
        "<wg>  = {:.0f}\n"
        "------------------------\n"
        "wg\n{}\nAk\n{}\nConfiguration\n{}\nJ-J{}\n"
        "Boundaries\n{}"
    )

    def report(self):
        """Generate band report as plain text"""
        report_text = self.REPORT_TEMP.format(
            self.name,
            self.cw,
            np.mean(self.energ_high),
            self.average_Ak,
            sum(self.stat_weights),
            self.stat_weights,
            ["%.2e" % i for i in self.Ak],
            self.configs,
            self.Js,
            f"{self.bounds[0]:.2f} - {self.bounds[-1]:.2f} (nm)",
        )
        return report_text

    HTML_REPORT_TEMP = (
        "<div class='report'>"
        "<span style='font-size: 20px; color: red;'>{}</span><br>"
        "central wavelengths, nm<br>"
        "{}<br>"
        "------------------------<br>"
        "Eterm = {:.2f} [eV]<br>"
        "&lt;Ak&gt; = {:.2e}<br>"
        "&lt;wg&gt; = {:.0f}<br>"
        "------------------------<br>"
        "wg<br>"
        "{}<br>"
        "Ak<br>"
        "{}<br>"
        "<span style='font-size: 14px; color: #3b9cd4;'>"
        "Configuration<br>"
        "</span>"
        "{}<br>"
        "<span style='font-size: 14px; color: #3b9cd4;'>"
        "J-J<br>"
        "</span>"
        "{}<br>"
        "<span style='font-size: 14px; color: #3b9cd4;'>"
        "Boundaries, nm<br>"
        "</span>"
        "{}<br>"
        "</div>"
    )

    def report_html(self):
        """Generate band report in html friendly format"""
        return self.HTML_REPORT_TEMP.format(
            self.name,
            self.cw,
            np.mean(self.energ_high),
            self.average_Ak,
            sum(self.stat_weights),
            self.stat_weights,
            ["%.2e" % i for i in self.Ak],
            self.configs,
            self.Js,
            f"{self.bounds[0]:.2f} - {self.bounds[-1]:.2f}",
        )

    def make_model(self):
        """Make lmfit model to fit the experimental data"""
        from lmfit.models import GaussianModel
        from lmfit.lineshapes import gaussian

        Np = len(self.cw)
        dlam = [i - self.cw[0] for i in self.cw]
        gauss = []
        for i in range(Np):
            gauss.append(GaussianModel(prefix="g{}_".format(i)))
            if i == 0:
                pars = gauss[i].make_params()
            else:
                pars.update(gauss[i].make_params())

        pars["g0_center"].set(self.cw[0])
        # Constrain center wavelength
        for i in range(1, Np):
            pars["g{}_center".format(i)].set(expr="g0_center+{}".format(dlam[i]))

        # Amplitueds, set minimum amplitude 0
        # [pars['g{}_amplitude'.format(i)].set(min=0) for i in range(Np)]
        # Ratio of the lines in the band proportional to Ak*wi
        ratio = self.Ak * self.stat_weights
        ratio = ratio / ratio[0]
        if "H-a" in self.name or "H-be" in self.name or "H-g" in self.name:
            ratio = [6.0e-2, 1]
        [
            pars.add("ratio{}".format(i), value=j, vary=False)
            for i, j in zip(range(Np), ratio)
        ]
        # delta
        pars.add("delta", value=1, vary=True, min=0)
        for i in range(Np):
            # pars.add('delta{}'.format(i), value = 1, min=1.2,max=2, vary=True)
            pars["g{}_amplitude".format(i)].set(expr="delta*ratio{}".format(i))

        # Sigma
        pars["g0_sigma"].set(0.06 / 2.355, min=0.03 / 2.355, max=0.2 / 2.355)
        [pars["g{}_sigma".format(i)].set(expr="g0_sigma") for i in range(1, Np)]
        self.pars = pars
        self.model = np.sum(gauss)

    def fitb(self, frame, s, **kws):
        """fit band"""
        doblin = kws.get("doblin", True)
        x, y = banddata(self, frame, s, doblin=doblin)
        ymax = y.max()
        # Detect if there is no line peak and the intensity is close
        # to the noise level, for W/(m2 nm)
        detect = kws.get("detect", False)
        if detect:
            if y.max() < y.std() + y.mean() + 0.11 and y.max() < 0.3:
                return {"fit": None, "exp": [x, y], "out": None}

        out = self.model.fit(y / ymax, self.pars, x=x)

        v = out.best_values
        xx = np.linspace(x[0], x[-1], len(x) * 10)
        ys = np.array(
            [
                gaussian(
                    xx,
                    v["g{}_amplitude".format(i)] * ymax,
                    v["g{}_center".format(i)],
                    v["g{}_sigma".format(i)],
                )
                for i in range(len(self.cw))
            ]
        )
        return {"fit": [xx, ys], "exp": [x, y], "out": out}


class FitResult:
    """ """

    def __init__(self, name):
        """ """
        self.name = name
        self.frames = None
        self.spectra = None
        self.x = None

        self.x_fit = None
        self.spectra_fit = None
        self.fit_out = None

        self.intensities = None
        self.intensities_fit = None
        self.fit_frames = None

    def set_frames(self, frames):
        self.frames = frames

    def set_info(self, n):
        """set number of frames from SIF image"""
        self.NumberOfFrames = n
        # make a mask for the full frame range, True if intensities
        #  are calculated for the frame
        if self.frames is None:
            msg = """self.frames are None, but should be an array of frames 
                     used to calculate intensities"""
            raise ValueError(msg)

        self.mask_frames = np.zeros(n, dtype=bool)
        self.mask_frames[self.frames] = True

    def init_spectra(self):
        """ """
        if self.frames is None:
            msg = """self.frames are None, but should be an array of frames 
                         used to calculate intensities"""
            raise ValueError(msg)

        self.spectra = []
        self.spectra_fit = {i: None for i in self.frames}
        self.fit_out = {}

    def integrate_spectra(self):
        """ """
        intensities = np.array([(np.trapz(s, self.x)) for s in self.spectra])

        self.intensities = np.empty(self.NumberOfFrames)
        self.intensities[:] = np.nan
        self.intensities[self.mask_frames] = intensities

    def integrate_fit(self):
        """ """
        self.intensities_fit = np.empty(self.NumberOfFrames)
        self.intensities_fit[:] = np.nan

        for frame in list(self.spectra_fit):
            if self.spectra_fit[frame] is not None:
                self.intensities_fit[frame] = np.trapz(
                    self.spectra_fit[frame].sum(axis=0), self.x_fit
                )

    def fill_nans(self):
        """fill missing values in the fitted intensities from
        numerically integrated values of unfitted data
        """
        self.intensities_combined = self.intensities_fit
        m = np.isnan(self.intensities_fit)
        self.intensities_combined[m] = self.intensities[m]
