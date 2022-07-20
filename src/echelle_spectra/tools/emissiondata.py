from . import emissionbands as eb

# ===============================================================================
# Define Emission Bands
# ===============================================================================

c444 = eb.EmissionBand(
    "CIV-444",
    # wavelengths = [444.14,444.255,444.279],# TALL
    wavelengths=[444.034, 444.150, 444.174],  # NIST
    configs=["1s2.5p-1s2.6d", "1s2.5p-1s2.6d", "1s2.5p-1s2.6d"],
    Js=["1/2-3/2", "3/2-5/2", "3/2-3/2"],
    Ak=[8.72e7, 1.05e8, 1.74e7],
    energ_low=[55.650849, 55.651717, 55.651717],
    energ_high=[58.442404, 58.442553, 58.44244],
)
c465 = eb.EmissionBand(
    "CIV-465",
    # wavelengths = [465.886,465.95],
    wavelengths=[465.751, 465.846],  # pfr
    # wavelengths = [465.751, 465.846,4660.4], #test
    # wavelengths = [465.86,465.95],
    configs=["1s2.5f-1s2.6g", "1s2.5g-1s2.6h"],
    Js=["*-*", "*-*"],
    Ak=[2.83e8, 4.21e8],
    energ_low=[55.785632, 55.786116],
    energ_high=[58.446892, 58.447016],
)
c547 = eb.EmissionBand(
    "CIV-547",
    # wavelengths = [547.15,547.26,547.26,547.27,547.27,547.27,547.29],
    wavelengths=[547.00, 547.11, 547.11, 547.11, 547.12, 547.12, 547.14],
    configs=[
        "1s2.7f-1s2.10g",
        "1s2.7h-1s2.10i",
        "1s2.7g-1s2.10h",
        "1s2.7i-1s2.10k",
        "1s2.7i-1s2.10h",
        "1s2.7h-1s2.10g",
        "1s2.7g-1s2.10f",
    ],
    Js=["*-*", "*-*", "*-*", "*-*", "*-*", "*-*", "*-*"],
    Ak=[1.82e7, 1.86e7, 2.02e7, 1.16e7, 5.52e4, 2.42e5, 6.73e5],
    energ_low=[
        60.050765,
        60.051248,
        60.051224,
        60.051273,
        60.051273,
        60.051248,
        60.051224,
    ],
    energ_high=[
        62.31675,
        62.316788,
        62.316775,
        62.316788,
        62.316775,
        62.31675,
        62.316651,
    ],
)
c580 = eb.EmissionBand(
    "CIV-580",
    wavelengths=[580.131],
    configs=["1s2.3s-1s2.3p"],
    Js=["1/2-3/2"],
    Ak=[3.17e7],
    energ_low=[37.548356],
    energ_high=[39.684927],
)
# First line from this CIV-706 band is removed to fit the LHD Echelle data
c706 = eb.EmissionBand(
    "CIV-706",
    wavelengths=[706.42, 706.43, 706.44, 706.5, 706.51, 706.55, 706.87],
    configs=[
        "1s2.7g-1s2.9h",
        "1s2.7h-1s2.9i",
        "1s2.7i-1s2.9h",
        "1s2.7i-1s2.9k",
        "1s2.7h-1s2.9g",
        "1s2.7g-1s2.9f",
        "1s2.7f-1s2.9d",
    ],
    Js=["*-*", "*-*", "*-*", "*-*", "*-*", "*-*", "*-*"],
    Ak=[3.18e7, 3.37e7, 1.16e5, 2.69e7, 4.46e5, 1.13e6, 2.41e6],
    energ_low=[
        60.051224,
        60.051248,
        60.051273,
        60.051273,
        60.051248,
        60.051224,
        60.050765,
    ],
    energ_high=[
        61.80632,
        61.80632,
        61.80632,
        61.806183,
        61.806121,
        61.805997,
        61.804758,
    ],
)
""" All lines for C706, but the 706.32 is always missing from LHD Echelle, trying to fit without it
c706 = eb.EmissionBand('CIV-706',
              wavelengths = [706.32,706.42,706.43,706.44,706.5,706.51,706.55,706.87],
              configs = ['1s2.7f-1s2.9g','1s2.7g-1s2.9h','1s2.7h-1s2.9i','1s2.7i-1s2.9h',
                         '1s2.7i-1s2.9k','1s2.7h-1s2.9g','1s2.7g-1s2.9f','1s2.7f-1s2.9d'],
              Js = ['*-*','*-*','*-*','*-*','*-*','*-*','*-*','*-*'],
              Ak = [2.64e7,3.18e7,3.37e7,1.16e5,2.69e7,4.46e5,1.13e6,2.41e6],
              energ_low = [60.050765,60.051224,60.051248,60.051273,
                           60.051273,60.051248,60.051224,60.050765],
              energ_high = [61.806121,61.80632,61.80632,61.80632,
                            61.806183,61.806121,61.805997,61.804758]
              )
"""  # """
c772 = eb.EmissionBand(
    "CIV-772",
    # wavelengths = [772.61,772.8,772.85], #,772.87],
    wavelengths=[772.430, 772.591, 772.621],
    configs=["1s2.6f-1s2.7g", "1s2.6g-1s2.7h", "1s2.6h-1s2.7i"],  #'1s2.6h-1s2.7g'],
    Js=["*-*", "*-*", "*-*"],  #'*-*'],
    Ak=[9.64e07, 1.36e08, 1.90e08],  # 8.60E+05],
    energ_low=[58.446471, 58.446892, 58.447016],  # 58.447016],
    energ_high=[60.051224, 60.051248, 60.051273],  # ,60.051224]
)

# CII
c515 = eb.EmissionBand(
    "CII-515",
    wavelengths=[515.1085],
    configs=["2s2p3s-2s2p3p"],
    Js=["5/2-3/2"],
    Ak=[4.16e7],
    energ_low=[20.709788],
    energ_high=[23.116071],
)
# CIII
c464 = eb.EmissionBand(
    "CIII-464",
    wavelengths=[464.7418],
    configs=["1s2.2s3s-1s2.2s3p"],
    Js=["1-2"],
    Ak=[7.26e7],
    energ_low=[29.534648],
    energ_high=[32.201709],
)

xe461 = eb.EmissionBand(
    "Xe-461",
    wavelengths=[461.18882],
    configs=["5p4.6s-5p2.7p"],
    Js=["2-1"],
    Ak=[1e6],  # NIST has no data!
    energ_low=[8.315316],
    energ_high=[11.0029229],
)

xe462 = eb.EmissionBand(
    "Xe-462",
    wavelengths=[462.42756],
    configs=["5p4.6s-5p2.7p"],
    Js=["2-2"],
    Ak=[7.2e5],
    energ_low=[8.3153160],
    energ_high=[10.9687844],
)

xe467 = eb.EmissionBand(
    "Xe-467",
    wavelengths=[467.12258],
    configs=["5p4.6s-5p2.7p"],
    Js=["2-3"],
    Ak=[1e6],
    energ_low=[8.3153160],
    energ_high=[10.9687844],
)

# Helium lines
he447 = eb.EmissionBand(
    "He-447",
    wavelengths=[447.150],
    configs=["1s2p-1s4d"],
    Js=["2-3"],
    terms=["2.3P-4.3D"],
    Ak=[0.24574e8],
    energ_low=[20.96408703],
    energ_high=[23.73609031],
)

he492 = eb.EmissionBand(
    "He-492",
    wavelengths=[492.193],
    configs=["2p-4d"],
    Js=["*-*"],
    terms=["2.1P-4.1D"],
    Ak=[0.19855e8],
    energ_low=[0],
    energ_high=[0],
)

he501 = eb.EmissionBand(
    "He-501",
    wavelengths=[501.568],
    configs=["1s2s-1s3p"],
    Js=["0-1"],
    terms=["2.1S-3.1P"],
    Ak=[0.13368e8],
    energ_low=[20.61577496],
    energ_high=[23.08701866],
)

he504 = eb.EmissionBand(
    "He-504",
    wavelengths=[504.774],
    configs=["1s2p-1s4s"],
    Js=["1-0"],
    terms=["2.3S-3.3P"],
    Ak=[0.06768e8],
    energ_low=[21.21802284],
    energ_high=[23.67357071],
)

he587 = eb.EmissionBand(
    "He-587",
    wavelengths=[587.566],
    configs=["1s2p-1s3d"],
    Js=["2-3"],
    terms=["2.3P-3.3D"],
    Ak=[0.70693e8],
    energ_low=[20.96408703],
    energ_high=[23.07365084],
)

he667 = eb.EmissionBand(
    "He-667",
    wavelengths=[667.815],
    configs=["1s2p-1s3d"],
    Js=["1-2"],
    terms=["2.1P-3.1D"],
    Ak=[0.63676e8],
    energ_low=[21.21802284],
    energ_high=[23.07407493],
)

he706 = eb.EmissionBand(
    "He-706",
    wavelengths=[706.525],
    configs=["1s2p-1s3s"],
    Js=["2-1"],
    terms=["2.3P-3.3S"],
    Ak=[0.27849e8],
    energ_low=[20.96408703],
    energ_high=[22.71846655],
)

he728 = eb.EmissionBand(
    "He-728",
    wavelengths=[728.135],
    configs=["2p-3s"],
    Js=["*-*"],
    terms=["2.1P-3.1S"],
    Ak=[0.18291e8],
    energ_low=[21.21802284],
    energ_high=[22.92031749],
)

# H and D alpha
halpha = eb.EmissionBand(
    "H-alpha",
    wavelengths=[656.279, 656.098629],
    configs=["2s-3s", "2s-3s"],
    Js=["*-*", "*-*"],
    terms=["2-3", "2-3"],
    Ak=[4.4101e07, 4.4114e07],
    energ_low=[12.0875051, 12.09079461],
    energ_high=[10.1988357, 10.20160487],
)

# Dbeta 486.00013 nm (NIST)
hbeta = eb.EmissionBand(
    "H-beta",
    wavelengths=[486.135, 486.00013],
    configs=["2s-4s", "2s-4s"],
    Js=["*-*", "*-*"],
    terms=["2-4", "2-4"],
    Ak=[8.4193e06, 8.4217e06],
    energ_low=[12.7485392, 12.75200754],
    energ_high=[10.1988357, 10.20160487],
)

hgamma = eb.EmissionBand(
    "H-gamma",
    wavelengths=[434.0472, 433.92833],
    configs=["2s-5s", "2s-5s"],
    Js=["*-*", "*-*"],
    terms=["2-5", "2-5"],
    Ak=[2.5304e06, 2.5311e06],
    energ_low=[13.0545016, 13.05805386],
    energ_high=[10.1988357, 10.20160487],
)

hdelta = eb.EmissionBand(
    "H-delta",
    wavelengths=[410.1734],
    configs=["2s-6s"],
    Js=["*-*"],
    terms=["2-6"],
    Ak=[9.7320e05],
    energ_low=[13.22070378],
    energ_high=[10.1988357],
)

crhe_names = [
    "HeI-447",
    "HeI-492",
    "HeI-501",
    "HeI-504",
    "HeI-587",
    "HeI-667",
    "HeI-706",
    "HeI-728",
]
crhe_bands = {
    "HeI-447": he447,
    "HeI-492": he492,
    "HeI-501": he501,
    "HeI-504": he504,
    "HeI-587": he587,
    "HeI-667": he667,
    "HeI-706": he706,
    "HeI-728": he728,
}
crc_names = ["CIV-580", "CIV-444", "CIV-465", "CIV-772", "CIV-706", "CIV-547"]
crc_bands = {
    "CIV-580": c580,
    "CIV-444": c444,
    "CIV-465": c465,
    "CIV-772": c772,
    "CIV-706": c706,
    "CIV-547": c547,
}

c_bands = {"CII-515": c515, "CIII-464": c464}


def set_bounds(dw=0.5):
    """ Expand bounds of crc bands by dw in both directions
    #dw = .23 - This value was used before I upgraded this code in Jan 2018
    """
    crc_centers = {
        "CIV-444": 444.1,
        "CIV-465": 465.8,
        "CIV-547": 547.05,
        "CIV-580": 580.18,
        "CIV-706": 706.35,
        "CIV-772": 772.62,
    }
    for name in crc_names:
        band = crc_bands[name]
        cw = crc_centers[name]
        band.bounds = [cw - dw, cw + dw]

    c580.bounds = [579.95, 580.35]
    c706.bounds = [706.3, 706.845]
    c464.bounds = [464.55, 464.92]
    c465.bounds = [465.50, 466.10]

    xe467.bounds = [467.0, 467.25]
    xe461.bounds = [461.05, 461.3]
    xe462.bounds = [462.3, 462.55]

    he501.bounds = [501.4, 501.8]
    he706.bounds = [706.4, 706.75]
    he728.bounds = [727.92, 728.4]

    halpha.bounds = [655.8, 656.6]
    hbeta.bounds = [485.8, 486.35]
    hgamma.bounds = [433.75, 434.2]
    hdelta.bounds = [409.8, 410.4]
    # CH.bounds = [429,431.05]


set_bounds()
