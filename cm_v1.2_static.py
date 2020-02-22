import math

import PySimpleGUI as sg

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasAgg
import matplotlib.backends.tkagg as tkagg
import tkinter as Tk

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from random import random
from PIL import Image

# INITIALIZE GLOBAL VARIABLES

pH = 5
conc_O2 = 9
biofilm = True
conc_Clm = 50
conc_Cl2 = 50
temp_c = 28
fluid_velocity = 10
water_resistivity = 200

pct_Cr = 13.115
pct_Mo = 0.029
pct_N = 0.019
pct_Mn = 0.786
pct_S = 0.017

steel_austenitic = False
steel_duplex = False
steel_ferritic = False

caco3 = 10
ca2p = 10
tds = 10
alkM = 10

conc_SO4 = 10

# Seawater exposure time in years
sea_exposure_time = 10

# BEGIN PITTING FUNCTIONS
##########################################


def func_cathodic_potential():
    # Where conc_O2 is the oxygen concentration in ppm, conc_Clm is the CHLORIDE concentration in ppm
    # biofilm is a boolean, where true = existence of biofilm
    # Function returns the cathodic potential in volts
    if not biofilm:
        cathodic_potential_O2 = 1.23 - (0.33 * math.log10(50 / conc_O2)) - (0.059 * pH)
    elif biofilm:
        cathodic_potential_O2 = 1.23 - (0.33 * math.log10(50 / conc_O2)) - (0.059 * pH) + 0.3
    else:
        cathodic_potential_O2 = 1.23 - (0.33 * math.log10(50 / conc_O2)) - (0.059 * pH) + 0.3

    cathodic_potential_Cl2 = 1.36 + (0.6 * math.log10(conc_Clm / 36))

    return max(cathodic_potential_O2, cathodic_potential_Cl2)


def func_pitting_potential():
    if steel_austenitic:
        k = 16
    elif steel_duplex:
        k = 30
    elif steel_ferritic:
        k = 0
    else:
        k = 16

    pren = pct_Cr + (3.3 * pct_Mo) + (k * pct_N)
    pin = pren - (0.1 * pct_Mn) - (100 * pct_S)     # Pitting Initiation Number

    pitting_potential = (298 / (temp_c + 273)) * ((pin/25) + math.log10(1 + fluid_velocity) -
                                                  0.25 * math.log10((conc_Clm / 36) + 1) + ((pH - 7) / 25))

    return pitting_potential, pren


def func_pitting_induction_time():
    cathodic_potential = func_cathodic_potential()
    pitting_potential, pren = func_pitting_potential()
    driving_voltage = cathodic_potential - pitting_potential

    if driving_voltage > 0:
        pitting_induction_time = math.pow(10, ((pren / (2 * math.log10(conc_Clm))) * (1 - driving_voltage)))
        # pit is the pitting induction time in hours
        return pitting_induction_time
    else:
        return "No pitting"


def func_pitting_corrosion_rate():
    # Returns pitting corrosion rate in mm/year
    oxygen_limiting_current = 10 * 2**((temp_c - 25)/25) * conc_O2 * (1 + math.sqrt(fluid_velocity))
    pitting_corrosion_rate = 20 * math.sqrt((oxygen_limiting_current / 1000) / (water_resistivity + 0.8))

    return pitting_corrosion_rate


# Shaft length is 1.87 m, diameter is 70 mm


# END PITTING FUNCTIONS
##########################################


# BEGIN LANGELIER FUNCTIONS
##########################################
# Source of equations: https://www.mae.gov.nl.ca/waterres/quality/drinkingwater/pdf/calculation_langelier_index.pdf

# caco3 is the alkalinity, in mg/L
# ca2 is the calcium ion concentration, in mg/L
# tds is the total dissolved solids, in mg/L
# temp is the water temperature, in Celsius

def ph_saturation():
    # Determine the Ionic Strength (I), in moles per litre (M), of the water:
    i = (2.5*10**(-5)) * tds

    # Convert temperature to Kelvins
    temp_kelvin = temp_c + 273.15
    # Dialectric constant for water
    d = 78.3
    # a is some constant we need later
    a = (1.82*(10**6)) * ((d*temp_kelvin)**(-1.5))
    # Oxidation number for monovalent ions
    z_mono = 1
    # Oxidation number for divalent ions
    z_di = 2

    # Determine gamma_m, the activity coefficient of monovalent ions (ions that are able to form only one covalent
    # or ionic bond – having only one valence) using the Davies relationship:
    if i <= 0.5:
        log_gamma_m = -a * z_mono**2 * ((math.sqrt(i)/(1 + math.sqrt(i))) - 0.2*i)
    elif i > 0.5:
        log_gamma_m = -a * z_mono**2 * (math.sqrt(i)/(1 + math.sqrt(i)))
    else:
        log_gamma_m = -a * z_mono ** 2 * (math.sqrt(i) / (1 + math.sqrt(i)))

    pk2 = 2902.39/temp_kelvin + 0.02379*temp_kelvin - 6.498
    k2 = 10**(-pk2)

    # Calculate gamma_d, the activity coefficient of divalent ions (ions having two valences):
    log_gamma_d = -a * z_di**2 * (math.sqrt(i)/(1 + math.sqrt(i)))
    gamma_d = 10**(log_gamma_d)

    k2_prime = k2 / gamma_d
    pk2_prime = math.log10(1/k2_prime)

    pks = 0.01183*temp_c + 8.03
    ks = 10**(-pks)
    ks_prime = ks/(gamma_d**2)
    pks_prime = math.log10(1/ks_prime)

    # Convert Ca2+ to moles/liter
    ca2_mol = (ca2p*10**(-3)) / 40
    pca2_mol = math.log10(1/ca2_mol)

    # Convert alkalinity to moles/liter
    alk_mol = (caco3*10**(-3)) / 100

    ph_s = pk2_prime + pca2_mol - pks_prime - math.log10(2 * alk_mol) - log_gamma_m

    print(ph_s)

    return ph_s


def langelier_index():
    ph_s = ph_saturation()
    lsi = pH - ph_s

    return lsi


def func_LSI():
    if temp_c <= 25:
        pHs = 12.65 - 0.0142 * (1.8*temp_c + 32) - math.log10(ca2p) - math.log10(alkM) + 0.1 * math.log(tds)
    elif temp_c > 25:
        pHs = 12.27 - 0.00915 * (1.8*temp_c + 32) - math.log10(ca2p) - math.log10(alkM) + 0.1 * math.log(tds)
    else:
        pHs = 12.27 - 0.00915 * (1.8 * temp_c + 32) - math.log10(ca2p) - math.log10(alkM) + 0.1 * math.log(tds)

    LSI = pH - pHs

    return LSI


def func_LSI_corrosion_rate():
    if func_LSI() <= 0:
        scaling = False
        # Corrosion rate in um/year
        LSI_corrosion_rate = 12 * (conc_O2 + 0.04 * conc_Cl2) * (1 + math.sqrt(fluid_velocity)) * math.pow(2, (temp_c-25)/25)

    elif func_LSI() > 0:
        scaling = True
        # Corrosion rate in um/year
        LSI_porosity = (10 - (func_LSI()**2)) / 10
        LSI_corrosion_rate = 1.2 * (conc_O2 + 0.04 * conc_Cl2) * (10 - (func_LSI()**2)) * math.pow(2, (temp_c-25)/25)

    return LSI_corrosion_rate

# END LANGELIER FUNCTIONS
##########################################

# START SRB FUNCTIONS
##########################################


def func_SRB():
    # Sulphate reducing bacteria corrosion rate, mm/year
    SRB_corrosion_rate = conc_SO4**0.5 - 0.0016 * (temp_c - 38)**2

    return SRB_corrosion_rate

# END SRB FUNCTIONS
##########################################

# BEGIN SEAWATER FUNCTIONS
##########################################

def func_seawater_stagnant():
    stagnant_seawater_rate = (12 * conc_O2 * math.pow(2, (temp_c-25)/25)) / (math.sqrt(sea_exposure_time))
    return stagnant_seawater_rate

def func_seawater_flowing():
    flowing_seawater_rate = 12 * conc_O2 * (1 + math.sqrt(fluid_velocity)) * math.pow(2, (temp_c-25)/25)
    return flowing_seawater_rate

def func_seawater_galvanic():
    seawater_galvanic_rate = 1.4 * math.sqrt(conc_O2 * math.pow(2, (temp_c-25)/25) * (1 + math.sqrt(fluid_velocity)))
    return seawater_galvanic_rate

# END SEAWATER FUNCTIONS
##########################################

# BEGIN SIMULATION FUNCTIONS
##########################################

poisson_constant = 10

def func_rate_to_probability():
    if func_pitting_induction_time() != "No pitting":
        poisson_probability = poisson_constant / func_pitting_induction_time()
        return 0.0001 #poisson_probability
    else:
        return 0.0001

def func_image_array(height_px, width_px, rate, elasped_years):
    global image_array
    image_array = np.zeros((height_px, width_px), dtype=int)

    rows = image_array.shape[0]
    cols = image_array.shape[1]

    for month in range(int(elasped_years)):
        for i in range(0, rows):
            for j in range(0, cols):
                if random() < func_rate_to_probability():
                    image_array[i, j] += 1
                    print(image_array)


    plt.figure(num=2)
    plt.imshow(image_array, interpolation='none', cmap='Reds')

    global fig2
    fig2 = plt.figure(num=2)

    return fig2

def func_get_material_loss(height, width, rate, days):
    return

# DEFINE THE CANVAS SIZE FIRST! fig0 is useless apart from setting size
fig0 = plt.figure(num=0)
figure_x, figure_y, figure_w, figure_h = fig0.bbox.bounds

def draw_figure(canvas, figure, loc=(0, 0)):
    """ Draw a matplotlib figure onto a Tk canvas
    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py
    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = Tk.PhotoImage(master=canvas, width=figure_w, height=figure_h)
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)
    return photo

def plot_f(x, pitting_induction_time, pitting_corrosion_rate):
    if (x*24) <= float(pitting_induction_time):
        return 0
    else:
        depth = (x - pitting_induction_time/24) * (pitting_corrosion_rate / 365)
        print(depth)
        return depth


def annot_max(x,y, ax=None):
    xmax = x[np.argmax(y)]
    ymax = y.max()
    text= "Max pit depth of {:.3f} mm after {:.0f} days".format(ymax, xmax)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              bbox=bbox_props, ha="right", va="top", clip_on=True)
    ax.annotate(text, xy=(xmax, ymax), xytext=(0.94,0.96), **kw)


# END SIMULATION FUNCTIONS
##########################################

##########################################


##########################################

# BEGIN GUI CODE
##########################################

# TAB 1 - PITTING

left_col_pitting = [[sg.Text('Enter environmental parameters:')],
                    [sg.Text('pH', size=(15,1)), sg.Input('5', key='pitting_pH', size=(10,1))],
                    [sg.Text('[O2] (ppm)', size=(15,1)), sg.Input('10', key='pitting_O2', size=(10,1))],
                    [sg.Text('[Cl-] (ppm)', size=(15,1)), sg.Input('50', key='pitting_Clm', size=(10,1))],
                    [sg.Text('Temperature (°C)', size=(15,1)), sg.Input('28', key='pitting_tempC', size=(10,1))],
                    [sg.Text('Water velocity (m/s)', size=(15,1)), sg.Input('10', key='pitting_velocity', size=(10,1))],
                    [sg.Text('Bacteria biofilm', size=(15, 1)), sg.Drop(values=('Present', 'Absent'), key='pitting_biofilm', auto_size_text=True)],
                    [sg.Text('\nSteel alloy composition:')],
                    [sg.Text('% Chromium', size=(15,1)), sg.Input('13.115', key='pitting_Cr', size=(10,1))],
                    [sg.Text('% Molybdenum', size=(15,1)), sg.Input('0.029', key='pitting_Mo', size=(10,1))],
                    [sg.Text('% Nitrogen', size=(15,1)), sg.Input('0.019', key='pitting_N', size=(10,1))],
                    [sg.Text('% Manganese', size=(15,1)), sg.Input('0.786', key='pitting_Mn', size=(10,1))],
                    [sg.Text('% Sulfur', size=(15,1)), sg.Input('0.017', key='pitting_S', size=(10,1))],
                    [sg.Text('-'*77)],
                    [sg.Text('Pitting resistance equivalence number:', size=(28,1)), sg.Input('', size=(10,1), key='PREN_out')],
                    [sg.Text('Average pitting induction time:', size=(28,1)), sg.Input('', size=(10, 1), key='PIT_out')],
                    [sg.Text('Max pitting corrosion rate:', size=(28,1)), sg.Input('', size=(10, 1), key='pittingRate_out')],
                    [sg.Text('-' * 77)],
                    [sg.Button('Calculate'), sg.Button('User Guide and Governing Equations')]
                    ]

tab1a_maxpitdepth = [[sg.Canvas(size=(figure_w, figure_h), key='canvas_maxpitdepth')]]

tab1b_simulation = [[sg.Canvas(size=(figure_w, figure_h), key='canvas_simulation')]]

right_col_pitting = [[sg.Text('Surface height (mm)', auto_size_text=True), sg.Input('50', size=(5,1), key='surface_height'), sg.Text('×'),
                      sg.Text('Surface width (mm)', auto_size_text=True), sg.Input('50', size=(5,1), key='surface_width'),  sg.Text('for'),
                      sg.Input('500', size=(5,1), key='days_elasped'), sg.Text('days elasped: ', auto_size_text=True),
                      sg.Button('Simulate'), sg.Button('Clear')],
                     [sg.TabGroup([[sg.Tab(' Max pit depth against time ', tab1a_maxpitdepth),
                                    sg.Tab(' Stochastic simulation ', tab1b_simulation)]])]
                     ]

tab1_pitting = [[sg.Text('\nUse this to calculate pitting corrosion rate and mean time to pit formation\n')],
                [sg.Column(left_col_pitting), sg.Column(right_col_pitting)],
                ]

# TAB 2 - UNIFORM FRESHWATER

tab2_uniform_freshwater = [[sg.Txt('Use this to calculate uniform corrosion rate of carbon steel in freshwater')]]

tab3_uniform_seawater = [[sg.Txt('Enter parameter')]]

tab4_srb = [[sg.Txt('This is inside tab 2')],
            [sg.In(key='in1')]]

tab5_galvanic = [[sg.Txt('This is inside tab 2')],
                 [sg.In(key='in2')]]

layout = [[sg.TabGroup([[sg.Tab(' General Pitting Corrosion ', tab1_pitting),
                         sg.Tab(' Uniform Corrosion in Freshwater ', tab2_uniform_freshwater),
                         sg.Tab(' Uniform Corrosion in Seawater ', tab3_uniform_seawater),
                         sg.Tab(' Corrosion by Sulphate Reducing Bacteria ', tab4_srb),
                         sg.Tab(' Galvanic Corrosion ', tab5_galvanic)]])],
          ]

window = sg.Window('Corrosion Toolbox - PUB', layout, font=("Helvetica", 10)).Finalize()

while True:  # Event Loop
    event, values = window.Read()
    print(event, values)

    if event is None or event == 'Exit':
        break

    if event == 'Calculate':
        try:
            pH = float(values['pitting_pH'])
            conc_O2 = float(values['pitting_O2'])
            conc_Clm = float(values['pitting_Clm'])
            temp_c = float(values['pitting_tempC'])
            fluid_velocity = float(values['pitting_velocity'])

            if values['pitting_biofilm'] == 'Present':
                biofilm = True
            else:
                biofilm = False

            pct_Cr = float(values['pitting_Cr'])
            pct_Mo = float(values['pitting_Mo'])
            pct_N = float(values['pitting_N'])
            pct_Mn = float(values['pitting_Mn'])
            pct_S = float(values['pitting_S'])

        except:
            sg.Popup('One of the parameters is invalid')

        print(pH, conc_O2, conc_Clm, temp_c, fluid_velocity, biofilm, pct_Cr, pct_Mo, pct_N, pct_Mn, pct_S)

        pitting_potential, pren = func_pitting_potential()
        pitting_induction_time = func_pitting_induction_time()
        pitting_corrosion_rate = func_pitting_corrosion_rate()

        x = np.linspace(0, int(values['days_elasped']), num=5000)
        plot_f2 = np.vectorize(plot_f, otypes=[np.float])
        y = plot_f2(x, pitting_induction_time, pitting_corrosion_rate)

        # plot with various axes scales
        fig1 = plt.figure(num=1)

        # linear
        plt.plot(x, y)
        plt.title('Max pit depth against elasped days')
        plt.grid(True)
        annot_max(x, y)

        window.Element('PREN_out').Update(float(pren))
        window.Element('PIT_out').Update(float(pitting_induction_time))
        window.Element('pittingRate_out').Update(float(pitting_corrosion_rate))

        fig1_photo = draw_figure(window.FindElement('canvas_maxpitdepth').TKCanvas, fig1)

        calculated = True

    if event == 'Simulate' and calculated == True:
        fig2_photo = draw_figure(window.FindElement('canvas_simulation').TKCanvas,
                                 func_image_array(int(values['surface_height']),
                                                  int(values['surface_width']),50,
                                                  int(values['days_elasped'])))

window.Close()
