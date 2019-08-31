import math
import PySimpleGUI as sg

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

# BEGIN PITTING FUNCTIONS
############################################


def func_cathodic_potential():
    # Where conc_O2 is the oxygen concentration in ppm, conc_Cl is the CHLORIDE concentration in ppm
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
        k = 0

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
        pit = math.pow(10, ((pren / (2 * math.log10(conc_Clm))) * (1 - driving_voltage)))
        # pit is the pitting induction time in hours
        return pit, driving_voltage
    else:
        return "no pitting"


def func_pitting_corrosion_rate():
    # Returns pitting corrosion rate in mm/year
    oxygen_limiting_current = 10 * 2**((temp_c - 25)/25) * conc_O2 * (1 + math.sqrt(fluid_velocity))
    pitting_corrosion_rate = 20 * math.sqrt((oxygen_limiting_current / 1000) / (water_resistivity + 0.8))

    return pitting_corrosion_rate

print(func_pitting_induction_time(), func_pitting_corrosion_rate())

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


# BEGIN GUI CODE
##########################################


tab1_uniform_layout = [[sg.T('Enter parameter')]]

tab2_pitting_layout = [[sg.T('This is inside tab 2')],
                       [sg.In(key='in')]]

tab3_srb_layout = [[sg.T('This is inside tab 2')],
                   [sg.In(key='in')]]

layout = [[sg.TabGroup([[sg.Tab('Uniform', tab1_uniform_layout, tooltip='tip'), sg.Tab('Pitting', tab2_pitting_layout)]], tooltip='TIP2')],
          [sg.Button('Read')]]

window = sg.Window('Corrosion Toolbox - PUB', layout, font=("Helvetica", 10))

while True:
    event, values = window.Read()
    print(event,values)
    if event is None:           # always,  always give a way out!
        break
