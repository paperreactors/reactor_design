
import numpy as np
import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

neutron_energy = Q_(0.0253, 'eV')
element_data = {} #create a dictionary of elements

### All element data sourced from Lamarsh, Introduction to Nuclear Engineering
### Cross sections at 0.0253 eV (2200 m/s)
O = {}
O['A'] = Q_(15.9994,'g/mol') #atomic weight
O['sig_a'] = Q_(0.00027,'barn') #microscopic XS
O['sig_s'] = Q_(3.76,'barn') #microscopic XS
element_data['O'] = O #add each element to the dictionary

C = {}
C['A'] = Q_(12.01115,'g/mol')
C['sig_a'] = Q_(0.0034,'barn')
C['sig_s'] = Q_(4.75,'barn')
element_data['C'] = C

H = {}
H['A'] = Q_(1.00797,'g/mol')
H['sig_a'] = Q_(0.332,'barn')
H['sig_s'] = Q_(49.62,'barn')
element_data['H'] = H

D = {}
D['A'] = Q_(2.01410,'g/mol')
D['sig_a'] = Q_(0.00053,'barn')
D['sig_s'] = Q_(4.92,'barn')
element_data['D'] = D

class moderator:
    def __init__(self,element_data,element_dict,rho):
        self.element_dict = element_dict
        self.rho = rho
        Ao = Q_(6.02214E23,'1/mol')

        # Calculate ALL the things!
        self.calc_A()
        self.calc_sig_a()
        self.calc_sig_s()
        self.calc_SIG_a(Ao)
        self.calc_SIG_s(Ao)
        self.calc_ksi()
        self.calc_MR()

        print self.element_dict
        print self.SIG_a
        print self.SIG_s
        print self.ksi
        print self.MR


    def calc_A(self):
        A = 0
        for key,value in self.element_dict.iteritems():
            # sum weight of constituent elements
            A = A + element_data[key]['A'] * value
        self.A = A

    def calc_sig_a(self):
        sig_a = 0
        for key,value in self.element_dict.iteritems():
            # sum microscopic cross sections
            sig_a = sig_a + element_data[key]['sig_a']*value
        self.sig_a = sig_a

    def calc_sig_s(self):
        sig_s = 0
        for key,value in self.element_dict.iteritems():
            # sum microscopic cross sections
            sig_s = sig_s + element_data[key]['sig_s']*value
        self.sig_s = sig_s

    def calc_SIG_a(self,Ao):
        # calculate macroscopic cross section (atom density * microscopic XS)
        SIG_a = (self.rho * Ao / self.A) * self.sig_a
        SIG_a.ito(1/ureg.cm)
        self.SIG_a = SIG_a

    def calc_SIG_s(self,Ao):
        SIG_s = (self.rho * Ao / self.A) * self.sig_s
        SIG_s.ito(1/ureg.cm)
        self.SIG_s = SIG_s

    def calc_ksi(self):
        # Calculate avg logarithmic energy decrement
        # calculace ksi for each element, and weight them by the scattering XS
        # Thanks for this explanation Carl! https://www.quora.com/What-is-the-average-fractional-energy-loss-per-collision-of-a-neutron-of-high-density-polyethylene
        numerator = 0
        for key,value in self.element_dict.iteritems():
            A = element_data[key]['A'].magnitude
            ksi_element = 1 + np.log((A-1.)/(A+1.))*(A-1.)**2.0/(2.0*A)
            numerator = numerator + element_data[key]['sig_s']*value*ksi_element

        ksi = numerator / self.sig_s
        self.ksi = ksi

    def calc_MR(self):
        # Finally, calculate moderating ratio!
        MR = self.ksi * self.SIG_s / self.SIG_a
        self.MR = MR


# Some simple materials:
H2O = moderator(element_data,{'H':2,'O':1},Q_(1.0,'g/cc'))
D2O = moderator(element_data,{'D':2,'O':1},Q_(1.11,'g/cc'))
syrup = moderator(element_data,{'C':12,'H':22,'O':11},Q_(1.59,'g/cc'))
heavysyrup = moderator(element_data,{'C':12,'D':22,'O':11},Q_(1.59,'g/cc'))
graphite = moderator(element_data,{'C':1},Q_(1.6,'g/cc'))
biphenyl = moderator(element_data,{'C':12,'H':10},Q_(1.04,'g/cc'))
terphenyl = moderator(element_data,{'C':18,'H':14},Q_(1.24,'g/cc'))
# these last two are typical in organic moderated/cooled reactors!
# check out https://en.wikipedia.org/wiki/Piqua_Nuclear_Generating_Station
# shout out to Piqua, OH
