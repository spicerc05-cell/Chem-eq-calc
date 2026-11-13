
from dataclasses import dataclass, field
import re


# Molecular weights of elements (g/mol)
molecular_dictionary = {
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "Ne": 20.1797,
    "Na": 22.989768,
    "Mg": 24.3050,
    "Al": 26.981539,
    "Si": 28.0855,
    "Cl": 35.4527,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955910,
    "Ti": 47.88,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.93805,
    "Fe": 55.847,
    "Co": 58.93320,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.39,
    "Ga": 69.723,
    "Ge": 72.61,
    "As": 74.92159,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.80,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.94,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.82,
    "Sn": 118.710,
    "Sb": 121.757,
    "Te": 127.60,
    "Xe": 131.29,
    "Cs": 132.90543,
    "Ba": 137.327,
    "La": 138.9055,
    "Ce": 140.115,
    "Pr": 140.90765,
    "Nd": 144.24,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.965,
    "Gd": 157.25,
    "Tb": 158.92534,
    "Dy": 162.50,
    "Ho": 164.93032,
    "Er": 167.26,
    "Tm": 168.93421,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W": 183.85,
    "Re": 186.207,
    "Os": 190.2,
    "Ir": 192.22,
    "Pt": 195.08,
    "Au": 196.96654,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98037,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226.0254,
    "Ac": 227,
    "Th": 232.0381,
    "Pa": 232.0359,
    "U": 238.0289,
    "Np": 237.0482,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 260,
    "Rf": 261,
    "Db": 262,
    "Sg": 263,
    "Bh": 262,
    "Hs": 265,
    "Mt": 266,
    "Ds": 281,
    "Rg": 282,
    "Cn": 285,
    "Nh": 286,
    "Fl": 289,
    "Mc": 290,
    "Lv": 293,
    "Ts": 294,
    "Og": 294, 
    "N": 14.00674,
    "B": 10.811,
    "C": 12.011,
    "O": 15.9994,
    "F": 18.9984032,
    "I": 126.90447,
    "S": 32.066,
    "P": 30.973762,
    "H": 1.00794,
}


# Convert a chemical formula into element counts
def parse_formula(formula):
    pattern = r"([A-Z][a-z]?)(\d*)"
    matches = re.findall(pattern, formula)

    rebuilt_formula = "".join(
        element + number for element, number in matches
    )

    if rebuilt_formula != formula or not matches:
        raise ValueError("Invalid chemical formula.")

    atom_counts = {}

    for element, number in matches:
        if element not in molecular_dictionary:
            raise ValueError("Unknown element: " + element)

        count = int(number) if number else 1

        if count <= 0:
            raise ValueError("Subscripts must be positive.")

        atom_counts[element] = atom_counts.get(element, 0) + count

    return atom_counts


# Calculate molecular weight
def calculate_molar_mass(formula):
    atom_counts = parse_formula(formula)
    total_mass = 0

    for element, count in atom_counts.items():
        total_mass += molecular_dictionary[element] * count

    return total_mass


# Parent class for chemical compounds
@dataclass
class Compound:
    formula: str
    stoichiometric: int
    molecular_weight: float = field(init=False)

    def __post_init__(self):
        if self.stoichiometric <= 0:
            raise ValueError("Coefficient must be positive.")

        self.mw_calc()

    def mw_calc(self):
        self.molecular_weight = calculate_molar_mass(self.formula)


# Reactant class
@dataclass
class Reactant(Compound):
    weight: float

    def __post_init__(self):
        super().__post_init__()

        if self.weight <= 0:
            raise ValueError("Reactant mass must be positive.")

    @property
    def moles(self):
        return self.weight / self.molecular_weight


# Product class
@dataclass
class Product(Compound):
    pass


# Count all atoms on one side of the equation
def get_total_atoms(compounds):
    total_atoms = {}

    for compound in compounds:
        atoms = parse_formula(compound.formula)

        for element, count in atoms.items():
            total_atoms[element] = (
                total_atoms.get(element, 0)
                + count * compound.stoichiometric
            )

    return total_atoms


# Check whether the equation is balanced
def is_balanced(reactants, products):
    reactant_atoms = get_total_atoms(reactants)
    product_atoms = get_total_atoms(products)

    return reactant_atoms == product_atoms


def find_limiting_reactant(reactants):
    limiting_reactant = reactants[0]

    smallest_ratio = (
        limiting_reactant.moles
        / limiting_reactant.stoichiometric
    )

    for reactant in reactants[1:]:
        ratio = reactant.moles / reactant.stoichiometric

        if ratio < smallest_ratio:
            smallest_ratio = ratio
            limiting_reactant = reactant

    return limiting_reactant


def calculate_theoretical_yield(limiting_reactant, product):

    product_moles = (
        limiting_reactant.moles
        * product.stoichiometric
        / limiting_reactant.stoichiometric
    )

    theoretical_yield = product_moles * product.molecular_weight

    return theoretical_yield


