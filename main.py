from bakery import assert_equal
from dataclasses import dataclass
from drafter import *
from dataclasses import dataclass, field

from Chemistry import (
    Reactant,
    Product,
    is_balanced,
    find_limiting_reactant,
    calculate_theoretical_yield,
    get_total_atoms
)

hide_debug_information()
set_website_title("Calculations for balanced Chem equations")
set_site_information(
    author="Connor Spicer",
    description="Calculates limiting reactant, actual yield, and theoretical yield from balanced equation",
    sources="Drafter docs"
)
set_website_framed(False)


# Stores one completed reaction calculation
@dataclass
class ReactionRecord:
    equation: str
    balanced: bool
    reactant_atoms: str
    product_atoms: str
    limiting_reactant: str
    theoretical_yield: float
    actual_yield: float = 0.0
    percent_yield: float = 0.0
    has_actual_yield: bool = False


# Stores information for the whole website
@dataclass
class State:
    reactant_1_formula: str = "H2"
    reactant_2_formula: str = "O2"
    product_formula: str = "H2O"

    reactant_1_coefficient: str = "2"
    reactant_2_coefficient: str = "1"
    product_coefficient: str = "2"

    reactant_1_mass: str = "4"
    reactant_2_mass: str = "32"

    error_message: str = ""
    result_message: str = ""

    actual_yield_input: str = ""

    calculations: list = field(default_factory=list)


# Converts atom counts into readable text
def format_atom_counts(atom_counts):
    parts = []

    for element in sorted(atom_counts):
        parts.append(element + ": " + str(atom_counts[element]))

    return ", ".join(parts)


# HOME PAGE
@route
def index(state: State) -> Page:
    return Page(state, [
        "CHEMICAL REACTION CALCULATOR\n",
        "Welcome! This application helps you analyze chemical reactions.\n",
        "You can check whether an equation is balanced, find the limiting reactant, and calculate theoretical and percent yield.\n",
        Button("Start Calculating", "calculator"),
        Button("View Examples", "examples"),
        Button("Help", "help_page"),
        Button("About", "about")
    ])


# CALCULATOR PAGE
@route
def calculator(state: State) -> Page:
    return Page(state, [
        "CHEMICAL REACTION CALCULATOR",
        LineBreak(),
        "Enter the formulas and coefficients for your reaction.",
        LineBreak(),
        LineBreak(),

        "REACTANT 1",
        LineBreak(),
        "Formula (example: H2):",
        LineBreak(),
        TextBox("reactant_1_formula", state.reactant_1_formula),
        LineBreak(),
        "Coefficient:",
        LineBreak(),
        TextBox("reactant_1_coefficient",
                state.reactant_1_coefficient, "number"),
        LineBreak(),
        "Starting mass in grams:",
        LineBreak(),
        TextBox("reactant_1_mass", state.reactant_1_mass, "number"),
        LineBreak(),
        LineBreak(),

        "REACTANT 2",
        LineBreak(),
        "Formula (example: O2):",
        LineBreak(),
        TextBox("reactant_2_formula", state.reactant_2_formula),
        LineBreak(),
        "Coefficient:",
        LineBreak(),
        TextBox("reactant_2_coefficient",
                state.reactant_2_coefficient, "number"),
        LineBreak(),
        "Starting mass in grams:",
        LineBreak(),
        TextBox("reactant_2_mass", state.reactant_2_mass, "number"),
        LineBreak(),
        LineBreak(),

        "PRODUCT",
        LineBreak(),
        "Formula (example: H2O):",
        LineBreak(),
        TextBox("product_formula", state.product_formula),
        LineBreak(),
        "Coefficient:",
        LineBreak(),
        TextBox("product_coefficient",
                state.product_coefficient, "number"),

        LineBreak(),
        LineBreak(),
        Button("Check Equation and Calculate", "calculate"),
        Button("Home", "index")
    ])

# CALCULATE ROUTE
@route
def calculate(
    state: State,
    reactant_1_formula: str,
    reactant_1_coefficient: str,
    reactant_1_mass: str,
    reactant_2_formula: str,
    reactant_2_coefficient: str,
    reactant_2_mass: str,
    product_formula: str,
    product_coefficient: str
) -> Page:

    # Save the user's inputs
    state.reactant_1_formula = reactant_1_formula.strip()
    state.reactant_2_formula = reactant_2_formula.strip()
    state.product_formula = product_formula.strip()

    state.reactant_1_coefficient = reactant_1_coefficient
    state.reactant_2_coefficient = reactant_2_coefficient
    state.product_coefficient = product_coefficient

    state.reactant_1_mass = reactant_1_mass
    state.reactant_2_mass = reactant_2_mass

    state.error_message = ""
    state.result_message = ""

    try:
        # Check that no fields are empty
        if (not reactant_1_formula.strip()
                or not reactant_2_formula.strip()
                or not product_formula.strip()):
            raise ValueError("Please enter all chemical formulas.")

        # Convert coefficients and masses into numbers
        r1_coefficient = int(reactant_1_coefficient)
        r2_coefficient = int(reactant_2_coefficient)
        p_coefficient = int(product_coefficient)

        r1_mass = float(reactant_1_mass)
        r2_mass = float(reactant_2_mass)

        # Create chemistry objects
        reactants = [
            Reactant(reactant_1_formula.strip(),
                     r1_coefficient, r1_mass),
            Reactant(reactant_2_formula.strip(),
                     r2_coefficient, r2_mass)
        ]

        products = [
            Product(product_formula.strip(), p_coefficient)
        ]

        # Check whether the equation is balanced
        balanced = is_balanced(reactants, products)

        if not balanced:
            state.error_message = (
                "Your equation is NOT balanced. "
                "Please correct the coefficients and try again."
            )
            return results(state)

        # Find the limiting reactant
        limiting = find_limiting_reactant(reactants)

        # Calculate theoretical yield
        theoretical_yield = calculate_theoretical_yield(
            limiting, products[0]
        )

        # Get atom counts for the details page
        reactant_atoms = format_atom_counts(
            get_total_atoms(reactants)
        )

        product_atoms = format_atom_counts(
            get_total_atoms(products)
        )

        # Create a record of this calculation
        equation = (
            str(r1_coefficient) + reactant_1_formula.strip()
            + " + "
            + str(r2_coefficient) + reactant_2_formula.strip()
            + " -> "
            + str(p_coefficient) + product_formula.strip()
        )

        record = ReactionRecord(
            equation=equation,
            balanced=True,
            reactant_atoms=reactant_atoms,
            product_atoms=product_atoms,
            limiting_reactant=limiting.formula,
            theoretical_yield=theoretical_yield
        )

        state.calculations.append(record)

        state.result_message = "Your equation is balanced!"

        return results(state)

    except ValueError as error:
        state.error_message = str(error)
        return results(state)


# RESULTS PAGE
@route
def results(state: State) -> Page:

    content = [
        "CALCULATION RESULTS\n"
    ]

    if state.error_message:
        content.append("ERROR: " + state.error_message + "\n")

    if state.calculations:
        latest = state.calculations[-1]

        content.extend([
            state.result_message + "\n",
            "Reaction: " + latest.equation + "\n",
            "Limiting reactant: "
            + latest.limiting_reactant + "\n",
            "Theoretical yield: "
            + str(round(latest.theoretical_yield, 2))
            + " grams\n"
        ])

        if latest.has_actual_yield:
            content.extend([
                "Actual yield: "
                + str(round(latest.actual_yield, 2))
                + " grams\n",
                "Percent yield: "
                + str(round(latest.percent_yield, 2))
                + "%\n"
            ])
        else:
            content.extend([
                "\nEnter your experimentally measured actual yield:\n",
                TextBox("actual_yield", state.actual_yield_input, "number"),
                Button("Calculate Percent Yield", "calculate_percent_yield")
            ])

    # Display all saved calculation records
    content.append("\nYOUR CALCULATIONS\n")

    for record in state.calculations:
        content.append(
            record.equation
            + " | Limiting reactant: "
            + record.limiting_reactant
            + " | Theoretical yield: "
            + str(round(record.theoretical_yield, 2))
            + " g\n"
        )

    content.extend([
        Button("Reaction Details", "reaction_details"),
        Button("Calculate Another Reaction", "calculator"),
        Button("Home", "index")
    ])

    return Page(state, content)


# CALCULATE ACTUAL AND PERCENT YIELD
@route
def calculate_percent_yield(
    state: State,
    actual_yield: str
) -> Page:

    if not state.calculations:
        state.error_message = "Calculate a reaction first."
        return results(state)

    try:
        actual = float(actual_yield)

        if actual < 0:
            raise ValueError("Actual yield cannot be negative.")

        latest = state.calculations[-1]

        if latest.theoretical_yield <= 0:
            raise ValueError("Theoretical yield must be positive.")

        latest.actual_yield = actual

        latest.percent_yield = (
            actual / latest.theoretical_yield
        ) * 100

        latest.has_actual_yield = True

        state.actual_yield_input = str(actual)
        state.error_message = ""
        state.result_message = "Percent yield calculated!"

    except ValueError as error:
        state.error_message = str(error)

    return results(state)


# REACTION DETAILS PAGE
@route
def reaction_details(state: State) -> Page:

    if not state.calculations:
        return Page(state, [
            "No reaction has been calculated yet.\n",
            Button("Go to Calculator", "calculator")
        ])

    latest = state.calculations[-1]

    return Page(state, [
        "REACTION DETAILS\n",
        "Equation: " + latest.equation + "\n",
        "Balanced: " + str(latest.balanced) + "\n",

        "\nTotal atoms on the reactant side:\n",
        latest.reactant_atoms + "\n",

        "\nTotal atoms on the product side:\n",
        latest.product_atoms + "\n",

        "\nLimiting reactant: "
        + latest.limiting_reactant + "\n",

        "Theoretical yield: "
        + str(round(latest.theoretical_yield, 2))
        + " grams\n",

        Button("Back to Results", "results"),
        Button("Home", "index")
    ])


# EXAMPLES PAGE
@route
def examples(state: State) -> Page:
    return Page(state, [
        "CHEMICAL REACTION EXAMPLES\n",

        "Example 1: Formation of water\n",
        "2H2 + O2 -> 2H2O\n",

        "Reactant 1: H2\n",
        "Reactant 2: O2\n",
        "Product: H2O\n",

        "\nExample 2: Formation of carbon dioxide\n",
        "C + O2 -> CO2\n",

        "Reactant 1: C\n",
        "Reactant 2: O2\n",
        "Product: CO2\n",

        Button("Try the Calculator", "calculator"),
        Button("Home", "index")
    ])


# HELP PAGE
@route
def help_page(state: State) -> Page:
    return Page(state, [
        "HELP\n",

        "1. Enter the chemical formulas for your reactants and product.\n",
        "2. Enter the coefficient for each compound.\n",
        "3. Enter the starting mass of each reactant in grams.\n",
        "4. Click Check Equation and Calculate.\n",
        "5. If the equation is balanced, the calculator finds the limiting reactant and theoretical yield.\n",
        "6. Enter your experimentally measured actual yield.\n",
        "7. Click Calculate Percent Yield.\n",
        "8. EQUATIONS MUST BE BALANCED FIRST.\n",

        "\nRemember: coefficients are the numbers in front of chemical formulas.\n",
        "For example, in 2H2 + O2 -> 2H2O, the coefficients are 2, 1, and 2.\n",

        Button("Calculator", "calculator"),
        Button("Home", "index")
    ])


# ABOUT PAGE
@route
def about(state: State) -> Page:
    return Page(state, [
        "ABOUT THIS PROJECT\n",

        "The Chemical Reaction Calculator checks whether chemical equations are balanced, identifies limiting reactants, calculates theoretical yield, and calculates percent yield.\n",
        

        "The calculator is designed for simple chemical formulas and reactions with two reactants and one product.\n",

        Button("Home", "index"),
        Button("Calculator", "calculator")
    ])


# START THE WEBSITE
start_server(State())

