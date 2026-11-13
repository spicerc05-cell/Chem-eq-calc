from bakery import assert_equal
from drafter import *
from dataclasses import dataclass



hide_debug_information()
set_website_title("Calculations for balanced Chem equations")
set_site_information(
    author="Connor Spicer",
    description="Calculates limiting reactant, actual yield, and theoretical yield from balanced equation"
    sources="Drafter docs"
)
set_website_framed(False)


@dataclass
class State:
    pass


@route
def index(state: State) -> Page:
    return Page(state, ["Hello World!"])

assert_equal(index(State()), Page(State(), ["Hello World!"]))

start_server(State())
