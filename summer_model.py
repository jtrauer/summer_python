
import numpy


class EpiModel:
    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=None, output_connections=(),
                 tracked_quantities=()):

        # features not to be changed
        self.available_birth_approaches = ["add_crude_birth_rate", "replace_deaths", "no_births"]

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, report, reporting_sigfigs, entry_compartment,
            starting_population, default_starting_compartment, equilibrium_stopping_tolerance, output_connections,
            tracked_quantities)

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        self.times, self.compartment_types, self.initial_conditions, self.parameters, self.requested_flows, \
            self.initial_conditions_to_total, self.infectious_compartment, self.birth_approach, self.report, \
            self.reporting_sigfigs, self.entry_compartment, self.starting_population, \
            self.default_starting_population, self.equilibrium_stopping_tolerance, self.output_connections,\
            self.tracked_quantities = [None] * 16

        # convert input arguments to model attributes
        for attribute in ["times", "compartment_types", "initial_conditions", "parameters", "infectious_compartment",
                          "birth_approach", "report", "reporting_sigfigs", "entry_compartment",
                          "starting_population", "default_starting_compartment", "infectious_compartment",
                          "equilibrium_stopping_tolerance", "output_connections", "tracked_quantities"]:
            setattr(self, attribute, eval(attribute))

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.add_default_quantities()

    def output_to_user(self, comment):
        """
        short function to save the if statement in every call to output some information
        """
        if self.report:
            print(comment)

    def check_and_report_attributes(
            self, times, compartment_types, initial_conditions, parameters, requested_flows,
            initial_conditions_to_total, infectious_compartment, birth_approach, report, reporting_sigfigs,
            entry_compartment, starting_population, default_starting_compartment, equilibrium_stopping_tolerance,
            output_connections, tracked_quantities):
        """
        check all input data have been requested correctly
        """

        # check that variables are of the expected type
        for expected_numeric_variable in ["reporting_sigfigs", "starting_population"]:
            if not isinstance(eval(expected_numeric_variable), int):
                raise TypeError("expected integer for %s" % expected_numeric_variable)
        for expected_list in ["times", "compartment_types", "requested_flows"]:
            if not isinstance(eval(expected_list), list):
                raise TypeError("expected list for %s" % expected_list)
        for expected_string in \
                ["infectious_compartment", "birth_approach", "entry_compartment", "default_starting_compartment"]:
            if not isinstance(eval(expected_string), str):
                raise TypeError("expected string for %s" % expected_string)
        for expected_boolean in ["initial_conditions_to_total", "report"]:
            if not isinstance(eval(expected_boolean), bool):
                raise TypeError("expected boolean for %s" % expected_boolean)
        for expected_dict in ["output_connections", "tracked_quantities"]:
            if not isinstance(eval(expected_dict), tuple):
                raise TypeError("expected dictionary for %s" % expected_dict)

        # check some specific requirements
        if infectious_compartment not in compartment_types:
            ValueError("infectious compartment name is not one of the listed compartment types")
        if birth_approach not in self.available_birth_approaches:
            ValueError("requested birth approach unavailable")
        if sorted(times) != times:
            self.output_to_user("requested integration times are not sorted, now sorting")
            self.times = sorted(self.times)

    def set_initial_conditions(self, initial_conditions_to_total):
        """
        set starting compartment values
        """

        pass

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain value
        """

        pass

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total
        """

        pass

    def implement_flows(self, requested_flows):
        """
        add all flows to create data frames from input lists
        """

        pass

    def add_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed
        """

        pass


if __name__ == "__main__":
    sir_model = EpiModel(numpy.linspace(0, 60 / 365, 60).tolist(),
                         ["susceptible", "infectious", "recovered"],
                         {"infectious": 0.001},
                         {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
                         [{"type": "standard_flows", "parameter": "recovery", "from": "infectious", "to": "recovered"},
                          {"type": "infection_density", "parameter": "beta", "from": "susceptible", "to": "infectious"},
                          {"type": "compartment_death", "parameter": "infect_death", "from": "infectious"}])
