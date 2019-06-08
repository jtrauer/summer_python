
import numpy
from scipy.integrate import odeint
import matplotlib.pyplot


def find_stem(stratified_string):
    """
    find the stem of the compartment name as the text leading up to the first occurrence of "X"
    """
    first_x_location = stratified_string.find("X")
    if first_x_location == -1:
        return stratified_string
    else:
        return stratified_string[: first_x_location]


def increment_compartment(ode_equations, compartment_number, increment):
    """
    general method to increment the odes by a value specified as an argument
    """
    ode_equations[compartment_number] += increment
    return ode_equations


class EpiModel:
    """
    model construction methods
    """

    def find_parameter_value(self, parameter_name, time):
        """
        find the value of a parameter with time-variant values trumping constant ones
        """
        if parameter_name in self.time_variants:
            return self.time_variants[parameter_name](time)
        else:
            return self.parameters[parameter_name]

    def output_to_user(self, comment):
        """
        short function to save the if statement in every call to output some information, may be adapted later and was
        more important to the R version of the repository
        """
        if self.report:
            print(comment)

    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=None):

        # attributes with specific format that are independent of user inputs
        self.compartment_values, self.tracked_quantities, self.output_connections, self.time_variants = \
            [{} for _ in range(4)]
        self.derived_outputs = {"times": []}
        self.transition_flows, self.death_flows = [[] for _ in range(2)]

        # features that should not be changed
        self.available_birth_approaches = ["add_crude_birth_rate", "replace_deaths", "no_births"]

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, report, reporting_sigfigs, entry_compartment,
            starting_population, default_starting_compartment, equilibrium_stopping_tolerance)

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        self.times, self.compartment_types, self.initial_conditions, self.parameters, self.requested_flows, \
            self.initial_conditions_to_total, self.infectious_compartment, self.birth_approach, self.report, \
            self.reporting_sigfigs, self.entry_compartment, self.starting_population, \
            self.default_starting_compartment, self.default_starting_population, self.equilibrium_stopping_tolerance, \
            self.unstratified_flows, self.outputs = [None for _ in range(17)]

        # convert input arguments to model attributes
        for attribute in ["times", "compartment_types", "initial_conditions", "parameters",
                          "initial_conditions_to_total", "infectious_compartment", "birth_approach", "report",
                          "reporting_sigfigs", "entry_compartment", "starting_population",
                          "default_starting_compartment", "infectious_compartment", "equilibrium_stopping_tolerance"]:
            setattr(self, attribute, eval(attribute))

        # set initial conditions and implement flows
        self.set_initial_conditions(initial_conditions_to_total)

        # implement unstratified flows
        self.implement_flows(requested_flows)

        # add any missing quantities that will be needed
        self.add_default_quantities()

    def check_and_report_attributes(
            self, times, compartment_types, initial_conditions, parameters, requested_flows,
            initial_conditions_to_total, infectious_compartment, birth_approach, report, reporting_sigfigs,
            entry_compartment, starting_population, default_starting_compartment, equilibrium_stopping_tolerance):
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

        # check some specific requirements
        if infectious_compartment not in compartment_types:
            ValueError("infectious compartment name is not one of the listed compartment types")
        if birth_approach not in self.available_birth_approaches:
            ValueError("requested birth approach unavailable")
        if sorted(times) != times:
            self.output_to_user("requested integration times are not sorted, now sorting")
            self.times = sorted(self.times)

        # report on characteristics of inputs
        if report:
            print("integrating from time %s to %s"
                  % (round(times[0], reporting_sigfigs), round(times[-1], reporting_sigfigs)))
            print("unstratified requested initial conditions are:")
            for compartment in initial_conditions:
                print("\t%s: %s" % (compartment, initial_conditions[compartment]))
            print("infectious compartment is called '%s'" % infectious_compartment)
            print("birth approach is %s" % birth_approach)

    def set_initial_conditions(self, initial_conditions_to_total):
        """
        set starting compartment values
        """

        # set starting values of unstratified compartments to requested value, or zero if no value requested
        for compartment in self.compartment_types:
            if compartment in self.initial_conditions:
                self.compartment_values[compartment] = self.initial_conditions[compartment]
            else:
                self.output_to_user("no starting value requested for %s so set to zero" % compartment)
                self.compartment_values[compartment] = 0

        # sum to a total value if requested
        if initial_conditions_to_total:
            self.sum_initial_compartments_to_total()

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain value
        """
        compartment = self.find_remainder_compartment()
        if sum(self.compartment_values.values()) > self.starting_population:
            raise ValueError("total of requested compartment values is greater than the requested starting population")
        remaining_population = self.starting_population - sum(self.compartment_values.values())
        self.output_to_user("requested that total population sum to %s" % self.starting_population)
        self.output_to_user("remaining population of %s allocated to %s compartment"
                            % (remaining_population, compartment))
        self.compartment_values[compartment] = remaining_population

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total
        """
        if len(self.default_starting_compartment) > 0 and \
                self.default_starting_compartment not in self.compartment_values:
            raise ValueError("starting compartment to populate with initial values not found in available compartments")
        elif len(self.default_starting_compartment) > 0:
            return self.default_starting_compartment
        else:
            self.output_to_user("no default starting compartment requested for unallocated population, " +
                                "so will be allocated to entry compartment %s" % self.entry_compartment)
            return self.entry_compartment

    def implement_flows(self, requested_flows):
        """
        add all flows to create data frames from input lists
        """
        for flow in requested_flows:
            if flow["parameter"] not in self.parameters:
                raise ValueError("flow parameter not found in parameter list")
            if flow["from"] not in self.compartment_types:
                raise ValueError("from compartment name not found in compartment types")
            if "to" in flow and flow["to"] not in self.compartment_types:
                raise ValueError("to compartment name not found in compartment types")

            if flow["type"] == "compartment_death":
                self.add_death_flow(flow)
            else:
                self.add_transition_flow(flow)

            if "infection" in flow["type"]:
                self.tracked_quantities["infectious_population"] = 0
            if flow["type"] == "infection_frequency":
                self.tracked_quantities["total_population"] = 0

        # retain a copy of the original flows for the purposes of graphing, etc.
        self.unstratified_flows = self.transition_flows

    def add_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed
        """

        # universal death rate
        if "universal_death_rate" not in self.parameters:
            self.parameters["universal_death_rate"] = 0

        # birth approach-specific parameters
        if self.birth_approach == "add_crude_birth_rate" and "crude_birth_rate" not in self.parameters:
            self.parameters["crude_birth_rate"] = 0
        elif self.birth_approach == "replace_deaths":
            self.tracked_quantities["total_deaths"] = 0

        # for each derived output to be recorded, initialise a tracked quantities key to zero
        for output in self.output_connections:
            self.tracked_quantities[output] = 0

        # parameters essential for stratification
        self.parameters["entry_fractions"] = 1

    def add_transition_flow(self, flow):
        """
        simply add a flow to the data frame storing the flows
        """
        flow["implement"] = 0
        self.transition_flows.append(flow)

    def add_death_flow(self, flow):
        """
        similarly for compartment-specific death flows
        """
        flow["implement"] = 0
        self.death_flows.append(flow)

    """
    methods for model running
    """

    def run_model(self):
        """
        integrate model odes
        """
        self.output_to_user("now integrating")
        self.prepare_stratified_parameter_calculations()

        def make_model_function(compartment_values, time):
            self.update_tracked_quantities(compartment_values)
            return self.apply_all_flow_types_to_odes([0] * len(self.compartment_values), compartment_values, time)

        self.outputs = odeint(make_model_function, list(self.compartment_values.values()), self.times)
        self.output_to_user("integration complete")

    def prepare_stratified_parameter_calculations(self):
        """
        for use in the stratified version
        """
        pass

    def set_stopping_conditions(self):
        """
        this and the birth rate method are the only two outstanding methods to the base class still to develop
        """
        pass

    def apply_all_flow_types_to_odes(self, ode_equations, compartment_values, time):
        """
        apply all flow types to a vector of zeros (deaths must come before births in case births replace deaths)
        """
        ode_equations = self.apply_transition_flows(ode_equations, compartment_values, time)
        if len(self.death_flows) > 0:
            self.apply_compartment_death_flows(ode_equations, compartment_values, time)
        ode_equations = self.apply_universal_death_flow(ode_equations, compartment_values, time)
        ode_equations = self.apply_birth_rate(ode_equations, compartment_values)
        return ode_equations

    def apply_transition_flows(self, ode_equations, compartment_values, time):
        """
        add fixed or infection-related flow to odes
        """

        # find adjusted parameter value
        for flow in self.transition_flows:
            if flow["implement"] == 0:
                adjusted_parameter = self.get_parameter_value(flow["parameter"], time)

                # find from compartment and "infectious population", which is 1 for standard flows
                infectious_population = self.find_infectious_multiplier(flow["type"])

                # calculate the flow and apply to the odes
                from_compartment = list(self.compartment_values.keys()).index(flow["from"])
                net_flow = adjusted_parameter * compartment_values[from_compartment] * infectious_population
                ode_equations = increment_compartment(ode_equations, from_compartment, -net_flow)
                ode_equations = increment_compartment(
                    ode_equations, list(self.compartment_values.keys()).index(flow["to"]), net_flow)

                # track any quantities dependent on flow rates
                self.track_derived_outputs(flow, net_flow)

        # add another element to the derived outputs vector
        self.extend_derived_outputs(time)

        # return flow rates
        return ode_equations

    def track_derived_outputs(self, flow, net_flow):
        """
        calculate derived quantities to be tracked
        """
        for output_type in self.output_connections:
            if self.output_connections[output_type]["from"] in flow["from"] \
                    and self.output_connections[output_type]["to"] in flow["to"]:
                self.tracked_quantities[output_type] += net_flow

    def extend_derived_outputs(self, time):
        """
        add the derived quantities being tracked to the end of the tracking vector
        """
        self.derived_outputs["times"].append(time)
        for output_type in self.output_connections:
            self.derived_outputs[output_type].append(self.tracked_quantities[output_type])

    def apply_compartment_death_flows(self, ode_equations, compartment_values, time):
        """
        equivalent method to for transition flows above, but for deaths
        """
        for flow in self.death_flows:
            if flow["implement"] == 0:
                adjusted_parameter = self.get_parameter_value(flow["parameter"], time)
                from_compartment = list(self.compartment_values.keys()).index(flow["from"])
                net_flow = adjusted_parameter * compartment_values[from_compartment]
                ode_equations = increment_compartment(ode_equations, from_compartment, -net_flow)
                if "total_deaths" in self.tracked_quantities:
                    self.tracked_quantities["total_deaths"] += net_flow
        return ode_equations

    def apply_universal_death_flow(self, ode_equations, compartment_values, time):
        """
        apply the population-wide death rate to all compartments
        """
        for compartment in self.compartment_values:
            adjusted_parameter = self.get_parameter_value("universal_death_rate", time)
            from_compartment = list(self.compartment_values.keys()).index(compartment)
            net_flow = adjusted_parameter * compartment_values[from_compartment]
            ode_equations = increment_compartment(ode_equations, from_compartment, -net_flow)

            # track deaths in case births are meant to replace deaths
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return ode_equations

    def apply_birth_rate(self, ode_equations, compartment_values):
        """
        apply a population-wide death rate to all compartments
        """
        return increment_compartment(
            ode_equations, list(self.compartment_values.keys()).index(self.entry_compartment),
            self.find_total_births(compartment_values))

    def find_total_births(self, compartment_values):
        """
        work out the total births to apply dependent on the approach requested
        """
        if self.birth_approach == "add_crude_birth_rate":
            return self.parameters["crude_birth_rate"] * sum(compartment_values)
        elif self.birth_approach == "replace_deaths":
            return self.tracked_quantities["total_deaths"]
        else:
            return 0.0

    def find_infectious_multiplier(self, flow_type):
        """
        find the multiplier to account for the infectious population in dynamic flows
        """
        if flow_type == "infection_density":
            return self.tracked_quantities["infectious_population"]
        elif flow_type == "infection_frequency":
            return self.tracked_quantities["infectious_population"] / \
                   self.tracked_quantities["total_population"]
        else:
            return 1

    def update_tracked_quantities(self, compartment_values):
        """
        update quantities that emerge during model running (not pre-defined functions of time)
        """
        for quantity in self.tracked_quantities:
            self.tracked_quantities[quantity] = 0
            if quantity == "infectious_population":
                self.find_infectious_population(compartment_values)
            elif quantity == "total_population":
                self.tracked_quantities["total_population"] = sum(compartment_values)

    def find_infectious_population(self, compartment_values):
        """
        calculations to find the effective infectious population
        """
        for compartment in self.compartment_values:
            if find_stem(compartment) == self.infectious_compartment:
                self.tracked_quantities["infectious_population"] += \
                    compartment_values[list(self.compartment_values.keys()).index(self.infectious_compartment)]

    def get_parameter_value(self, parameter, time):
        """
        need to split this out as a function in order to allow stratification later
        """
        return self.find_parameter_value(parameter, time)


class StratifiedModel(EpiModel):
    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=None):
        EpiModel.__init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=None)

        self.strata, self.removed_compartments, self.overwrite_parameter, self.compartment_types_to_stratify, \
            self.parameter_components = [[]] * 5
        self.heterogeneous_infectiousness = False

    """
    pre-integration methods
    """

    def stratify(self, stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests=(),
                 requested_proportions=(), infectiousness_adjustments=(), report=True):
        """
        initial preparation and checks
        """
        strata_names = self.prepare_and_check_stratification(
            stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, report)

        # stratify the compartments
        requested_proportions = self.tidy_starting_proportions(strata_names, requested_proportions, report)
        self.stratify_compartments(
            stratification_name, strata_names, adjustment_requests, requested_proportions, report)

        # stratify the flows
        self.stratify_transition_flows(stratification_name, strata_names, adjustment_requests, report)
        self.stratify_entry_flows(stratification_name, strata_names, requested_proportions, report)
        if len(self.death_flows) > 0:
            self.stratify_death_flows(stratification_name, strata_names, adjustment_requests, report)
        self.stratify_universal_death_rate(stratification_name, strata_names, adjustment_requests, report)

        # heterogeneous infectiousness adjustments
        self.apply_heterogeneous_infectiousness(stratification_name, strata_request, infectiousness_adjustments)

        # work out ageing flows (comes first so that the compartment names are still in the unstratified form)
        if stratification_name == "age":
            self.set_ageing_rates(strata_names, report)

    def prepare_and_check_stratification(self, stratification_name, strata_request, compartment_types_to_stratify,
                                         adjustment_requests, report):
        """
        initial preparation and checks
        """
        self.report = report
        if stratification_name == "age":
            strata_request = self.check_age_stratification(strata_request, compartment_types_to_stratify)
        else:
            self.output_to_user("\nimplementing stratification for: %s" % stratification_name)

        # make sure all stratification names are characters
        if type(stratification_name) != str:
            stratification_name = str(stratification_name)
            self.output_to_user("converting stratification name %s to character" % stratification_name)

        # record stratification as model attribute, find the names to apply strata and check requests
        self.strata.append(stratification_name)
        strata_names = self.find_strata_names_from_input(strata_request)
        self.check_compartment_request(compartment_types_to_stratify)
        self.check_parameter_adjustment_requests(adjustment_requests, strata_names)
        return strata_names

    def check_age_stratification(self, strata_request, compartment_types_to_stratify):
        """
        check that request meets the requirements for stratification by age
        """
        self.output_to_user("implementing age-specific stratification with specific behaviour")
        if "age" in self.strata:
            raise ValueError("requested age stratification, but has specific behaviour and can only be applied once")
        elif len(compartment_types_to_stratify) > 0:
            raise ValueError("requested age stratification, but compartment request should be passed as empty vector " +
                             "in order to apply to all compartments")
        elif any([type(stratum) != int and type(stratum) != float for stratum in strata_request]):
            raise ValueError("inputs for age strata breakpoints are not numeric")
        if 0 not in strata_request:
            self.output_to_user("adding age stratum called '0' as not requested, to represent those aged less than %s"
                                % min(strata_request))
            strata_request.append(0)
        if strata_request != sorted(strata_request):
            strata_request = sorted(strata_request)
            self.output_to_user("requested age strata not ordered, so have been sorted to: %s" % strata_request)
        return strata_request

    def find_strata_names_from_input(self, strata_request):
        """
        find the names of the stratifications from a particular user request
        """
        if type(strata_request) == list and (strata_request) == 0:
            raise ValueError("requested to stratify, but no strata provided")
        elif type(strata_request) == int:
            strata_names = numpy.arange(1, strata_request + 1)
            self.output_to_user("integer provided as strata labels for stratification, hence strata implemented " +
                                "are integers from 1 to %s" % strata_request)
        elif type(strata_request) != str:
            raise ValueError("number passed as request for strata labels, but not an integer greater than one, " +
                             "unclear what to do, stratification failed")
        else:
            strata_names = strata_request
        for name in strata_names:
            self.output_to_user("adding stratum: %s" % name)
        return strata_names

    def check_compartment_request(self, compartment_types_to_stratify):
        """
        check the requested compartments to be stratified has been requested correctly
        """

        # if vector of length zero passed, stratify all the compartment types in the model
        if len(compartment_types_to_stratify) == 0:
            self.output_to_user("no compartment names specified for this stratification, " +
                                "so stratification applied to all model compartments")
            self.compartment_types_to_stratify = self.compartment_types

        # otherwise check all the requested compartments are available and implement the user request
        elif any([compartment not in self.compartment_types for compartment in self.compartment_types_to_stratify]):
            raise ValueError("requested compartment or compartments to be stratified are not available in this model")
        else:
            self.compartment_types_to_stratify = compartment_types_to_stratify

    def check_parameter_adjustment_requests(self, adjustment_requests, strata_names):
        """
        check parameter adjustments have been requested appropriately
        """
        pass

    def tidy_starting_proportions(self, strata_names, requested_proportions, report):
        """
        prepare user inputs for starting proportions as needed
        """
        return 0

    def stratify_compartments(self,
                              stratification_name, strata_names, adjustment_requests, requested_proportions, report):
        """
        compartment stratification
        """
        pass

    def stratify_transition_flows(self, stratification_name, strata_names, adjustment_requests, report):
        """
        stratify flows depending on whether inflow, outflow or both need replication
        """
        pass

    def stratify_entry_flows(self, stratification_name, strata_names, requested_proportions, report):
        """
        stratify entry/recruitment/birth flows
        """
        pass

    def stratify_death_flows(self, stratification_name, strata_names, adjustment_requests, report):
        """
        add compartment-specific death flows to death data frame
        """
        pass

    def stratify_universal_death_rate(self, stratification_name, strata_names, adjustment_requests, report):
        """
        stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
        """
        pass

    def add_adjusted_parameter(self, unadjusted_parameter, stratification_name, stratum, adjustment_requests):
        """
        find the adjustment request that is relevant to a particular unadjusted parameter and stratum, otherwise allow return of null
        """
        pass

    def apply_heterogeneous_infectiousness(self, stratification_name, strata_request, infectiousness_adjustments):
        """
        work out infectiousness adjustments and set as model attributes
        """
        pass

    def set_ageing_rates(self, strata_names, report):
        """
        set intercompartmental flows for ageing from one stratum to the next
        """
        pass

    def add_stratified_flows(self, flow, stratification_name, strata_names, stratify_from, stratify_to,
                             adjustment_requests, report):
        """
        add additional stratified flow to flow data frame
        """
        pass

    def sort_absent_parameter_request(self, stratification_name, strata_names, stratum, stratify_from, stratify_to,
                                      flow):
        """
        work out what to do if a specific parameter adjustment has not been requested
        """
        pass


if __name__ == "__main__":
    sir_model = StratifiedModel(numpy.linspace(0, 60 / 365, 61).tolist(),
                         ["susceptible", "infectious", "recovered"],
                         {"infectious": 0.001},
                         {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
                         [{"type": "standard_flows", "parameter": "recovery", "from": "infectious", "to": "recovered"},
                          {"type": "infection_density", "parameter": "beta", "from": "susceptible", "to": "infectious"},
                          {"type": "compartment_death", "parameter": "infect_death", "from": "infectious"}],
                         report=False)
    sir_model.stratify("potatoes", 10, [],
                       [{"recovery": {"adjustments": {"negative": 0.7, "positive": 0.5}}},
                        {"infect_death": {"adjustments": {"negative": 0.5}}}],
                       {"negative": 0.6, "positive": 0.4}, report=True)
    sir_model.run_model()
    outputs_plot = matplotlib.pyplot.plot(sir_model.times, sir_model.outputs[:, 1])
    # matplotlib.pyplot.show()
    # print(sir_model.times)
    #
    # print(sir_model.outputs[:, 0])


