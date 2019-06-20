import numpy
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot
import copy
import pandas
from graphviz import Digraph
from sqlalchemy import create_engine
import os

# set path - sachin
os.environ["PATH"] += os.pathsep + 'C:/Users/swas0001/graphviz-2.38/release/bin'

# set path - james desktop
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin'


def find_stem(stratified_string):
    """
    find the stem of the compartment name as the text leading up to the first occurrence of "X"
    """
    first_x_location = stratified_string.find("X")
    return stratified_string if first_x_location == -1 else stratified_string[: first_x_location]


def create_stratified_name(stem, stratification_name, stratum_name):
    """
    function to generate a standardised stratified compartment name
    """
    return stem + create_stratum_name(stratification_name, stratum_name)


def create_stratum_name(stratification_name, stratum_name):
    """
    generate the name just for the particular stratification
    """
    return "X%s_%s" % (stratification_name, str(stratum_name))


def extract_x_positions(parameter):
    """
    find the positions within a string which are X and return as list reversed, including length of list
    """
    result = [loc for loc in range(len(parameter)) if parameter[loc] == "X"]
    result.append(len(parameter))
    return result


def extract_reversed_x_positions(parameter):
    """
    find the positions within a string which are X and return as list reversed, including length of list
    """
    result = extract_x_positions(parameter)
    result.reverse()
    return result


def increment_compartment(ode_equations, compartment_number, increment):
    """
    general method to increment the odes by a value specified as an argument
    """
    ode_equations[compartment_number] += increment
    return ode_equations


def normalise_dict(value_dict):
    """
    simple function to normalise the values from a list
    """
    return {key: value_dict[key] / sum(value_dict.values()) for key in value_dict}


def num_str(f):
    """
    currently unused function to convert number to a readable string
    """
    abs_f = abs(f)
    if abs_f > 1e9:
        return "%.1fB" % (f / 1e9)
    if abs_f > 1e6:
        return "%.1fM" % (f / 1e6)
    if abs_f > 1e3:
        return "%.1fK" % (f / 1e3)
    if abs_f > 100.:
        return "%.0f" % f
    if abs_f > 0.5:
        return "%.1f" % f
    if abs_f > 5e-2:
        return "%.2f" % f
    if abs_f > 5e-4:
        return "%.4f" % f
    if abs_f > 5e-6:
        return "%.6f" % f
    return str(f)


def create_flowchart(model_object, strata=-1, stratify=True, name="flow_chart"):
    """
    use graphviz module to create flow diagram of compartments and intercompartmental flows.
    """

    # set styles for graph
    styles = {"graph": {"label": "",
                        "fontsize": "16", },
              "nodes": {"fontname": "Helvetica",
                        "style": "filled",
                        "fillcolor": "#CCDDFF", },
              "edges": {"style": "dotted",
                        "arrowhead": "open",
                        "fontname": "Courier",
                        "fontsize": "10", }}

    def apply_styles(graph, styles):
        graph.graph_attr.update(("graph" in styles and styles["graph"]) or {})
        graph.node_attr.update(("nodes" in styles and styles["nodes"]) or {})
        graph.edge_attr.update(("edges" in styles and styles["edges"]) or {})
        return graph

    # find input nodes and edges
    if stratify:
        input_nodes = model_object.compartment_names
        type_of_flow = model_object.transition_flows[model_object.transition_flows.implement == strata] if strata != -1 else \
            model_object.transition_flows[model_object.transition_flows.implement == len(model_object.strata)]
    else:
        input_nodes = model_object.compartment_types
        type_of_flow = model_object.unstratified_flows

    # inputs are sectioned according to the stem value so colours can be added to each type, not yet implemented
    broken_down_nodes = []
    for stem_value in range(len(model_object.compartment_types)):
        x_vector = []
        for stem_type in range(len(input_nodes)):
            if model_object.compartment_types[stem_value] == find_stem(input_nodes[stem_type]):
                x_vector.append(input_nodes[stem_type])
        broken_down_nodes.append(stem_value)
        broken_down_nodes[stem_value] = x_vector

    model_object.flow_diagram = Digraph(format="png")
    if strata != -1:

        # find the compartment names that will be used to make the graph
        new_labels = list(set().union(type_of_flow["origin"].values, type_of_flow["to"].values))
        for label in new_labels:
            model_object.flow_diagram.node(label)
    else:
        for label in model_object.compartment_names:
            model_object.flow_diagram.node(label)

    # build the graph edges
    for row in type_of_flow.iterrows():
        model_object.flow_diagram.edge(row[1]["origin"], row[1]["to"], row[1]["parameter"])
    model_object.flow_diagram = apply_styles(model_object.flow_diagram, styles)
    model_object.flow_diagram.render(name)


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
                 default_starting_compartment="", equilibrium_stopping_tolerance=1e-6, integration_type="odeint"):

        # set flow attributes as pandas dataframes with set column names
        self.transition_flows = pandas.DataFrame(columns=["type", "parameter", "origin", "to", "implement"])
        self.death_flows = pandas.DataFrame(columns=["type", "parameter", "origin", "implement"])

        # attributes with specific format that are independent of user inputs
        self.tracked_quantities, self.output_connections, self.time_variants = \
            [{} for _ in range(3)]
        self.derived_outputs = {"times": []}
        self.compartment_values, self.compartment_names = \
            [[] for _ in range(2)]

        # features that should not be changed
        self.available_birth_approaches = ["add_crude_birth_rate", "replace_deaths", "no_births"]

        # ensure requests are fed in correctly
        self.check_and_report_attributes(
            times, compartment_types, initial_conditions, parameters, requested_flows, initial_conditions_to_total,
            infectious_compartment, birth_approach, report, reporting_sigfigs, entry_compartment,
            starting_population, default_starting_compartment, equilibrium_stopping_tolerance, integration_type)

        # stop ide complaining about attributes being defined outside __init__, even though they aren't
        self.times, self.compartment_types, self.initial_conditions, self.parameters, self.requested_flows, \
            self.initial_conditions_to_total, self.infectious_compartment, self.birth_approach, self.report, \
            self.reporting_sigfigs, self.entry_compartment, self.starting_population, \
            self.default_starting_compartment, self.default_starting_population, self.equilibrium_stopping_tolerance, \
            self.unstratified_flows, self.outputs, self.integration_type, self.flow_diagram = [None for _ in range(19)]

        # convert input arguments to model attributes
        for attribute in ["times", "compartment_types", "initial_conditions", "parameters",
                          "initial_conditions_to_total", "infectious_compartment", "birth_approach", "report",
                          "reporting_sigfigs", "entry_compartment", "starting_population",
                          "default_starting_compartment", "infectious_compartment", "equilibrium_stopping_tolerance",
                          "integration_type"]:
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
            entry_compartment, starting_population, default_starting_compartment, equilibrium_stopping_tolerance,
            integration_type):
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
                ["infectious_compartment", "birth_approach", "entry_compartment", "default_starting_compartment",
                 "integration_type"]:
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
        self.compartment_names = copy.copy(self.compartment_types)
        self.compartment_values = [0.0] * len(self.compartment_names)
        for compartment in [comp for comp in self.compartment_types if comp in self.initial_conditions]:
            self.compartment_values[self.compartment_names.index(compartment)] = self.initial_conditions[compartment]

        # sum to a total value if requested
        if initial_conditions_to_total:
            self.sum_initial_compartments_to_total()

    def sum_initial_compartments_to_total(self):
        """
        make initial conditions sum to a certain value
        """
        compartment = self.find_remainder_compartment()
        remaining_population = self.starting_population - sum(self.compartment_values)
        if remaining_population < 0.0:
            raise ValueError("total of requested compartment values is greater than the requested starting population")
        self.output_to_user("requested that total population sum to %s" % self.starting_population)
        self.output_to_user("remaining population of %s allocated to %s compartment"
                            % (remaining_population, compartment))
        self.compartment_values[self.compartment_names.index(compartment)] = remaining_population

    def find_remainder_compartment(self):
        """
        find the compartment to put the remaining population that hasn't been assigned yet when summing to total
        """
        if len(self.default_starting_compartment) > 0 and \
                self.default_starting_compartment not in self.compartment_types:
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

        # check flow requested correctly
        for flow in requested_flows:
            if flow["parameter"] not in self.parameters:
                raise ValueError("flow parameter not found in parameter list")
            if flow["origin"] not in self.compartment_types:
                raise ValueError("from compartment name not found in compartment types")
            if "to" in flow and flow["to"] not in self.compartment_types:
                raise ValueError("to compartment name not found in compartment types")

            # add flow to appropriate dataframe
            if flow["type"] == "compartment_death":
                self.add_death_flow(flow)
            else:
                self.add_transition_flow(flow)

            # add any tracked quantities that will be needed for calculating flow rates during integration
            if "infection" in flow["type"]:
                self.tracked_quantities["infectious_population"] = 0.0
            if flow["type"] == "infection_frequency":
                self.tracked_quantities["total_population"] = 0.0

    def add_default_quantities(self):
        """
        add parameters and tracked quantities that weren't requested but will be needed
        """

        # universal death rate
        if "universal_death_rate" not in self.parameters:
            self.parameters["universal_death_rate"] = 0.0

        # birth approach-specific parameters
        if self.birth_approach == "add_crude_birth_rate" and "crude_birth_rate" not in self.parameters:
            self.parameters["crude_birth_rate"] = 0.0
        elif self.birth_approach == "replace_deaths":
            self.tracked_quantities["total_deaths"] = 0.0

        # for each derived output to be recorded, initialise a tracked quantities key to zero
        for output in self.output_connections:
            self.tracked_quantities[output] = 0.0

        # parameters essential for stratification
        self.parameters["entry_fractions"] = 1.0

    def add_transition_flow(self, flow):
        """
        simply add a flow to the pandas dataframe storing the flows
        """
        flow["implement"] = 0
        self.transition_flows = self.transition_flows.append(flow, ignore_index=True)

    def add_death_flow(self, flow):
        """
        similarly for compartment-specific death flows
        """
        flow["implement"] = 0
        self.death_flows = self.death_flows.append(flow, ignore_index=True)

    """
    methods for model running
    """

    def run_model(self):
        """
        integrate model odes
        """
        self.output_to_user("now integrating")
        self.prepare_stratified_parameter_calculations()

        # basic default integration method
        if self.integration_type == "odeint":
            def make_model_function(compartment_values, time):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)
            self.outputs = odeint(make_model_function, self.compartment_values, self.times)

        # alternative integration method
        elif self.integration_type == "solve_ivp":

            # solve_ivp requires arguments to model function in the reverse order
            def make_model_function(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                return self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)

            # add a stopping condition, which was the original purpose of using this integration approach
            def set_stopping_conditions(time, compartment_values):
                self.update_tracked_quantities(compartment_values)
                net_flows = \
                    self.apply_all_flow_types_to_odes([0.0] * len(self.compartment_names), compartment_values, time)
                return max(list(map(abs, net_flows))) - self.equilibrium_stopping_tolerance
            set_stopping_conditions.terminal = True

            # solve_ivp returns more detailed structure, with (transposed) outputs (called "y") being just one component
            self.outputs = solve_ivp(
                make_model_function,
                (self.times[0], self.times[-1]), self.compartment_values, t_eval=self.times,
                events=set_stopping_conditions)["y"].transpose()

        else:
            raise ValueError("integration approach requested not available")
        self.output_to_user("integration complete")

    def prepare_stratified_parameter_calculations(self):
        """
        for use in the stratified version only
        """
        pass

    def apply_all_flow_types_to_odes(self, ode_equations, compartment_values, time):
        """
        apply all flow types to a vector of zeros (deaths must come before births in case births replace deaths)
        """
        ode_equations = self.apply_transition_flows(ode_equations, compartment_values, time)
        ode_equations = self.apply_compartment_death_flows(ode_equations, compartment_values, time)
        ode_equations = self.apply_universal_death_flow(ode_equations, compartment_values, time)
        return self.apply_birth_rate(ode_equations, compartment_values, time)

    def apply_transition_flows(self, ode_equations, compartment_values, time):
        """
        add fixed or infection-related flow to odes
        """
        for f in self.transition_flows[self.transition_flows.implement == len(self.strata)].index:

            # find adjusted parameter value
            adjusted_parameter = self.get_parameter_value(self.transition_flows.parameter[f], time)

            # find from compartment and "infectious population", which is 1 for standard flows
            infectious_population = self.find_infectious_multiplier(self.transition_flows.type[f])

            # calculate the flow and apply to the odes
            from_compartment = self.compartment_names.index(self.transition_flows.origin[f])
            net_flow = adjusted_parameter * compartment_values[from_compartment] * infectious_population
            ode_equations = increment_compartment(ode_equations, from_compartment, -net_flow)
            ode_equations = increment_compartment(
                ode_equations, self.compartment_names.index(self.transition_flows.to[f]), net_flow)

            # track any quantities dependent on flow rates
            self.track_derived_outputs(f, net_flow)

        # add another element to the derived outputs vector
        self.extend_derived_outputs(time)

        # return flow rates
        return ode_equations

    def track_derived_outputs(self, n_flow, net_flow):
        """
        calculate derived quantities to be tracked
        """
        for output_type in self.output_connections:
            if self.output_connections[output_type]["origin"] in self.transition_flows.origin[n_flow] \
                    and self.output_connections[output_type]["to"] in self.transition_flows.to[n_flow]:
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
        for f in self.death_flows[self.death_flows.implement == len(self.strata)].index:
            adjusted_parameter = self.get_parameter_value(self.death_flows.parameter[f], time)
            from_compartment = self.compartment_names.index(self.death_flows.origin[f])
            net_flow = adjusted_parameter * compartment_values[from_compartment]
            ode_equations = increment_compartment(ode_equations, from_compartment, -net_flow)
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return ode_equations

    def apply_universal_death_flow(self, ode_equations, compartment_values, time):
        """
        apply the population-wide death rate to all compartments
        """
        for compartment in self.compartment_names:
            adjusted_parameter = self.get_parameter_value("universal_death_rate", time)
            from_compartment = self.compartment_names.index(compartment)
            net_flow = adjusted_parameter * compartment_values[from_compartment]
            ode_equations = increment_compartment(ode_equations, from_compartment, -net_flow)

            # track deaths in case births are meant to replace deaths
            if "total_deaths" in self.tracked_quantities:
                self.tracked_quantities["total_deaths"] += net_flow
        return ode_equations

    def apply_birth_rate(self, ode_equations, compartment_values, time):
        """
        apply a birth rate to the entry compartments
        """
        return increment_compartment(ode_equations, self.compartment_names.index(self.entry_compartment),
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
            return self.tracked_quantities["infectious_population"] / self.tracked_quantities["total_population"]
        else:
            return 1

    def update_tracked_quantities(self, compartment_values):
        """
        update quantities that emerge during model running (not pre-defined functions of time)
        """
        for quantity in self.tracked_quantities:
            self.tracked_quantities[quantity] = 0.0
            if quantity == "infectious_population":
                self.find_infectious_population(compartment_values)
            elif quantity == "total_population":
                self.tracked_quantities["total_population"] = sum(compartment_values)

    def find_infectious_population(self, compartment_values):
        """
        calculations to find the effective infectious population
        """
        for compartment in [comp for comp in self.compartment_names if find_stem(comp) == self.infectious_compartment]:
            self.tracked_quantities["infectious_population"] += \
                compartment_values[self.compartment_names.index(compartment)]

    def get_parameter_value(self, parameter, time):
        """
        very simple, but need to split this out as a function in order to allow stratification later
        """
        return self.find_parameter_value(parameter, time)

    def store_database(self):
        """
        store outputs from the model in sql database for use in producing outputs later
        """
        engine = create_engine("sqlite:///outputs.db", echo=False)
        output_df = pandas.DataFrame(self.outputs, columns=self.compartment_names)
        output_df.insert(0, 'times', self.times)
        output_df.to_sql("outputs", con=engine, if_exists="replace", index=False)


class StratifiedModel(EpiModel):
    def add_compartment(self, new_compartment_name, new_compartment_value):
        """
        add a compartment by specifying its name and value to take
        """
        self.compartment_names.append(new_compartment_name)
        self.compartment_values.append(new_compartment_value)
        self.output_to_user("adding compartment: %s" % new_compartment_name)

    def remove_compartment(self, compartment):
        """
        remove a compartment by taking the element out of the compartment values attribute
        """
        self.removed_compartments.append(compartment)
        del self.compartment_values[self.compartment_names.index(compartment)]
        del self.compartment_names[self.compartment_names.index(compartment)]
        self.output_to_user("removing compartment: %s" % compartment)

    def __init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                 initial_conditions_to_total=True, infectious_compartment="infectious", birth_approach="no_birth",
                 report=False, reporting_sigfigs=4, entry_compartment="susceptible", starting_population=1,
                 default_starting_compartment="", equilibrium_stopping_tolerance=1e-6, integration_type="odeint"):
        EpiModel.__init__(self, times, compartment_types, initial_conditions, parameters, requested_flows,
                          initial_conditions_to_total=initial_conditions_to_total,
                          infectious_compartment=infectious_compartment, birth_approach=birth_approach,
                          report=report, reporting_sigfigs=reporting_sigfigs, entry_compartment=entry_compartment,
                          starting_population=starting_population,
                          default_starting_compartment=default_starting_compartment,
                          equilibrium_stopping_tolerance=equilibrium_stopping_tolerance,
                          integration_type=integration_type)

        self.strata, self.removed_compartments, self.overwrite_parameters, self.compartment_types_to_stratify = \
            [[] for _ in range(4)]
        self.heterogeneous_infectiousness = False
        self.infectiousness_adjustments, self.parameter_components = [{} for _ in range(2)]

    """
    pre-integration methods
    """

    def stratify(self, stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests=(),
                 requested_proportions={}, infectiousness_adjustments=(), report=True):
        """
        initial preparation and checks
        """
        strata_names, adjustment_requests = self.prepare_and_check_stratification(
            stratification_name, strata_request, compartment_types_to_stratify, adjustment_requests, report)

        # work out ageing flows (comes first so that the compartment names are still in the unstratified form)
        if stratification_name == "age":
            self.set_ageing_rates(strata_names)

        # stratify the compartments
        requested_proportions = self.tidy_starting_proportions(strata_names, requested_proportions)
        self.stratify_compartments(stratification_name, strata_names, requested_proportions)

        # stratify the flows
        self.stratify_transition_flows(stratification_name, strata_names, adjustment_requests)
        self.stratify_entry_flows(stratification_name, strata_names, requested_proportions)
        if self.death_flows.shape[0] > 0:
            self.stratify_death_flows(stratification_name, strata_names, adjustment_requests)
        self.stratify_universal_death_rate(stratification_name, strata_names, adjustment_requests)

        # heterogeneous infectiousness adjustments
        self.apply_heterogeneous_infectiousness(stratification_name, strata_request, infectiousness_adjustments)

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
            self.output_to_user("converting stratification name %s to string" % stratification_name)

        # record stratification as model attribute, find the names to apply strata and check requests
        self.strata.append(stratification_name)
        strata_names = self.find_strata_names_from_input(strata_request)
        adjustment_requests = self.alternative_adjustment_request(adjustment_requests)
        self.check_compartment_request(compartment_types_to_stratify)
        adjustment_requests = self.check_parameter_adjustment_requests(adjustment_requests, strata_names)
        return strata_names, adjustment_requests

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
        find the names of the strata to be implemented from a particular user request
        """
        if type(strata_request) == list and len(strata_request) == 0:
            raise ValueError("requested to stratify, but no strata provided")
        elif type(strata_request) == int:
            strata_names = numpy.arange(1, strata_request + 1)
            self.output_to_user("integer provided as strata labels for stratification, hence strata implemented " +
                                "are integers from 1 to %s" % strata_request)
        elif type(strata_request) == float:
            raise ValueError("number passed as request for strata labels, but not an integer greater than one, " +
                             "unclear what to do, stratification failed")
        else:
            strata_names = strata_request
        for name in range(len(strata_names)):
            strata_names[name] = str(strata_names[name])
            self.output_to_user("adding stratum: %s" % strata_names[name])
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

    def alternative_adjustment_request(self, adjustment_requests):
        """
        alternative approach to working out which parameters to overwrite - can put a capital W at the string's end
        """
        revised_adjustments = {}
        for parameter in adjustment_requests:
            revised_adjustments[parameter] = {}
            revised_adjustments[parameter]["overwrite"] = \
                adjustment_requests[parameter]["overwrite"] if "overwrite" in adjustment_requests[parameter] else []
            for stratum in adjustment_requests[parameter]:
                if stratum == "overwrite":
                    continue
                elif stratum[-1] == "W":
                    revised_adjustments[parameter][stratum[: -1]] = adjustment_requests[parameter][stratum]
                    revised_adjustments[parameter]["overwrite"].append(stratum[: -1])
                else:
                    revised_adjustments[parameter][stratum] = adjustment_requests[parameter][stratum]
            adjustment_requests[parameter] = revised_adjustments[parameter]
        return revised_adjustments

    def check_parameter_adjustment_requests(self, adjustment_requests, strata_names):
        """
        check parameter adjustments have been requested appropriately
        """
        for parameter in adjustment_requests:

            # check all the requested strata for parameter adjustments were strata that were requested
            if any([requested_stratum not in strata_names + ["overwrite"]
                    for requested_stratum in adjustment_requests[parameter]]):
                raise ValueError("stratum requested in adjustments but unavailable")

            # if any strata were not requested, assume a value of one
            for stratum in strata_names:
                if stratum not in adjustment_requests[parameter]:
                    adjustment_requests[parameter][stratum] = 1
                    self.output_to_user("no request made for adjustment to %s within stratum " % parameter +
                                        "%s so accepting parent value by default" % stratum)
        return adjustment_requests

    def tidy_starting_proportions(self, strata_names, requested_proportions):
        """
        prepare user inputs for starting proportions as needed
        """

        # assuming an equal proportion of the total for the compartment if not otherwise specified
        for stratum in strata_names:
            if stratum not in requested_proportions:
                starting_proportion = 1.0 / len(strata_names)
                requested_proportions[stratum] = starting_proportion
                self.output_to_user("no starting proportion requested for stratum %s" % stratum +
                                    " so allocated %s of total" % round(starting_proportion, self.reporting_sigfigs))
        return normalise_dict(requested_proportions)

    def stratify_compartments(self, stratification_name, strata_names, requested_proportions):
        """
        compartment stratification
        """

        # find the existing compartments that need stratification
        for compartment in \
                [comp for comp in self.compartment_names if find_stem(comp) in self.compartment_types_to_stratify]:

            # add and remove compartments
            for stratum in strata_names:
                self.add_compartment(
                    create_stratified_name(compartment, stratification_name, stratum),
                    self.compartment_values[self.compartment_names.index(compartment)] * requested_proportions[stratum])
            self.remove_compartment(compartment)

    def stratify_transition_flows(self, stratification_name, strata_names, adjustment_requests):
        """
        stratify flows depending on whether inflow, outflow or both need replication
        """
        for flow in self.transition_flows[self.transition_flows.implement == len(self.strata) - 1].index:
            self.add_stratified_flows(
                flow, stratification_name, strata_names,
                find_stem(self.transition_flows.origin[flow]) in self.compartment_types_to_stratify,
                find_stem(self.transition_flows.to[flow]) in self.compartment_types_to_stratify,
                adjustment_requests)
        self.output_to_user("stratified transition flows matrix:\n%s" % self.transition_flows)

    def stratify_entry_flows(self, stratification_name, strata_names, requested_proportions):
        """
        stratify entry/recruitment/birth flows
        """
        entry_fractions = {}

        # work out parameter values for stratifying the entry proportion adjustments
        if self.entry_compartment in self.compartment_types_to_stratify:
            for stratum in strata_names:
                entry_fraction_name = create_stratified_name("entry_fraction", stratification_name, stratum)
                if stratification_name == "age" and str(stratum) == "0":
                    entry_fractions[entry_fraction_name] = 1.0
                    continue
                elif stratification_name == "age":
                    entry_fractions[entry_fraction_name] = 0.0
                    continue
                elif "adjustments" in requested_proportions and stratum in requested_proportions["adjustments"]:
                    entry_fractions[entry_fraction_name] = requested_proportions["adjustments"][stratum]
                    self.output_to_user("assigning specified proportion of starting population to %s" % stratum)
                else:
                    entry_fractions[entry_fraction_name] = 1.0 / len(strata_names)
                    self.output_to_user("assuming %s " % entry_fractions[entry_fraction_name] +
                                        "of starting population to be assigned to %s stratum by default" % stratum)
        self.parameters.update(normalise_dict(entry_fractions))

    def stratify_death_flows(self, stratification_name, strata_names, adjustment_requests):
        """
        add compartment-specific death flows to death data frame
        """
        for flow in self.death_flows[self.death_flows.implement == len(self.strata) - 1].index:
            for stratum in strata_names:
                parameter_name = self.add_adjusted_parameter(
                    self.death_flows.parameter[flow], stratification_name, stratum, adjustment_requests)
                if not parameter_name:
                    parameter_name = self.death_flows.parameter[flow]
                self.death_flows = self.death_flows.append(
                    {"type": self.death_flows.type[flow],
                     "parameter": parameter_name,
                     "origin": create_stratified_name(self.death_flows.origin[flow], stratification_name, stratum),
                     "implement": len(self.strata)},
                    ignore_index=True)

    def stratify_universal_death_rate(self, stratification_name, strata_names, adjustment_requests):
        """
        stratify the approach to universal, population-wide deaths (which can be made to vary by stratum)
        adjust each parameter that refers to the universal death rate according to user request
        """
        for parameter in [param for param in self.parameters if find_stem(param) == "universal_death_rate"]:
            for stratum in strata_names:
                self.add_adjusted_parameter(parameter, stratification_name, stratum, adjustment_requests)

    def add_adjusted_parameter(self, unadjusted_parameter, stratification_name, stratum, adjustment_requests):
        """
        find the adjustment request that is relevant to a particular unadjusted parameter and stratum
        otherwise allow return of None
        """
        parameter_adjustment_name = None

        # find the adjustment request that is an extension of the base parameter type being considered
        for parameter_request in [req for req in adjustment_requests if unadjusted_parameter.startswith(req)]:
            self.output_to_user(
                "modifying %s for %s stratum of %s" % (unadjusted_parameter, stratum, stratification_name))
            parameter_adjustment_name = create_stratified_name(unadjusted_parameter, stratification_name, stratum)

            # implement user request (otherwise parameter will be left out and assumed to be 1 during integration)
            if stratum in adjustment_requests[parameter_request]:
                self.parameters[parameter_adjustment_name] = adjustment_requests[parameter_request][stratum]

            # overwrite parameters higher up the tree by listing which ones to be overwritten
            if "overwrite" in adjustment_requests[parameter_request] and \
                    stratum in adjustment_requests[parameter_request]["overwrite"]:
                self.overwrite_parameters.append(parameter_adjustment_name)
        return parameter_adjustment_name

    def apply_heterogeneous_infectiousness(self, stratification_name, strata_request, infectiousness_adjustments):
        """
        work out infectiousness adjustments and set as model attributes
        """
        if len(infectiousness_adjustments) == 0:
            self.output_to_user("heterogeneous infectiousness not requested for this stratification")
        elif self.infectious_compartment not in self.compartment_types_to_stratify:
            raise ValueError("request for infectiousness stratification does not apply to the infectious compartment")
        else:
            self.heterogeneous_infectiousness = True
            for stratum in infectiousness_adjustments:
                if stratum not in strata_request:
                    raise ValueError("stratum to have infectiousness modified not found within requested strata")
                adjustment_name = create_stratified_name("", stratification_name, stratum)
                self.infectiousness_adjustments[adjustment_name] = infectiousness_adjustments[stratum]

    def set_ageing_rates(self, strata_names):
        """
        set intercompartmental flows for ageing from one stratum to the next
        """
        for stratum_number in range(len(strata_names[: -1])):
            start_age = int(strata_names[stratum_number])
            end_age = int(strata_names[stratum_number + 1])
            ageing_parameter_name = "ageing%sto%s" % (start_age, end_age)
            ageing_rate = 1.0 / (end_age - start_age)
            self.output_to_user("ageing rate from age group %s to %s is %s"
                                % (start_age, end_age, round(ageing_rate, self.reporting_sigfigs)))
            self.parameters[ageing_parameter_name] = ageing_rate
            for compartment in self.compartment_names:
                self.transition_flows = self.transition_flows.append(
                    {"type": "standard_flows",
                     "parameter": ageing_parameter_name,
                     "origin": create_stratified_name(compartment, "age", start_age),
                     "to": create_stratified_name(compartment, "age", end_age),
                     "implement": len(self.strata)},
                    ignore_index=True)

    def add_stratified_flows(self, flow, stratification_name, strata_names, stratify_from, stratify_to,
                             adjustment_requests):
        """
        add additional stratified flow to flow data frame
        """
        if stratify_from or stratify_to:
            self.output_to_user(
                "for flow from %s to %s in stratification %s"
                % (self.transition_flows.origin[flow], self.transition_flows.to[flow], stratification_name))

            # loop over each stratum in the requested stratification structure
            for stratum in strata_names:

                # find parameter name
                parameter_name = self.add_adjusted_parameter(
                    self.transition_flows.parameter[flow], stratification_name, stratum, adjustment_requests)
                if not parameter_name:
                    parameter_name = self.sort_absent_parameter_request(
                        stratification_name, strata_names, stratum, stratify_from, stratify_to, flow)
                self.output_to_user("\tadding parameter %s" % parameter_name)

                # determine whether to and/or from compartments are stratified
                from_compartment = \
                    create_stratified_name(self.transition_flows.origin[flow], stratification_name, stratum) if \
                        stratify_from else self.transition_flows.origin[flow]
                to_compartment = \
                    create_stratified_name(self.transition_flows.to[flow], stratification_name, stratum) if \
                        stratify_to else self.transition_flows.to[flow]

                # add the new flow
                self.transition_flows = self.transition_flows.append(
                    {"type": self.transition_flows.type[flow],
                     "parameter": parameter_name,
                     "origin": from_compartment,
                     "to": to_compartment,
                     "implement": len(self.strata)},
                    ignore_index=True)

    def sort_absent_parameter_request(self, stratification_name, strata_names, stratum, stratify_from, stratify_to,
                                      flow):
        """
        work out what to do if a specific parameter adjustment has not been requested
        """

        # default behaviour if not specified is to split the parameter into equal parts if to compartment is split
        if not stratify_from and stratify_to:
            self.output_to_user("\tsplitting existing parameter value %s into %s equal parts"
                                % (self.transition_flows.parameter[flow], len(strata_names)))
            self.parameters[
                create_stratified_name(self.transition_flows.parameter[flow], stratification_name, stratum)] = \
                1.0 / len(strata_names)

        # otherwise if no request, retain the existing parameter
        else:
            parameter_name = self.transition_flows.parameter[flow]
            self.output_to_user("\tretaining existing parameter value %s" % parameter_name)
        return parameter_name

    """
    methods to be called during the process of model running
    """

    def prepare_stratified_parameter_calculations(self):
        """
        prior to integration commencing, work out what the components are of each parameter being implemented
        """

        # create list of all the parameters that we need to find the list of adjustments for
        parameters_to_adjust = []
        for f in range(self.transition_flows.shape[0]):
            if self.transition_flows.implement[f] == len(self.strata) and \
                    self.transition_flows.parameter[f] not in parameters_to_adjust:
                parameters_to_adjust.append(self.transition_flows.parameter[f])
        for f in range(self.death_flows.shape[0]):
            if self.death_flows.implement[f] == len(self.strata) and \
                    self.death_flows.parameter[f] not in parameters_to_adjust:
                parameters_to_adjust.append(self.death_flows.parameter[f])
        parameters_to_adjust.append("universal_death_rate")
        for parameter in parameters_to_adjust:
            self.find_parameter_components(parameter)

    def find_parameter_components(self, parameter):
        """
        extract the components of the stratified parameter into a list structure
        """

        # collate all the parameter components into time-variant or constant
        self.parameter_components[parameter] = {"time_variants": [], "constants": [], "constant_value": 1}
        for x_instance in extract_reversed_x_positions(parameter):
            component = parameter[: x_instance]
            is_time_variant = component in self.time_variants
            if component in self.overwrite_parameters and is_time_variant:
                self.parameter_components[parameter] = \
                    {"time_variants": [component], "constants": [], "constant_value": 1}
                break
            elif component in self.overwrite_parameters and not is_time_variant:
                self.parameter_components[parameter] = \
                    {"time_variants": [], "constants": [component], "constant_value": 1}
                break
            elif is_time_variant:
                self.parameter_components[parameter]["time_variants"].append(component)
            elif component in self.parameters:
                self.parameter_components[parameter]["constants"].append(component)

        # pre-calculate the constant component by multiplying through all the constant values
        for constant_parameter in self.parameter_components[parameter]["constants"]:
            self.parameter_components[parameter]["constant_value"] *= self.parameters[constant_parameter]

    def get_parameter_value(self, parameter, time):
        """
        calculate adjusted parameter value from pre-calculated product of constant components and time variants
        """
        adjusted_parameter = self.parameter_components[parameter]["constant_value"]
        for time_variant in self.parameter_components[parameter]["time_variants"]:
            adjusted_parameter *= self.time_variants[time_variant](time)
        return adjusted_parameter

    def find_infectious_population(self, compartment_values):
        """
        calculations to find the effective infectious population
        """

        # loop through all compartments and find the ones representing active infectious disease
        for compartment in self.compartment_names:
            if find_stem(compartment) == self.infectious_compartment:

                # assume homogeneous infectiousness until requested otherwise
                infectiousness_modifier = 1.0

                # haven't yet done heterogeneous infectiousness
                if self.heterogeneous_infectiousness:
                    for adjustment in [adj for adj in self.infectiousness_adjustments if adj in compartment]:
                        infectiousness_modifier = self.infectiousness_adjustments[adjustment]

                self.tracked_quantities["infectious_population"] += \
                    compartment_values[self.compartment_names.index(compartment)] * \
                    infectiousness_modifier

    def apply_birth_rate(self, ode_equations, compartment_values, time):
        """
        apply a population-wide death rate to all compartments
        """
        total_births = self.find_total_births(compartment_values)

        # split the total births across entry compartments
        for compartment in self.compartment_names:
            if find_stem(compartment) == self.entry_compartment:

                # calculate adjustment to original stem entry rate
                entry_fraction = 1.0
                x_positions = extract_x_positions(compartment)

                if len(x_positions) > 1:
                    for x_instance in range(len(x_positions) - 1):
                        entry_fraction *= \
                            self.parameters["entry_fractionX%s"
                                            % compartment[x_positions[x_instance] + 1: x_positions[x_instance + 1]]]
                compartment_births = entry_fraction * total_births
                ode_equations = increment_compartment(
                    ode_equations, self.compartment_names.index(compartment), compartment_births)
        return ode_equations


if __name__ == "__main__":
    sir_model = StratifiedModel(numpy.linspace(0, 60 / 365, 61).tolist(),
                         ["susceptible", "infectious", "recovered"],
                         {"infectious": 0.001},
                         {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
                         [{"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
                          {"type": "infection_density", "parameter": "beta", "origin": "susceptible", "to": "infectious"},
                          {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}],
                         report=False, integration_type="solve_ivp")
    sir_model.stratify("hiv", ["negative", "positive"], [],
                       {"recovery": {"negative": 0.7, "positive": 0.5},
                        "infect_death": {"negative": 0.5}},
                       {"negative": 0.6}, report=False)
    # sir_model.stratify("age", [1, 10, 3], [], {}, report=False)

    sir_model.run_model()

    create_flowchart(sir_model, strata=len(sir_model.strata))

    outputs_plot = matplotlib.pyplot.plot(sir_model.times, sir_model.outputs[:, 2] + sir_model.outputs[:, 3])

    # print(len(sir_model.times))
    # print(sir_model.outputs[1])

    # print(sir_model.outputs)

    print(sir_model.transition_flows)

    matplotlib.pyplot.show()
    # print(sir_model.times)
    #
    # print(sir_model.outputs[:, 0])


