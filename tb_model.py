
import summer_model
import matplotlib.pyplot
import os


def add_vtp_latency_parameters(parameters_, change_time_unit=365.25):
    """
    function to add the latency parameters estimated by Ragonnet et al from our paper in Epidemics to the existing
    parameter dictionary
    """
    vtp_latency_parameters = \
        {"early_progression": 1.1e-3,
         "stabilisation": 1.0e-2,
         "late_progression": 5.5e-6}
    parameters_.update({key: value * change_time_unit for key, value in vtp_latency_parameters.items()})
    return parameters_


def get_age_specific_latency_parameters(parameter, unit_change=365.25):
    """
    get the age-specific latency parameters estimated by Ragonnet et al
    """
    age_stratified_parameters = \
        {"early_progression": {"0W": 6.6e-3, "5W": 2.7e-3, "15W": 2.7e-4},
         "stabilisation": {"0W": 1.2e-2, "5W": 1.2e-2, "15W": 5.4e-3},
         "late_progression": {"0W": 1.9e-11, "5W": 6.4e-6, "15W": 3.3e-6}}
    return {key: value * unit_change for key, value in age_stratified_parameters[parameter].items()}


def get_all_age_specific_latency_parameters(parameters_=("early_progression", "stabilisation", "late_progression")):
    """
    collate all the latency parameters together from the previous function
    """
    return {parameter: get_age_specific_latency_parameters(parameter) for parameter in parameters_}


def add_standard_latency_flows(flows_):
    """
    adds our standard latency flows to the list of flows to be implemented in the model
    """
    flows_ += [
        {"type": "standard_flows", "parameter": "early_progression", "origin": "early_latent", "to": "infectious"},
        {"type": "standard_flows", "parameter": "stabilisation", "origin": "early_latent", "to": "late_latent"},
        {"type": "standard_flows", "parameter": "late_progression", "origin": "late_latent", "to": "infectious"}]
    return flows_


def sinusoidal_scaling_function(start_time, baseline_value, end_time, final_value):
    """
    with a view to implementing scale-up functions over time, use the cosine function to produce smooth scale-up
    functions from one point to another
    """
    def sinusoidal_function(x):
        if not isinstance(x, float):
            raise ValueError("value fed into scaling function not a float")
        elif start_time > end_time:
            raise ValueError("start time is later than end time")
        elif x < start_time:
            return baseline_value
        elif start_time < x < end_time:
            return baseline_value + \
                   (final_value - baseline_value) * \
                   (0.5 - 0.5 * numpy.cos((x - start_time) * numpy.pi / (end_time - start_time)))
        else:
            return final_value
    return sinusoidal_function


def convert_competing_proportion_to_rate(competing_flows):
    """
    convert a proportion to a rate dependent on the other flows coming out of a compartment
    """
    return lambda proportion: proportion * competing_flows / (1.0 - proportion)


def return_function_of_function(inner_function, outer_function):
    """
    general method to return a chained function from two functions
    """
    return lambda value: outer_function(inner_function(value))


if __name__ == "__main__":

    # set basic parameters, flows and times, except for latency flows and parameters, then functionally add latency
    case_fatality_rate = 0.4
    untreated_disease_duration = 3.0
    parameters = \
        {"beta": 10.0,
         "recovery": case_fatality_rate / untreated_disease_duration,
         "infect_death": (1.0 - case_fatality_rate) / untreated_disease_duration,
         "universal_death_rate": 1.0 / 50.0,
         "case_detection": 0.0}
    parameters = add_vtp_latency_parameters(parameters)

    times = numpy.linspace(1800., 2020.0, 201).tolist()
    flows = [{"type": "infection_frequency", "parameter": "beta", "origin": "susceptible", "to": "early_latent"},
             {"type": "infection_frequency", "parameter": "beta", "origin": "recovered", "to": "early_latent"},
             {"type": "standard_flows", "parameter": "recovery", "origin": "infectious", "to": "recovered"},
             {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"}]
    flows = add_standard_latency_flows(flows)

    tb_model = summer_model.StratifiedModel(
        times, ["susceptible", "early_latent", "late_latent", "infectious", "recovered"], {"infectious": 1e-3},
        parameters, flows, birth_approach="replace_deaths")

    tb_model.add_transition_flow(
        {"type": "standard_flows", "parameter": "case_detection", "origin": "infectious", "to": "recovered"})

    cdr_scaleup = sinusoidal_scaling_function(1950.0, 0.0, 2010.0, 0.6)
    prop_to_rate = convert_competing_proportion_to_rate(1.0 / untreated_disease_duration)
    detect_rate = return_function_of_function(cdr_scaleup, prop_to_rate)

    tb_model.time_variants["case_detection"] = detect_rate

    print(get_all_age_specific_latency_parameters())
    tb_model.stratify("age", [5, 15], [],
                      adjustment_requests=get_all_age_specific_latency_parameters(),
                      report=False)
    tb_model.run_model()

    print(os.getcwd())
    # tb_model.transition_flows.to_csv("_.csv")

    # get outputs
    infectious_population = tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_0")] + \
                            tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_5")] + \
                            tb_model.outputs[:, tb_model.compartment_names.index("infectiousXage_15")]

    matplotlib.pyplot.plot(times, infectious_population * 1e5)
    # print(infectious_population * 1e5)

    # tb_model.death_flows.to_csv("tb_model_deaths.csv")

    matplotlib.pyplot.xlim((1950., 2010.))
    matplotlib.pyplot.ylim((0.0, 2000.0))
    # matplotlib.pyplot.show()
