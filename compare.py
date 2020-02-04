import numpy as np
from matplotlib import pyplot as plt
from scipy import constants
import ballistic


def get_mass_ratio_ballistic(R, I_sp):
    """Get the mass ratio m_initial / m_final for a vehicle
    flying a suborbital ballistic trajectory.

    Arguments:
        R: Ground distance of trajectory [units: meter].
        I_sp: Specific impulse of the propulsion system [units: s]
    """
    # Range angle of the free-flight portion of the trajectory.
    # [units: rad]
    psi = R / ballistic.R_earth
    # Assumed altitude at burnout
    # [units: m]
    alt_bo = 100e3

    # Burnout velocity
    # [units: m s^-1]
    v_bo = ballistic.get_burnout_velocity_max_range(psi, alt_bo)
    # Assumed extra delta-v needed for losses, landing, etc.
    # [units: m s^-1]
    dv_extra = 500.
    # Total delta-v
    dv = v_bo + dv_extra

    # rocket equation
    mass_ratio = np.exp(dv / (I_sp * constants.g))
    return mass_ratio


def get_flight_time_ballistic(R):
    """Get the flight time for a vehicle
    flying a suborbital ballistic trajectory.

    Arguments:
        R: Ground distance of trajectory [units: meter].
    """
    # Range angle of the free-flight portion of the trajectory.
    # [units: rad]
    psi = R / ballistic.R_earth
    # Assumed altitude at burnout
    # [units: m]
    alt_bo = 100e3

    return ballistic.get_time_of_free_flight_max_range(psi, alt_bo)


def get_mass_ratio_aircraft(R, I_sp, LD, v_cruise):
    """Get the mass ratio m_initial / m_final for a vehicle
    flying in the atmosphere in steady, level flight.

    Arguments:
        R: Ground distance of trajectory [units: meter].
        I_sp: Specific impulse of the propulsion system [units: s].
        LD: Lift/drag ratio [units: dimensionless].
        v_cruise: Cruise velocity [units: m s^-1].
    """
    # Assumed extra 'range' to account for climb, fuel reserves
    R_extra = 200e3
    R_eff = R + R_extra

    # Breguet range equation
    mass_ratio = np.exp(R / (v_cruise * I_sp * LD))
    return mass_ratio


def get_flight_time_aircraft(R, v_cruise):
    """Get the flight time for a vehicle
    flying in the atmosphere in steady, level flight.

    Arguments:
        R: Ground distance of trajectory [units: meter].

        v_cruise: Cruise velocity [units: m s^-1].
    """
    return R / v_cruise


def main():
    # Distance examples [units: km]
    distance_examples = {
        'New York - London': 5571,
        'SF - Tokyo': 8266,
        'New York - Shanghai': 11860,
        'New York - Sydney': 15989,
        'London - Auckland': 18337,
    }
    # Range
    R = np.linspace(500e3, np.pi * ballistic.R_earth)
    # Fuel mass per mass of vehicle at landing
    # [units: dimensionless]
    m_fuel = {}
    # Flight time [units: s]
    t_flight = {}

    ### Sub-orbital ballistic rocket ###
    # Values for a methane/oxygen engine
    I_sp = 340.
    o_f_ratio = 3.8
    mass_ratio_rocket = get_mass_ratio_ballistic(R, I_sp)
    # Propellant mass per mass of vehicle at landing
    # [units: dimensionless]
    m_prop_rocket = mass_ratio_rocket - 1
    m_fuel['rocket'] = m_prop_rocket / (1 + o_f_ratio)
    t_flight['rocket'] = get_flight_time_ballistic(R)

    ### Subsonic aircraft ###
    I_sp = 6000.
    LD = 20.
    v_cruise = 240.
    mass_ratio = get_mass_ratio_aircraft(R, I_sp, LD, v_cruise)
    m_fuel['subsonic aircraft'] = mass_ratio - 1
    t_flight['subsonic aircraft'] = get_flight_time_aircraft(R, v_cruise)

    ### Mach 2 supersonic aircraft ###
    # based on Concorde
    I_sp = 3010.
    LD = 7.14
    v_cruise = 599.
    mass_ratio = get_mass_ratio_aircraft(R, I_sp, LD, v_cruise)
    m_fuel['M2 aircraft'] = mass_ratio - 1
    t_flight['M2 aircraft'] = get_flight_time_aircraft(R, v_cruise)


    ### Plot ###
    colors = {
        'rocket': 'tab:orange',
        'subsonic aircraft': 'tab:blue',
        'M2 aircraft': 'tab:purple',
    }
    labels = {
        'rocket': 'Rocket, suborbital',
        'subsonic aircraft': 'Aircraft, subsonic',
        'M2 aircraft': 'Aircraft, Mach 2',
    }
    fig, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True, figsize=(6, 6))
    ax_fuel, ax_time = axes

    for key in colors:
        ax_fuel.plot(
            1e-3 * R, m_fuel[key],
            color=colors[key], label=labels[key])
        ax_time.plot(
            1e-3 * R, t_flight[key] / 3600.,
            color=colors[key], label=labels[key])

    # Plot Concorde max range point for comparison
    concorde_max_to_mass = 185070.
    concorde_max_fuel_mass = 95680.
    concorde_m_fuel_m_land = concorde_max_fuel_mass / (concorde_max_to_mass - concorde_max_fuel_mass)
    concorde_range = 7222e3
    ax_fuel.scatter(
        1e-3 * concorde_range, concorde_m_fuel_m_land,
        color=colors['M2 aircraft'], marker='x')
    ax_fuel.text(
        x=1e-3 * concorde_range, y=concorde_m_fuel_m_land + 0.1,
        s='Concorde', color=colors['M2 aircraft'])

    ax_fuel.text(
        x=6500, y=2.1, s='Only fuel, not oxidizer',
        color=colors['rocket'])

    ax_fuel.set_ylabel('Mass fuel / mass of vehicle at landing [-]')
    ax_fuel.legend()
    ax_fuel.set_ylim([0, ax_fuel.get_ylim()[1]])

    ax_time.set_ylabel('Flight time [hr]')
    ax_time.set_xlabel('Range [km]')
    ax_time.set_ylim([0, ax_time.get_ylim()[1]])
    for key, dist in distance_examples.items():
        ax_time.text(
            x=dist, y=2, s=key,
            color='grey', rotation=90)
    ax_time.legend()

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)


if __name__ == '__main__':
    main()
    plt.show()
