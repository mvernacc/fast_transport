"""Calculations for ballistic suborbital trajectories.

References:
  [1] Roger Bate, Donald Mueller and Jerry White. "Fundamentals of Astrodynamics"
    Dover Publications, 1971.
"""
import numpy as np
from matplotlib import pyplot as plt

# Radius of the Earth
# [units: m]
R_earth = 6371e3
# Standard gravitational parameter of the Earth
# [units: m^3 s^-2]
mu_earth = 3.986e14


def get_burnout_velocity_max_range(psi, alt_bo):
    # Eq. 6.2-20 from [1]:
    # Q is dimensionless.
    Q_bo = 2 * np.sin(psi / 2) / (1  + np.sin(psi / 2))
    # Burnout radius from Earth's center
    # [units: m]
    r_bo = alt_bo + R_earth
    # velocity at burnout [units: m s^-1].
    v_bo = (Q_bo * mu_earth / r_bo)**.5
    return v_bo


def get_time_of_free_flight_max_range(psi, alt_bo):
    # Eq. 6.2-20 from [1]:
    Q_bo = 2 * np.sin(psi / 2) / (1  + np.sin(psi / 2))
    # velocity at burnout [units: m s^-1].
    v_bo = get_burnout_velocity_max_range(psi, alt_bo)
    # Burnout radius from Earth's center
    # [units: m]
    r_bo = alt_bo + R_earth
    # Flight path angle at burnout
    # Equation 6.2-18 from [1]
    # [units: radian]
    phi_bo = 0.25 * (np.pi - psi)
    # Eccentricity of the orbit
    # Equation 6.2-11 from [1]
    # [units: dimensionless]
    e = (1 + Q_bo * (Q_bo - 2) * np.cos(phi_bo)**2)**0.5
    # Semi-major axis of the orbit
    # Equation 6.2-3 from [1]
    # [units: m]
    a = r_bo / (2 - Q_bo)
    # Equation 6.2-21 from [1]
    E_1 = (e - np.cos(psi / 2)) / (1 - e * np.cos(psi / 2))
    # Time of flight
    # Equation 6.2-22 from [1]
    # [units: s]
    t_ff = 2 * (a**3 / mu_earth)**0.5 * (np.pi - E_1 + e * np.sin(E_1))
    return t_ff


def demo():
    psi = np.linspace(0.01, np.pi)
    alt_bo = 100e3
    v_bo = get_burnout_velocity_max_range(psi, alt_bo)
    t_ff = get_time_of_free_flight_max_range(psi, alt_bo)

    fig, axes = plt.subplots(
        nrows=2, ncols=1, sharex=True,
        figsize=(6, 6))
    axes[0].plot(psi, 1e-3 * v_bo)
    axes[0].set_ylabel('Burnout velocity $v_{bo}$ [km s$^{-1}$]')
    axes[0].set_title('Max-range ballistic trajectories')
    axes[0].set_ylim([0, axes[0].get_ylim()[1]])

    axes[1].plot(psi, t_ff / 60)
    axes[1].set_ylabel('Time of free flight $t_{ff}$ [min]')
    axes[1].set_xlabel('Range angle $\\Psi$ [rad]')
    axes[1].set_ylim([0, axes[1].get_ylim()[1]])
    fig.tight_layout()

if __name__ == '__main__':
    demo()
    plt.show()
