import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker


# Configure matplotlib for LaTeX-style fonts (without unicode issues)
nice_fonts = {
    "text.usetex": False,  # Disable to avoid Unicode/LaTeX issues
    "font.family": "serif",
    "axes.labelsize": 14,
    "font.size": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}
matplotlib.rcParams.update(nice_fonts)

# Set up CoolProp for nitrogen
fluid = 'N2'
Water = CP.AbstractState("HEOS", fluid)

pc = Water.keyed_output(CP.iP_critical)
Tc = Water.keyed_output(CP.iT_critical)
pmax = Water.keyed_output(CP.iP_max)
print(f"pc = {pc:.0f} ; Tc= {Tc:.0f}")

# Temperature range within valid saturation limits
Ts = np.linspace(64, Tc, 500)

# Calculate saturated liquid and vapor densities
rhoL = np.array([CP.CoolProp.PropsSI('D', 'T', T, 'Q', 0, fluid) for T in Ts])
rhoV = np.array([CP.CoolProp.PropsSI('D', 'T', T, 'Q', 1, fluid) for T in Ts])
ps = np.array([CP.CoolProp.PropsSI('P', 'T', T, 'Q', 0, fluid) for T in Ts])

# Set up plot for the diagram without the dots
fig = plt.figure(figsize=(12, 9))  # Increase figure size for a larger graph
ax = fig.add_subplot(111)
lw = 2

# Plot saturation dome (P vs rho)
ax.plot(rhoL, ps, color='blue', lw=lw, label="Saturated Liquid")
ax.plot(rhoV, ps, color='red', lw=lw, label="Saturated Vapor")

# Critical point
dc = CP.CoolProp.PropsSI('D', 'T', Tc, 'P', pc, fluid)
ax.plot(dc, pc, 'ko', label='Critical Point')

# Axes styling
ax.set_xlabel("Density $\\rho$ [kg/m$^3$]")
ax.set_ylabel("Pressure $P$ [bar]")  # Change unit to bar
ax.grid(True)

# Format Y-axis ticks in bar
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x/1e5:.1f}'))

ax.grid(True)

# Critical lines
ax.axhline(pc, color='gray', linestyle='--', lw=1)
ax.axvline(dc, color='gray', linestyle='--', lw=1)

# Simulated Points
Ts_refuel = [76, 76, 76, 76]
Ps_refuel = [3.1e5, 21.1e5, 31.1e5, 39.1e5]
# Calculate density for these points using CoolProp
rho_refuel = [CP.CoolProp.PropsSI('D', 'T', T, 'P', P, fluid) for T, P in zip(Ts_refuel, Ps_refuel)]

# Plot the refuelling points as orange dots
ax.plot(rho_refuel, Ps_refuel, "o", ms=6, color="orange", label="LN2 refuelling start")

# Add refuelling end point (T_end = 293 K, P_end = 1.1e5 Pa)
T_end = 293  # K
P_end = 1.1e5  # Pa
rho_end = CP.CoolProp.PropsSI('D', 'T', T_end, 'P', P_end, fluid)

# Plot the refuelling end point as a red dot
ax.plot(rho_end, P_end, "o", ms=6, color="red", label="Refuelling end")

# Position the legend inside the graph (top left)
ax.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95), fontsize=12)  # Legend inside the plot

# Save the plot with dots (including refuelling end)
plt.tight_layout()
plt.savefig("N2_with_Ranges_And_AllPoints.png", dpi=1000)

# Show the plot
plt.show()
