import matplotlib
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib.patches as patches 
import matplotlib.lines as mlines
import matplotlib.patches as mpatches 


nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 14,
        "font.size": 14,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
}

matplotlib.rcParams.update(nice_fonts)

#Water = CP.AbstractState("HEOS", "H2")
Water = CP.AbstractState("HEOS", "N2")

pc = Water.keyed_output(CP.iP_critical) # critical pressure
Tc = Water.keyed_output(CP.iT_critical) # critical temperature
print("pc = %d ; Tc= %d" %(pc,Tc))
Tmin = 50
Tmax = 300
pmax = Water.keyed_output(CP.iP_max) 


fig = plt.figure(figsize = (11,8))
ax = fig.add_subplot(111)
lw = 2

# ----------------
# Saturation curve
# ----------------
Ts = np.linspace(5, Tc, 1000)
ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,'N2')

# ------
# Labels this is the saturation line
# ------
plt.plot(Ts,ps, 'grey', zorder=1, lw = lw, solid_capstyle = 'round')

# Critical lines
plt.axvline(Tc, dashes = [2, 2])
plt.axhline(pc, dashes = [2, 2])

# Labels for regions
plt.text(1.9*Tc, 3.2*pc, 'supercritical gas',ha= 'center')
plt.text(1.9*Tc, 0.6*pc, 'superheated gas', rotation = 0,ha= 'center')#,color="magenta"
plt.text(0.55*Tc, 3.2*pc, 'supercritical liquid', rotation = 0, ha = 'center')
plt.text(0.55*Tc, 0.2*pc, 'liquid', rotation = 45)
plt.text(0.75*Tc, 0.05*pc, 'gas', rotation = 45) 

# ------
# Plot
# ------
plt.ylim(5e3,5e7)
plt.gca().set_yscale('log')
plt.gca().set_xlim(Tmin, Tmax)
plt.gca().set_xscale('log')
plt.ylabel('Pressure [Pa]')
plt.xlabel('Temperature [K]')

# ------
# Add regions of interest (operating regions for Gas turbines, etc.....)
# ------

plt.tight_layout()
matplotlib.pyplot.savefig(
    'N2_phase_diagram.png',
    dpi=1000
)


from matplotlib.patches import Polygon
vertices = [(120, 33e5), (132, 33e5), (132, 35e5), (120, 35e5)]
polygon = patches.Polygon(vertices, closed=True, zorder=2, color='#BD5FA7', alpha=0.5)
patches , colors_p = [polygon], ['#BD5FA7']

from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
our_cmap = ListedColormap(colors_p)

p = PatchCollection(patches, match_original=True)

ax.add_collection(p)

ax.legend(handles=patches, bbox_to_anchor=(0.99, 0.01), loc='lower right',ncols=1,  borderaxespad=0.)

# ------

plt.tight_layout()
matplotlib.pyplot.savefig(
    'N2_with_Ranges.png',
    dpi=1000
)

# ------
# Add Simulated Points
# ------
from matplotlib.legend_handler import HandlerLine2D
xs = [76, 76, 76, 76]
ys = [3.1e5, 21.1e5, 31.1e5, 39.1e5]
plt.plot(xs,ys,"o",ms=3,color ="orange" , label = "LN2 refuelling start")

T_gas = 293  # kg/m^3
P_gas = 1.1e5  # Pa

# Find the rho using CoolProp
rho_gas = CP.CoolProp.PropsSI('D', 'T', T_gas, 'P', P_gas, 'N2')

print("T_gas = %d ; P_gas = %d ; rho_gas = %.3f" % (T_gas, P_gas, rho_gas))

# Given temperature (same for all cases)
T_in = 76  # Kelvin
P_in = [3.1e5, 21.1e5, 31.1e5, 39.1e5]  # Pascal

# Loop over the cases and compute rho
for i, P in enumerate(P_in):
    rho_in = CP.CoolProp.PropsSI('D', 'T', T_in, 'P', P, 'N2')
    dynamic_viscosity = CP.CoolProp.PropsSI('V', 'T', T_in, 'P', P, 'N2')  # Dynamic viscosity (Pa.s)
    print(f"T_in{i+1} = {T_in} K, P_in{i+1} = {P:.1f} Pa, rho_in{i+1} = {rho_in:.3f} kg/m³, Dynamic Viscosity (T_in{i+1}) = {dynamic_viscosity:.6e} Pa.s")
    
    


xt = [293, 293, 293, 293]
yt = [1.1e5, 1.1e5, 1.1e5, 1.1e5]
plt.plot(xt,yt,"o",ms=3,color ="red")

ax.add_patch(polygon)

start_point_legend = mlines.Line2D([], [], color='orange', marker='o', linestyle='None', markersize=6, alpha=0.6, label='LN2 refuelling start')
end_point_legend = mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=6, alpha=0.6, label='LN2 refuelling end')


ax.add_patch(polygon)
gtinj_legend = mpatches.Patch(color='#BD5FA7', alpha=0.5, label='Supercritical region of focus')

#####################################################################################################################

#arrow_near1 = mpatches.FancyArrowPatch((T_in, P_in[1]), (xt[1], yt[1]), mutation_scale=10,
#                                   color="blue",alpha=0.6) 
#arrow_near2 = mpatches.FancyArrowPatch((T_in, P_in[2]), (xt[2], yt[2]), mutation_scale=10,
#                                   color="blue",alpha=0.6) 
#arrow_near3 = mpatches.FancyArrowPatch((T_in, P_in[3]), (xt[3], yt[3]), mutation_scale=10,
#                                   color="blue",alpha=0.6) 
#arrow_near4 = mpatches.FancyArrowPatch((T_in, P_in[0]), (xt[0], yt[0]), mutation_scale=10,
#                                   color="blue",alpha=0.6) 
# Define known values
Trr = 76  # temperature in Kelvin
Prr = 21.1e5  # pressure in Pa (3.1 bar)

# Get surface tension at T=76K, using quality to define liquid/vapor boundary
surface_tension = CP.CoolProp.PropsSI("SURFACE_TENSION", "T", Trr, "Q", 0, "Nitrogen")

print(f"Surface tension of LN2 at T = {Trr} K: {surface_tension:.6f} N/m")

T_aa = 76  # kg/m^3
P_aa = 1.1e5  # Pa

# Find the rho using CoolProp
rho_aa = CP.CoolProp.PropsSI('D', 'T', T_aa, 'P', P_aa, 'N2')

print("T_aa = %d ; P_aa = %d ; rho_aa = %.3f" % (T_aa, P_aa, rho_aa))

#ax.add_patch(arrow_near1)
#ax.add_patch(arrow_near2)
#ax.add_patch(arrow_near3)
#ax.add_patch(arrow_near4)

#ax.add_patch(polygon)
#arrow_near_legend = mpatches.Patch(color='blue', alpha=0.6, label='Near critical')
#arrow_near_legend = mlines.Line2D([], [], color='blue', linestyle='-', linewidth=2, alpha=0.6, label='N2 spray')
######

######################################################################################################################
#show the legend
ax.legend(handles=[gtinj_legend, start_point_legend, end_point_legend], loc='lower right')

#ax.legend()
plt.tight_layout()
matplotlib.pyplot.savefig(
    'N2_with_Ranges_And_AllPoints.png',
    dpi=1000
)
# plt.show()
plt.close(fig)



print(f"Tc = {Tc} K, Pc = {pc} Pa") 