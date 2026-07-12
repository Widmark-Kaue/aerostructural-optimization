#%% libs
import numpy as np
import matplotlib.pyplot as plt
from pickle import load
from pathlib import Path

from asalib import ASAOptimization
from utils import set_aiaa_style

#%% Plot settings
set_aiaa_style(16)
dpi = 600
format = 'pdf'
saveflag = True
#%% Path settings
rootdir = Path('.')
if not rootdir.absolute().name == 'exam':
    rootdir = Path('.', 'exam')

datadir = rootdir.joinpath('data')
datadir.mkdir(exist_ok=True, parents=True)
print("Data will be load in:", datadir)

imagdir = rootdir.joinpath('images')
imagdir.mkdir(exist_ok=True)
print("Images will be saved in:", imagdir)

#%% Load results
cases = [
    {'obj': 'FB', 'AR': 6.0, 'file': 'FB_results_AR6.0.pkl'},
    {'obj': 'FB', 'AR': 10.0, 'file': 'FB_results_AR10.0.pkl'},
    {'obj': 'Weight', 'AR': 6.0, 'file': 'Weight_results_AR6.0.pkl'},
    {'obj': 'Weight', 'AR': 10.0, 'file': 'Weight_results_AR10.0.pkl'}
]

loaded_data = []

for case in cases:
    filepath = datadir / case['file']
    
    assert filepath.exists(), 'File not exist'
    
    with open(filepath, 'rb') as f:
        res = load(f)
    
    # Extract the ASA object 
    asa:ASAOptimization = res.pop('asa')
    desVarsOpt = res.x
    
    # Solve physics for the optimized design variables
    out, stateVars = asa.solve_asa(desVarsOpt)
    
    loaded_data.append({
        'obj': case['obj'],
        'AR': case['AR'],
        'res': res,
        'asa': asa,
        'desVarsOpt': desVarsOpt,
        'out': out,
        'stateVars': stateVars
    })

npanels = loaded_data[0]['asa'].npanels

#%% Print comparison table
print("\n" + "="*80)
print(f"{'COMPARISON OF OPTIMIZED DESIGNS':^80}")
print("="*80)
print(f"{'Case':<20} | {'Obj Value':<12} | {'Fuel Burn (N)':<14} | {'Ws (N)':<10} | {'W0 (N)':<10} |")
print("-"*80)
for d in loaded_data:
    obj_val = d['res'].fun
    fb = d['out']['fb']
    w0 = d['out']['weight']
    
    # Compute structural weight 
    ws = w0 - fb - d['asa'].fixedMass
    label = f"{d['obj']} (AR={d['AR']})"
    
    print(f"{label:<20} | {obj_val:<12.4f} | {fb:<14.4f} | {ws:<10.4f} | {w0:<10.4f} |")
print("="*80 + "\n")

#%% 5.8 Compare Design Variables
plt.close('all')
plt.figure(figsize=(12, 4.5))

colors = {
    ('FB', 6.0): 'b',
    ('FB', 10.0): 'c',
    ('Weight', 6.0): 'r',
    ('Weight', 10.0): 'darkorange'
}

markers = {
    ('FB', 6.0): 'o-',
    ('FB', 10.0): 's-',
    ('Weight', 6.0): '^-',
    ('Weight', 10.0): 'd-'
}

# Subplot 1: Twist
plt.subplot(1, 2, 1)
for d in loaded_data:
    asa = d['asa']
    b_wing = asa.span
    y_norm = asa.ypanels / b_wing
    twist_opt = np.rad2deg(d['desVarsOpt'][:npanels])
    lbl = f"{d['obj']} (AR={d['AR']})"
    key = (d['obj'], d['AR'])
    plt.plot(y_norm, twist_opt, markers[key], color=colors[key], label=lbl)

plt.title('(a)')
plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$\tau$ [deg]')
plt.xlim(-0.5, 0.5)
plt.grid(True)
plt.legend()

# Subplot 2: Thickness
plt.subplot(1, 2, 2)
for d in loaded_data:
    asa = d['asa']
    b_wing = asa.span
    y_norm = asa.ypanels / b_wing
    t_opt = d['desVarsOpt'][npanels:] / 100 * 1000 # Convert to mm
    lbl = f"{d['obj']} (AR={d['AR']})"
    key = (d['obj'], d['AR'])
    plt.plot(y_norm, t_opt, markers[key], color=colors[key], label=lbl)

plt.title('(b)')
plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$t$ [mm]')
plt.xlim(-0.5, 0.5)
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig(imagdir / f'q5.8_1_comp_desVars.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% 5.8 Compare Circulation Distributions
plt.figure(figsize=(8, 5))

# Plot the normalized elliptical reference distribution
y_dense = np.linspace(-0.5, 0.5, 200)
gamma_elliptical_dense_norm = 2.0 * np.sqrt((0.5 + y_dense) * (0.5 - y_dense))
plt.plot(y_dense, gamma_elliptical_dense_norm, 'k-', label='Elíptica')

# Plot the markers at the panel positions 
y_panels_norm = loaded_data[0]['asa'].ypanels / loaded_data[0]['asa'].span
gamma_elliptical_panels_norm = 2.0 * np.sqrt((0.5 + y_panels_norm) * (0.5 - y_panels_norm))
plt.plot(y_panels_norm, gamma_elliptical_panels_norm, 'ko')

for d in loaded_data:
    asa = d['asa']
    span = asa.span
    vinf = asa.Vinfm
    y_norm = asa.ypanels / span
    
    gamma_opt = d['stateVars'][:npanels]
    cl_val = d['out']['cl']
    
    # Analytical Elliptical lift distribution scaling parameter
    gamma_0 = 2 * vinf * 16.0 / (np.pi * span) * cl_val
    
    # Normalize optimized circulation
    gamma_opt_norm = gamma_opt / gamma_0
    
    lbl = f"{d['obj']} (AR={d['AR']})"
    key = (d['obj'], d['AR'])
    color = colors[key]
    
    plt.plot(y_norm, gamma_opt_norm, markers[key], color=color, label=lbl)

plt.xlabel(r'$y/b$ [-]')
plt.ylabel(r'$\Gamma / \Gamma_0$ [-]')
plt.xlim(-0.5, 0.5)
plt.ylim(-0.2, 1.8)
plt.grid(True)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(imagdir / f'q5.8_2_comp_lift.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% 5.8 Weight Fractions 
plt.figure(figsize=(8, 5))

cases_labels = [f"{d['obj']}\n(AR={d['AR']})" for d in loaded_data]
w_f_fracs = []
w_s_fracs = []
w_fb_fracs = []

for d in loaded_data:
    asa = d['asa']
    w0 = d['out']['weight']
    fb = d['out']['fb']
    ws = w0 - fb - asa.fixedMass
    
    w_f_fracs.append(asa.fixedMass / w0 * 100)
    w_s_fracs.append(ws / w0 * 100)
    w_fb_fracs.append(fb / w0 * 100)

w_f_fracs = np.array(w_f_fracs)
w_s_fracs = np.array(w_s_fracs)
w_fb_fracs = np.array(w_fb_fracs)

# Stacked bar plot
p1 = plt.bar(cases_labels, w_f_fracs, color='b', edgecolor='k', width=0.4, label=r'$W_f$')
p2 = plt.bar(cases_labels, w_s_fracs, bottom=w_f_fracs, color='r', edgecolor='k', width=0.4, label=r'$W_s$')
p3 = plt.bar(cases_labels, w_fb_fracs, bottom=w_f_fracs + w_s_fracs, color='g', edgecolor='k', width=0.4, label=r'$W_{FB}$')

# Add labels inside the bars
for idx in range(len(cases_labels)):
    plt.text(idx, w_f_fracs[idx]/2, f"{w_f_fracs[idx]:.1f}%", ha='center', va='center', color='white', fontweight='bold')
    plt.text(idx, w_f_fracs[idx] + w_s_fracs[idx]/2, f"{w_s_fracs[idx]:.1f}%", ha='center', va='center', color='white', fontweight='bold')
    plt.text(idx, w_f_fracs[idx] + w_s_fracs[idx] + w_fb_fracs[idx]/2, f"{w_fb_fracs[idx]:.1f}%", ha='center', va='center', color='white', fontweight='bold')

plt.ylabel('Fração de Peso [%]')
plt.ylim(0, 115)
plt.grid(axis='y', linestyle=':', alpha=0.5)
plt.legend(ncol=3, loc='upper center')
plt.tight_layout()
plt.savefig(imagdir / f'q5.8_3_weight_fractions.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
# plt.show()

#%% 5.8 Induced Drag Relative Error (Optimized vs Elliptical)
plt.figure(figsize=(8, 4))

drag_errors = []

for d in loaded_data:
    asa = d['asa']
    w0 = d['out']['weight']
    fb = d['out']['fb']
    cl_val = d['out']['cl']
    
    # Calculate L/D from fuel burn formula
    l_over_d = -(asa.endurance * asa.TSFC) / np.log(1 - fb/w0)
    
    # Calculate total CD and induced CDi
    cd_val = cl_val / l_over_d
    cdi_val = cd_val - asa.cd0
    
    # Analytical elliptical drag coefficient
    cd_an = (cl_val**2) * 16.0 / (np.pi * (asa.span**2))
    
    # Relative error in percentage
    rel_error = (cdi_val - cd_an) / cd_an * 100.0
    drag_errors.append(rel_error)

drag_errors = np.array(drag_errors)

# Plot bar chart
bars = plt.bar(cases_labels, drag_errors, color='orange', edgecolor='k', width=0.4)

# Add line at y=0
plt.axhline(0, color='black', linewidth=1.2)

# Label bars with appropriate sign and positioning
for bar in bars:
    height = bar.get_height()
    if height >= 0:
        plt.text(bar.get_x() + bar.get_width()/2.0, height + 3.0, f"+{height:.1f}%", ha='center', va='bottom', fontweight='bold')
    else:
        plt.text(bar.get_x() + bar.get_width()/2.0, height - 3.0, f"{height:.1f}%", ha='center', va='top', fontweight='bold')

plt.ylabel(r'$(C_{Di} - C_{D,\text{eli}}) / C_{D,\text{eli}}$')

# Adjust y limits to leave space for labels
ymin, ymax = min(drag_errors), max(drag_errors)
plt.ylim(ymin - 15.0 if ymin < 0 else -10.0, ymax + 20.0)

plt.grid(axis='y', linestyle=':', alpha=0.5)
plt.tight_layout()
plt.savefig(imagdir / f'q5.8_4_comp_drag_ratio.{format}', dpi=dpi, bbox_inches='tight') if saveflag else None
plt.show()
