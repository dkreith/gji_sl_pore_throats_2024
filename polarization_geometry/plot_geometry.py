import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 14})

# Pore radii
r1x_mesh = np.loadtxt('geometry/r1_r2/r1x_mesh.txt')
r2_mesh = np.loadtxt('geometry/r1_r2/r2_mesh.txt')

sig_r1_r2_sl = np.loadtxt('geometry/r1_r2/SL_sig_r1_r2.txt')
tau_r1_r2_sl = np.loadtxt('geometry/r1_r2/SL_tau_r1_r2.txt')

sig_r1_r2_mm = np.loadtxt('geometry/r1_r2/MM_sig_r1_r2.txt')
tau_r1_r2_mm = np.loadtxt('geometry/r1_r2/MM_tau_r1_r2.txt')

fig_r1_r2, [[ax_sig_r1_r2_sl, ax_tau_r1_r2_sl],
            [ax_sig_r1_r2_mm, ax_tau_r1_r2_mm]] = plt.subplots(2,2)
fig_r1_r2.set_figwidth(13)
fig_r1_r2.set_figheight(10)

div_r1_r2_sig_sl = make_axes_locatable(ax_sig_r1_r2_sl)
div_r1_r2_tau_sl = make_axes_locatable(ax_tau_r1_r2_sl)
cax_r1_r2_sig_sl = div_r1_r2_sig_sl.append_axes('right', size="5%", pad=0.1)
cax_r1_r2_tau_sl = div_r1_r2_tau_sl.append_axes('right', size="5%", pad=0.1)

r1_r2_sig_sl = ax_sig_r1_r2_sl.pcolormesh(r1x_mesh, r2_mesh,
                                          np.log10(sig_r1_r2_sl),
                                          shading='nearest', cmap='viridis',
                                          vmin=-14, vmax=0)
ax_sig_r1_r2_sl.set_xscale('log')
ax_sig_r1_r2_sl.set_yscale('log')
ax_sig_r1_r2_sl.set_xlabel('$r_1$ (m)')
ax_sig_r1_r2_sl.set_ylabel('$r_2$ (m)')
ax_sig_r1_r2_sl.set_title('(a)', loc='left')
cb_r1_r2_sig_sl = plt.colorbar(r1_r2_sig_sl, cax=cax_r1_r2_sig_sl)
cb_r1_r2_sig_sl.set_label('log($\sigma_{max}\'\'/\sigma_0$)')

r1_r2_tau_sl = ax_tau_r1_r2_sl.pcolormesh(r1x_mesh, r2_mesh,
                                          np.log10(tau_r1_r2_sl),
                                          shading='nearest', cmap='viridis',
                                          vmin=-5, vmax=-0.5)
ax_tau_r1_r2_sl.set_xscale('log')
ax_tau_r1_r2_sl.set_yscale('log')
ax_tau_r1_r2_sl.set_xlabel('$r_1$ (m)')
ax_tau_r1_r2_sl.set_ylabel('$r_2$ (m)')
ax_tau_r1_r2_sl.set_title('(b)', loc='left')
cb_r1_r2_tau_sl = plt.colorbar(r1_r2_tau_sl, cax=cax_r1_r2_tau_sl)
cb_r1_r2_tau_sl.set_label('log($\u03C4$)')

div_r1_r2_sig_mm = make_axes_locatable(ax_sig_r1_r2_mm)
div_r1_r2_tau_mm = make_axes_locatable(ax_tau_r1_r2_mm)
cax_r1_r2_sig_mm = div_r1_r2_sig_mm.append_axes('right', size="5%", pad=0.1)
cax_r1_r2_tau_mm = div_r1_r2_tau_mm.append_axes('right', size="5%", pad=0.1)

r1_r2_sig_mm = ax_sig_r1_r2_mm.pcolormesh(r1x_mesh, r2_mesh,
                                          np.log10(sig_r1_r2_mm),
                                          shading='nearest', cmap='viridis',
                                          vmin=-14, vmax=0)
ax_sig_r1_r2_mm.set_xscale('log')
ax_sig_r1_r2_mm.set_yscale('log')
ax_sig_r1_r2_mm.set_xlabel('$r_1$ (m)')
ax_sig_r1_r2_mm.set_ylabel('$r_2$ (m)')
ax_sig_r1_r2_mm.set_title('(c)', loc='left')
cb_r1_r2_sig_mm = plt.colorbar(r1_r2_sig_mm, cax=cax_r1_r2_sig_mm)
cb_r1_r2_sig_mm.set_label('log($\sigma_{max}\'\'/\sigma_0$)')

r1_r2_tau_mm = ax_tau_r1_r2_mm.pcolormesh(r1x_mesh, r2_mesh,
                                          np.log10(tau_r1_r2_mm),
                                          shading='nearest', cmap='viridis',
                                          vmin=-5, vmax=-0.5)
ax_tau_r1_r2_mm.set_xscale('log')
ax_tau_r1_r2_mm.set_yscale('log')
ax_tau_r1_r2_mm.set_xlabel('$r_1$ (m)')
ax_tau_r1_r2_mm.set_ylabel('$r_2$ (m)')
ax_tau_r1_r2_mm.set_title('(d)', loc='left')
cb_r1_r2_tau_mm = plt.colorbar(r1_r2_tau_mm, cax=cax_r1_r2_tau_mm)
cb_r1_r2_tau_mm.set_label('log($\u03C4$)')

ax_sig_r1_r2_sl.fill([1e-8, 1e-5, 1e-8, 1e-8], [1e-8, 1e-5, 1e-5, 1e-8],
                     color="silver")
ax_sig_r1_r2_mm.fill([1e-8, 1e-5, 1e-8, 1e-8], [1e-8, 1e-5, 1e-5, 1e-8],
                     color="silver")
ax_tau_r1_r2_sl.fill([1e-8, 1e-5, 1e-8, 1e-8], [1e-8, 1e-5, 1e-5, 1e-8],
                     color="silver")
ax_tau_r1_r2_mm.fill([1e-8, 1e-5, 1e-8, 1e-8], [1e-8, 1e-5, 1e-5, 1e-8],
                     color="silver")

fig_r1_r2.tight_layout()
fig_r1_r2.savefig("Figure_11.png",dpi=300,bbox_inches="tight")

# Pore lengths

L1_mesh = np.loadtxt('geometry/L1_L2/L1_mesh.txt')
L2_mesh = np.loadtxt('geometry/L1_L2/L2_mesh.txt')

sig_L1_L2_sl = np.loadtxt('geometry/L1_L2/SL_sig_L1_L2.txt')
tau_L1_L2_sl = np.loadtxt('geometry/L1_L2/SL_tau_L1_L2.txt')

sig_L1_L2_mm = np.loadtxt('geometry/L1_L2/MM_sig_L1_L2.txt')
tau_L1_L2_mm = np.loadtxt('geometry/L1_L2/MM_tau_L1_L2.txt')

fig_L1_L2, [[ax_sig_L1_L2_sl, ax_tau_L1_L2_sl],
            [ax_sig_L1_L2_mm, ax_tau_L1_L2_mm]] = plt.subplots(2,2)
fig_L1_L2.set_figwidth(13)
fig_L1_L2.set_figheight(10)

div_L1_L2_sig_sl = make_axes_locatable(ax_sig_L1_L2_sl)
div_L1_L2_tau_sl = make_axes_locatable(ax_tau_L1_L2_sl)
cax_L1_L2_sig_sl = div_L1_L2_sig_sl.append_axes('right', size="5%", pad=0.1)
cax_L1_L2_tau_sl = div_L1_L2_tau_sl.append_axes('right', size="5%", pad=0.1)

L1_L2_sig_sl = ax_sig_L1_L2_sl.pcolormesh(L1_mesh, L2_mesh,
                                          np.log10(sig_L1_L2_sl),
                                          shading='nearest', cmap='viridis',
                                          vmin=-6, vmax=-1)
ax_sig_L1_L2_sl.set_xscale('log')
ax_sig_L1_L2_sl.set_yscale('log')
ax_sig_L1_L2_sl.set_xlabel('$L_1$ (m)')
ax_sig_L1_L2_sl.set_ylabel('$L_2$ (m)')
ax_sig_L1_L2_sl.set_title('(a)', loc='left')
cb_L1_L2_sig_sl = plt.colorbar(L1_L2_sig_sl, cax=cax_L1_L2_sig_sl)
cb_L1_L2_sig_sl.set_label('log($\sigma_{max}\'\'/\sigma_0$)')

L1_L2_tau_sl = ax_tau_L1_L2_sl.pcolormesh(L1_mesh, L2_mesh,
                                          np.log10(tau_L1_L2_sl),
                                          shading='nearest', cmap='viridis',
                                          vmin=-5.5, vmax=2.5)
ax_tau_L1_L2_sl.set_xscale('log')
ax_tau_L1_L2_sl.set_yscale('log')
ax_tau_L1_L2_sl.set_xlabel('$L_1$ (m)')
ax_tau_L1_L2_sl.set_ylabel('$L_2$ (m)')
ax_tau_L1_L2_sl.set_title('(b)', loc='left')
cb_L1_L2_tau_sl = plt.colorbar(L1_L2_tau_sl, cax=cax_L1_L2_tau_sl)
cb_L1_L2_tau_sl.set_label('log($\u03C4$)')

div_L1_L2_sig_mm = make_axes_locatable(ax_sig_L1_L2_mm)
div_L1_L2_tau_mm = make_axes_locatable(ax_tau_L1_L2_mm)
cax_L1_L2_sig_mm = div_L1_L2_sig_mm.append_axes('right', size="5%", pad=0.1)
cax_L1_L2_tau_mm = div_L1_L2_tau_mm.append_axes('right', size="5%", pad=0.1)

L1_L2_sig_mm = ax_sig_L1_L2_mm.pcolormesh(L1_mesh, L2_mesh,
                                          np.log10(sig_L1_L2_mm),
                                          shading='nearest', cmap='viridis',
                                          vmin=-6, vmax=-1)
ax_sig_L1_L2_mm.set_xscale('log')
ax_sig_L1_L2_mm.set_yscale('log')
ax_sig_L1_L2_mm.set_xlabel('$L_1$ (m)')
ax_sig_L1_L2_mm.set_ylabel('$L_2$ (m)')
ax_sig_L1_L2_mm.set_title('(c)', loc='left')
cb_L1_L2_sig_mm = plt.colorbar(L1_L2_sig_mm, cax=cax_L1_L2_sig_mm)
cb_L1_L2_sig_mm.set_label('log($\sigma_{max}\'\'/\sigma_0$)')

L1_L2_tau_mm = ax_tau_L1_L2_mm.pcolormesh(L1_mesh, L2_mesh,
                                          np.log10(tau_L1_L2_mm),
                                          shading='nearest', cmap='viridis',
                                          vmin=-5.5, vmax=2.5)
ax_tau_L1_L2_mm.set_xscale('log')
ax_tau_L1_L2_mm.set_yscale('log')
ax_tau_L1_L2_mm.set_xlabel('$L_1$ (m)')
ax_tau_L1_L2_mm.set_ylabel('$L_2$ (m)')
ax_tau_L1_L2_mm.set_title('(d)', loc='left')
cb_L1_L2_tau_mm = plt.colorbar(L1_L2_tau_mm, cax=cax_L1_L2_tau_mm)
cb_L1_L2_tau_mm.set_label('log($\u03C4$)')

fig_L1_L2.tight_layout()
fig_L1_L2.savefig("Figure_A01.png",dpi=300,bbox_inches="tight")

# Pore length and radius of wide pore

L1_mesh_lr = np.loadtxt('geometry/L1_r1/L1_mesh.txt')
r1y_mesh = np.loadtxt('geometry/L1_r1/r1y_mesh.txt')

sig_L1_r1_sl = np.loadtxt('geometry/L1_r1/SL_sig_L1_r1.txt')
tau_L1_r1_sl = np.loadtxt('geometry/L1_r1/SL_tau_L1_r1.txt')

sig_L1_r1_mm = np.loadtxt('geometry/L1_r1/MM_sig_L1_r1.txt')
tau_L1_r1_mm = np.loadtxt('geometry/L1_r1/MM_tau_L1_r1.txt')

fig_L1_r1, [[ax_sig_L1_r1_sl, ax_tau_L1_r1_sl],
            [ax_sig_L1_r1_mm, ax_tau_L1_r1_mm]] = plt.subplots(2,2)
fig_L1_r1.set_figwidth(13)
fig_L1_r1.set_figheight(10)

div_L1_r1_sig_sl = make_axes_locatable(ax_sig_L1_r1_sl)
div_L1_r1_tau_sl = make_axes_locatable(ax_tau_L1_r1_sl)
cax_L1_r1_sig_sl = div_L1_r1_sig_sl.append_axes('right', size="5%", pad=0.1)
cax_L1_r1_tau_sl = div_L1_r1_tau_sl.append_axes('right', size="5%", pad=0.1)

L1_r1_sig_sl = ax_sig_L1_r1_sl.pcolormesh(L1_mesh_lr, r1y_mesh,
                                          np.log10(sig_L1_r1_sl),
                                          shading='nearest', cmap='viridis',
                                          vmin=-16, vmax=0)
ax_sig_L1_r1_sl.set_xscale('log')
ax_sig_L1_r1_sl.set_yscale('log')
ax_sig_L1_r1_sl.set_xlabel('$L_1$ (m)')
ax_sig_L1_r1_sl.set_ylabel('$r_1$ (m)')
ax_sig_L1_r1_sl.set_title('(a)', loc='left')
cb_L1_r1_sig_sl = plt.colorbar(L1_r1_sig_sl, cax=cax_L1_r1_sig_sl)
cb_L1_r1_sig_sl.set_label('log($\sigma_{max}\'\'/\sigma_0$)')

L1_r1_tau_sl = ax_tau_L1_r1_sl.pcolormesh(L1_mesh_lr, r1y_mesh,
                                          np.log10(tau_L1_r1_sl),
                                          shading='nearest', cmap='viridis',
                                          vmin=-6, vmax=3)
ax_tau_L1_r1_sl.set_xscale('log')
ax_tau_L1_r1_sl.set_yscale('log')
ax_tau_L1_r1_sl.set_xlabel('$L_1$ (m)')
ax_tau_L1_r1_sl.set_ylabel('$r_1$ (m)')
ax_tau_L1_r1_sl.set_title('(b)', loc='left')
cb_L1_r1_tau_sl = plt.colorbar(L1_r1_tau_sl, cax=cax_L1_r1_tau_sl)
cb_L1_r1_tau_sl.set_label('log($\u03C4$)')

div_L1_r1_sig_mm = make_axes_locatable(ax_sig_L1_r1_mm)
div_L1_r1_tau_mm = make_axes_locatable(ax_tau_L1_r1_mm)
cax_L1_r1_sig_mm = div_L1_r1_sig_mm.append_axes('right', size="5%", pad=0.1)
cax_L1_r1_tau_mm = div_L1_r1_tau_mm.append_axes('right', size="5%", pad=0.1)

L1_r1_sig_mm = ax_sig_L1_r1_mm.pcolormesh(L1_mesh_lr, r1y_mesh,
                                          np.log10(sig_L1_r1_mm),
                                          shading='nearest', cmap='viridis',
                                          vmin=-16, vmax=0)
ax_sig_L1_r1_mm.set_xscale('log')
ax_sig_L1_r1_mm.set_yscale('log')
ax_sig_L1_r1_mm.set_xlabel('$L_1$ (m)')
ax_sig_L1_r1_mm.set_ylabel('$r_1$ (m)')
ax_sig_L1_r1_mm.set_title('(c)', loc='left')
cb_L1_r1_sig_mm = plt.colorbar(L1_r1_sig_mm, cax=cax_L1_r1_sig_mm)
cb_L1_r1_sig_mm.set_label('log($\sigma_{max}\'\'/\sigma_0$)')

L1_r1_tau_mm = ax_tau_L1_r1_mm.pcolormesh(L1_mesh_lr, r1y_mesh,
                                          np.log10(tau_L1_r1_mm),
                                          shading='nearest', cmap='viridis',
                                          vmin=-6, vmax=3)
ax_tau_L1_r1_mm.set_xscale('log')
ax_tau_L1_r1_mm.set_yscale('log')
ax_tau_L1_r1_mm.set_xlabel('$L_1$ (m)')
ax_tau_L1_r1_mm.set_ylabel('$r_1$ (m)')
ax_tau_L1_r1_mm.set_title('(d)', loc='left')
cb_L1_r1_tau_mm = plt.colorbar(L1_r1_tau_mm, cax=cax_L1_r1_tau_mm)
cb_L1_r1_tau_mm.set_label('log($\u03C4$)')

fig_L1_r1.tight_layout()
fig_L1_r1.savefig("Figure_A02.png",dpi=300,bbox_inches="tight")