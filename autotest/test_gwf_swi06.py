"""
Test of the seawater intrusion (SWI) package for a sloping offshore ghb
boundary.

"""


import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "swi06a",
]

Lx = 10000.0  # x position of right edge of model domain
delr = 100.0
delc = 1.0
nlay = 1
nrow = 1
ncol = int(Lx / delr)
top_island = 50.0
botm_aquifer = -50.0
recharge = 0.0001
hydraulic_conductivity = 100.0
ghb_cond_fact = 1.0
perlen = 10 * 365.0
nper = 10
perioddata = nper * [(perlen, 1, 1.0)]
sealevel_start = 0.0
sealevel_stop = 40.0


def build_models(idx, test):
    ws = test.workspace
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=ws,
        exe_name="mf6",
        memory_print_option="all",
        print_input=True,
    )
    tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        no_ptcrecord=True,
        outer_maximum=500,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        newtonoptions=True,
    )
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top_island,
        botm=botm_aquifer,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        # alternative_cell_averaging="amt-hmk",  # "harmonic",
        icelltype=1,
        k=hydraulic_conductivity,
    )
    zeta_file = name + ".zta"
    swi = flopy.mf6.ModflowGwfswi(
        gwf,
        zeta_filerecord=zeta_file,
        tva6_filerecord=f"{name}.swi.tva",
    )

    # initialize the tva subpackage with the saltwater head
    sl = np.linspace(sealevel_start, sealevel_stop, nper)
    aux = {}
    for kper in range(nper):
        # this list has 2 entries, one for hsalt and one for myauxvar
        aux[kper] = [sl[kper], np.ones((nlay, nrow, ncol))]
    swi.tva.initialize(auxiliary=["saltwater_head", "myauxvar"], aux=aux)

    cghb = hydraulic_conductivity * ghb_cond_fact * delr * delc / 1
    ghbspd = {}
    for kper in range(nper):
        ghb_list = []
        for j in [0, ncol - 1]:
            freshwater_head = sl[kper]
            ghb_list.append([0, 0, j, freshwater_head, cghb])
        ghbspd[kper] = ghb_list

    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
    )
    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    budget_file = name + ".bud"
    head_file = name + ".hds"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        head_filerecord=head_file,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def plot_output(idx, test):
    import matplotlib.pyplot as plt

    sim = test.sims[0]
    ws = test.workspace
    gwf = sim.gwf[0]
    head = gwf.output.head().get_data().flatten()
    zeta = gwf.swi.output.zeta().get_data().flatten()
    volume_fresh = (head - zeta) * delr
    volume_fresh = volume_fresh.sum()
    pxs = flopy.plot.PlotCrossSection(gwf, line={"row": 0})
    ax = pxs.ax

    zobj = gwf.swi.output.zeta()
    times = zobj.times
    for t in times:
        zeta = zobj.get_data(totim=t)
        pxs.plot_surface(zeta)

    pxs.plot_grid()
    title = "Rising Sea Level"
    ax.set_title(title)
    ax.set_ylim(botm_aquifer, top_island)
    plt.savefig(ws / "zeta.png")
    plt.close("all")


def check_output(idx, test):
    # get the flopy sim object
    sim = test.sims[0]


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
