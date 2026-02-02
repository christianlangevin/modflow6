"""
Test of the seawater intrusion (SWI) package for a
simple 2-cell transient model.  Cell 1 has a GHB boundary
condition representing sea level.  Cell 2 starts with
head at 2.0 m above sea level.  After one time step,
the freshwater volume change in the aquifer should
equal the flow through the GHB boundary.

"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "swi07a",
]

newtonoptions = ["newton"]
icelltype = 1
iconvert = 0

# Lx = 10000  # meters
delr, delc = 100.0, 1.0
ncol = 2
nlay = 1
nrow = 1
top = [0.0, 0.0]
botm = -100.0
hydraulic_conductivity = 1.0
specific_yield = 0.2
specific_storage = 0.0  # 1.0e-3
sea_level = 0.0
h0 = (
    sea_level - top[0]
) * 1.025  # sea level converted to freshwater head at top of first cell
strt = [h0, 2.0]
perlen = 100.0
nstp = 1
tsmult = 1.0
perioddata = [(perlen, nstp, tsmult)]


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
    tdis = flopy.mf6.ModflowTdis(sim, perioddata=perioddata)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="all",
        outer_maximum=500,
        linear_acceleration="bicgstab",
        inner_dvclose=1e-8,
        outer_dvclose=1e-8,
        rcloserecord=1e-8,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        newtonoptions=newtonoptions[idx],
    )
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        alternative_cell_averaging=None,
        icelltype=icelltype,
        k=hydraulic_conductivity,
    )
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=iconvert,
        sy=specific_yield,
        ss=specific_storage,
        transient={0: True},
    )
    zeta_file = name + ".zta"
    swi = flopy.mf6.ModflowGwfswi(
        gwf,
        zeta_filerecord=zeta_file,
        saltwater_head=0.0,
    )
    cghb = hydraulic_conductivity / 10.0 * delr * delc / 1.0
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=[[0, 0, 0, h0, cghb]],
    )
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


def make_cross_section_plot(sim, idx, title):
    import matplotlib.pyplot as plt

    ws = pl.Path(sim.sim_path)
    gwf = sim.gwf[0]
    x = gwf.modelgrid.xcellcenters.flatten()
    botm = gwf.dis.botm.array.flatten()
    ws = sim.sim_path
    fpth = ws / "mymodel.zta"
    head = gwf.output.head().get_data().flatten()
    zeta = flopy.utils.HeadFile(fpth, text="zeta").get_data().flatten()
    tp = top
    pxs = flopy.plot.PlotCrossSection(gwf, line={"row": 0})
    ax = pxs.ax
    ax.plot(x, head, "k-")
    ax.plot(x, zeta, "k--")
    ax.fill_between(x, tp, zeta, color="cyan")
    ax.fill_between(x, zeta, botm, color="red")
    ax.set_title(title)
    ax.set_ylim(-400, 50)
    plt.savefig(ws / "zeta.png")
    plt.close("all")


def plot_output(idx, test):
    title = {
        0: "Case 7a",
    }
    sim = test.sims[0]
    make_cross_section_plot(sim, idx, title[idx])


def get_freshwater_volume(delr, delc, specific_yield, top, zeta):
    volume_fresh = 0.0
    for i in range(len(zeta)):
        dz_fresh = top[i] - zeta[i]
        volume_fresh += dz_fresh
    return volume_fresh * delr * delc * specific_yield


def check_output(idx, test):
    # get the flopy sim object
    sim = test.sims[0]
    gwf = sim.gwf[0]
    ws = sim.sim_path
    head = gwf.output.head().get_data().flatten()
    zeta = gwf.swi.output.zeta().get_data().flatten()
    ghb_flow = gwf.output.budget().get_data(text="GHB")[0]["q"][0]
    swi_storage_flow = gwf.output.budget().get_data(text="STORAGE")[0].flatten()
    print(f"swi_storage_flow = {swi_storage_flow}")
    print(f"head = {head}")
    print(f"zeta = {zeta}")

    # calculate calculate freshwater volumes
    zeta_start = [-40.0 * h for h in strt]
    print(f"zeta_start = {zeta_start}")
    volume_fresh_start = get_freshwater_volume(
        delr, delc, specific_yield, top, zeta_start
    )
    volume_fresh_end = get_freshwater_volume(
        delr, delc, specific_yield, top, zeta.flatten()
    )
    qout_calc = (volume_fresh_end - volume_fresh_start) / perlen

    qout_calc_alt = zeta.flatten() - np.array(zeta_start)
    qout_calc_alt = qout_calc_alt.sum() * delr * delc * specific_yield / perlen
    print(f"qout_calc_alt = {qout_calc_alt}")

    print(f"volume_fresh_start = {volume_fresh_start}")
    print(f"volume_fresh_end   = {volume_fresh_end}")
    print(f"qout_calc          = {qout_calc}")
    print(f"ghb_flow           = {ghb_flow}")
    print(f"qout - ghb_flow    = {qout_calc - ghb_flow}")
    print(f"swi_storage_flow   = {swi_storage_flow.sum()}")
    # assert np.isclose(qout, ghb_flow), f"qout {qout} != ghb_flow {ghb_flow}"
    # assert np.isclose(qout_calc, swi_storage_flow.sum()), (
    #     f"qout_calc {qout_calc} != swi_storage_flow {swi_storage_flow.sum()}"
    # )
    assert np.isclose(ghb_flow, -swi_storage_flow.sum()), (
        f"ghb_flow {ghb_flow} != swi_storage_flow {swi_storage_flow.sum()}"
    )
    assert np.allclose(-head * 40, zeta), f"head {-head * 40} != zeta {zeta}"


@pytest.mark.developmode
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
