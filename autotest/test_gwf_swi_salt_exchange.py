"""
Salt-exchange term of the SWI formulation, single-fluid mode.

A one-layer strip of ten 100-ft cells: five land cells (top +10 ft) with
recharge, then five seafloor cells whose tops step down below sea level.
Three seafloor cells carry a GHB at the equivalent freshwater head of the sea
at the cell top (-top/alphaf with sea level at zero); the fourth has only a
drain, below that head; the fifth has no boundary at all. The lens discharges
through the GHB cells; the last two cells hold no freshwater, so in the
freshwater-only formulation their head is not set by any flow of freshwater.
Without the salt-exchange term those rows are held only by the minimum-zone
conductance floor; the term holds such a cell at its zone-vanishing head, the
head that puts the interface exactly at the top of the freshwater zone, and
reports the flow needed to do so as exchange with the unsimulated saltwater.

swi-saltexch
  the term is active (the default).
swi-saltexch-off
  DEV_NO_SALT_EXCHANGE: the term is off, no SALT-EXCHANGE budget entry is
  written, and the run still conserves mass.

Checks for swi-saltexch:
  * the run converges and conserves mass;
  * the two cells without freshwater sit at their zone-vanishing head
    (-top/alphaf), so their interface is at the cell top;
  * the SALT-EXCHANGE budget term equals the drain discharge from the cell
    that holds no freshwater (the only place the term can be nonzero at
    steady state), and is zero in every cell that holds freshwater;
  * the lens itself is unaffected: the recharge leaves through the GHB cells.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["swi-saltexch", "swi-saltexch-off"]

delr = 100.0
delc = 1.0
ncol = 10
land_top = 10.0
sea_tops = [-20.0, -40.0, -60.0, -80.0, -100.0]
botm = -150.0
k = 10.0
recharge = 0.001
alphaf = 40.0
drn_elev = 0.5
drn_cond = 0.01
ghb_cond = 100.0


def build_models(idx, test):
    ws = test.workspace
    name = "vzone"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        csv_outer_output_filerecord=f"{name}.outer.csv",
        outer_maximum=200,
        inner_maximum=300,
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-9,
        linear_acceleration="bicgstab",
        backtracking_number=20,
        backtracking_tolerance=1.05,
        backtracking_reduction_factor=0.5,
        backtracking_residual_limit=0.002,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    top = np.array(5 * [land_top] + sea_tops)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=k)
    flopy.mf6.ModflowGwfswi(
        gwf,
        zeta_filerecord=f"{name}.zta",
        saltwater_head=0.0,
        dev_no_salt_exchange=cases[idx].endswith("-off"),
    )
    # seafloor GHB cells at the equivalent freshwater head of the sea
    ghb = [[(0, 0, j), -top[j] / alphaf, ghb_cond] for j in (5, 6, 7)]
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghb)
    # a drain in a cell that holds no freshwater, below its zone-vanishing head
    flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=[[(0, 0, 8), drn_elev, drn_cond]])
    flopy.mf6.ModflowGwfrcha(gwf, recharge={0: 5 * [recharge] + 5 * [0.0]})
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.bud",
        budgetcsv_filerecord=f"{name}.bud.csv",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim, None


def check_output(idx, test):
    sim = test.sims[0]
    gwf = sim.gwf[0]
    ws = sim.sim_path
    name = gwf.name
    top = gwf.dis.top.array.flatten()
    head = gwf.output.head().get_data().flatten()
    zeta = gwf.swi.output.zeta().get_data().flatten()
    hpin = -top / alphaf

    outer = pd.read_csv(ws / f"{name}.outer.csv")
    print(f"outer iterations {len(outer)}")

    # mass balance
    df = pd.read_csv(ws / f"{name}.bud.csv")
    imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
    peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
    assert imbalance < 1.0e-8 + 1.0e-6 * peak, f"budget imbalance {imbalance:.3e}"

    cbc = flopy.utils.CellBudgetFile(ws / f"{name}.bud")
    if cases[idx].endswith("-off"):
        # the term is off: no budget entry anywhere
        assert not any("SALT-EXCHANGE" in c for c in df.columns), (
            "SALT-EXCHANGE in the budget with the term disabled"
        )
        names = [n.decode().strip() for n in cbc.get_unique_record_names()]
        assert "SALT-EXCHANGE" not in names, "SALT-EXCHANGE in the budget file"
        return

    assert len(outer) < 60, f"{len(outer)} outer iterations; the term is not holding"

    # the two cells without freshwater sit at their zone-vanishing head; the
    # drain cell can only be held from below, so it sits a hair under it
    # the term ramps over ZONE_PIN_BAND of the cell thickness (as head), and the
    # floor conductance to the drain cell leaves the head a few percent into it
    assert abs(head[9] - hpin[9]) < 1.0e-3, f"isolated cell head {head[9]} != {hpin[9]}"
    assert head[8] <= hpin[8] and hpin[8] - head[8] < 5.0e-3, (
        f"drain cell head {head[8]} not held at {hpin[8]}"
    )
    assert np.allclose(zeta[8:], top[8:], atol=1.0e-3), "interface not at the cell top"
    # the lens cells hold freshwater: head above the zone-vanishing head
    assert np.all(head[:5] > hpin[:5] + 1.0e-3), "lens cells lost their freshwater"
    assert np.all(head[5:8] >= hpin[5:8] - 1.0e-6), "GHB cells pulled below the sea"

    # the SALT-EXCHANGE term is the drain discharge from the cell without
    # freshwater and nothing else
    qpin = np.asarray(cbc.get_data(text="SALT-EXCHANGE")[-1]).flatten()
    qdrn = drn_cond * (head[8] - drn_elev)
    print(f"salt-exchange flow {qpin}, drain discharge {qdrn:.6f}")
    assert np.allclose(qpin[:8], 0.0, atol=1.0e-10), (
        "term active where freshwater exists"
    )
    # the floor conductances to the neighbours carry a few tenths of a percent
    assert abs(qpin[8] - qdrn) < 1.0e-2 * qdrn, "term does not balance the drain"
    assert abs(qpin[9]) < 1.0e-2 * qdrn, "term active in the cell with no boundary"
    vz_in = df["SALT-EXCHANGE(SWI)_IN"].iloc[-1]
    assert abs(vz_in - qdrn) < 1.0e-2 * qdrn, "budget term does not match"

    # the lens discharges through the GHB cells
    qrch = recharge * 5 * delr * delc
    qghb = df["GHB(GHB_0)_OUT"].iloc[-1] - df["GHB(GHB_0)_IN"].iloc[-1]
    assert abs(qghb - qrch) < 1.0e-2 * qrch, "recharge not leaving via GHB"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()
