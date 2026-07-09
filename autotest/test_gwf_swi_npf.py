"""
Reproduction test: with the interface (zeta) held well below the model bottom,
the SWI Package must reproduce a plain NPF simulation. This confirms that SWI is
a true no-op when there is no saltwater in the domain, and it exercises the
horizontal and vertical SWI conductance in the all-freshwater limit -- the
freshwater slab is the full cell (S^s = 0) and the vertical interface gate is
fully open (factor = 1).

The model is confined, multi-layer, with recharge on the top layer and a fixed
head on the bottom layer, so there is real vertical flow through the layer faces.
The single-fluid SWI package uses the default saltwater head of zero, so
zeta = -alphaf*hf is on the order of -300, far below the model bottom. The same
model is run with and without SWI and the heads must match.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["swi-npf"]

nlay, nrow, ncol = 3, 1, 5


def build_sim(ws, exe, with_swi):
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim)  # single steady-state period
    flopy.mf6.ModflowIms(
        sim,
        print_option="all",
        linear_acceleration="bicgstab",
        outer_maximum=100,
        inner_maximum=100,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=100.0,
        delc=1.0,
        top=10.0,
        botm=[5.0, 0.0, -5.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0, k33=1.0)
    # recharge on the top layer + fixed head on the bottom layer -> vertical flow
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001)
    chd = [[(nlay - 1, 0, c), 5.0] for c in range(ncol)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    if with_swi:
        # default saltwater head is zero -> zeta = -alphaf*hf << model bottom
        flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=name + ".zta")
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=name + ".hds", saverecord=[("HEAD", "ALL")]
    )
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, test.targets["mf6"], with_swi=True), None


def check_output(idx, test):
    # -- build and run the equivalent NPF-only (no SWI) model in a sibling dir
    npf_ws = test.workspace / "npf"
    sim_npf = build_sim(npf_ws, test.targets["mf6"], with_swi=False)
    sim_npf.write_simulation(silent=True)
    success, _ = sim_npf.run_simulation(silent=True)
    assert success, "NPF reference model failed to run"
    # -- heads must match: SWI is a no-op when zeta is below the model bottom
    h_swi = test.sims[0].gwf[0].output.head().get_data()
    h_npf = sim_npf.gwf[0].output.head().get_data()
    dmax = np.abs(h_swi - h_npf).max()
    print(f"max |h_swi - h_npf| = {dmax:.3e}")
    assert np.allclose(h_swi, h_npf, atol=1.0e-9, rtol=0.0), (
        f"SWI must reproduce NPF when zeta is below the model bottom "
        f"(max head difference {dmax:.3e})"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()
