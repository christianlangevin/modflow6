"""
Test for the Linear Network Flow (LNF) model.

This is a simple 3-node linear network test with:
- Node 1: CHD boundary at head = 1.0 (upstream)
- Node 2: Active cell (no boundary)
- Node 3: CHD boundary at head = 0.0 (downstream)
- Connections: 1-2 and 2-3, each with conductance = 1.0

Expected:
- Head at node 2 = 0.5 (midpoint between h1=1.0 and h3=0.0)
- Flow through network = C * dh = 1.0 * 0.5 = 0.5
- CHD at node 1: outflow = -0.5
- CHD at node 3: inflow = +0.5
"""

import flopy
import numpy as np
import pytest

from framework import TestFramework


cases = [
    "lnf_01",
]


def build_models(idx, test):
    name = "lnf01"
    sim_ws = test.workspace

    # Create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=sim_ws,
        exe_name="mf6",
        print_input=True,
    )

    # Create TDIS package
    flopy.mf6.ModflowTdis(
        sim,
        nper=1,
        perioddata=[(1.0, 1, 1.0)],
    )

    # Create IMS package
    flopy.mf6.ModflowIms(
        sim,
        complexity="SIMPLE",
        print_option="SUMMARY",
        outer_maximum=100,
        inner_maximum=100,
        outer_dvclose=1e-9,
        inner_dvclose=1e-12,
    )

    # Create LNF model
    lnf = flopy.mf6.ModflowLnf(
        sim,
        modelname=name,
        save_flows=True,
        print_flows=True,
    )

    # Create DISU package for 3 nodes with 2 connections
    # Network: 1 -- 2 -- 3
    #
    # Node connectivity (CSR format):
    # Node 1: connected to 1 (diag), 2
    # Node 2: connected to 1, 2 (diag), 3
    # Node 3: connected to 2, 3 (diag)
    #
    # IAC = [2, 3, 2] (node 1 has 2, node 2 has 3, node 3 has 2)
    # JA = [1, 2, 1, 2, 3, 2, 3] (1-based) or [0, 1, 0, 1, 2, 1, 2] (0-based)
    nodes = 3
    nja = 7
    iac = [2, 3, 2]
    ja = [
        0, 1, 
        1, 0, 2,
        2, 1
    ]  # 0-based node numbers

    flopy.mf6.ModflowLnfdisu(
        lnf,
        nodes=nodes,
        top=1.0,
        bot=0.0,
        area=1.0,
        nja=nja,
        iac=iac,
        ja=ja,
        ihc=1,
        cl12=1.0,
        hwva=1.0,
    )

    # Create CND package with conductance for each connection
    # Conductance array has nja entries (CSR format)
    # Position: 0=diag1, 1=1->2, 2=2->1, 3=diag2, 4=2->3, 5=3->2, 6=diag3
    # Diagonal entries are not used, off-diagonal = 1.0
    conductance = [0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0]

    flopy.mf6.ModflowLnfcnd(
        lnf,
        nja=nja,
        conductance=conductance,
    )

    # Create CHD package
    # Node 1: head = 1.0 (upstream)
    # Node 3: head = 0.0 (downstream)
    chd_spd = [
        (0, 1.0),  # node 1 (0-based), head = 1.0
        (2, 0.0),  # node 3 (0-based), head = 0.0
    ]

    flopy.mf6.ModflowLnfchd(
        lnf,
        stress_period_data=chd_spd,
    )

    # Create OC package
    flopy.mf6.ModflowLnfoc(
        lnf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    name = "lnf01"
    sim_ws = test.workspace

    # Read head file
    hds_file = sim_ws / f"{name}.hds"
    hds = flopy.utils.HeadFile(hds_file)
    head = hds.get_data()

    # Check heads
    # Node 1: CHD = 1.0
    # Node 2: Should be 0.5 (midpoint, since equal conductance on both sides)
    # Node 3: CHD = 0.0
    expected_head = np.array([1.0, 0.5, 0.0])
    np.testing.assert_allclose(
        head.flatten(),
        expected_head,
        rtol=1e-6,
        err_msg="Head values do not match expected values",
    )

    # Read budget file
    cbc_file = sim_ws / f"{name}.cbc"
    cbc = flopy.utils.CellBudgetFile(cbc_file)

    # Check CHD flows
    # Flow through network = C * dh = 1.0 * (1.0 - 0.5) = 0.5
    # CHD at node 1: outflow = -0.5 (water leaving CHD into network)
    # CHD at node 3: inflow = +0.5 (water entering CHD from network)
    chd_flow = cbc.get_data(text="CHD")[0]
    expected_chd = np.array([-0.5, 0.5])
    np.testing.assert_allclose(
        chd_flow["q"],
        expected_chd,
        rtol=1e-6,
        err_msg="CHD flow values do not match expected values",
    )


@pytest.mark.developmode
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
