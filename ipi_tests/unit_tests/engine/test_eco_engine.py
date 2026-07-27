"""Engine-level tests for economised (eco) path integrals: input parsing,
binding, the omegak depend chain, spring/estimator consistency, propagation
and the guards against incompatible options."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pytest
import numpy as np

from ipi.engine.beads import Beads
from ipi.engine.cell import Cell
from ipi.inputs.beads import InputBeads
from ipi.inputs.cell import InputCell
from ipi.inputs.normalmodes import InputNormalModes
from ipi.scripting import InteractiveSimulation
from ipi.utils import nmtransform
from ipi.utils.depend import dstrip
from ipi.utils.io.inputs.io_xml import xml_parse_string
from ipi.utils.units import Constants

NBEADS = 8
NATOMS = 2


def eco_xml(
    frequencies="<frequencies style='eco' units='inversecm'> [2500] </frequencies>",
    propagator="exact",
    dynamics="nvt",
    thermostat="<thermostat mode='pile_l'> <tau units='femtosecond'> 10 </tau> </thermostat>",
    nm_extra="",
    motion=None,
):
    """Builds the XML input for a small harmonic-potential PIMD simulation."""

    if motion is None:
        motion = f"""<motion mode='dynamics'>
      <dynamics mode='{dynamics}'>
        <timestep units='femtosecond'> 0.25 </timestep>
        {thermostat}
      </dynamics>
    </motion>"""

    rng = np.random.RandomState(27182)
    beads = Beads(nbeads=NBEADS, natoms=NATOMS)
    # spreads the beads so that the internal modes are excited
    beads.q = np.array([0.0, 0.0, 0.0, 2.0, 0.0, 0.0] * NBEADS).reshape(
        NBEADS, 3 * NATOMS
    ) + 0.1 * rng.normal(size=(NBEADS, 3 * NATOMS))
    beads.m = np.full(NATOMS, 1837.36)
    beads.names = np.full(NATOMS, "H")
    input_beads = InputBeads()
    input_beads.store(beads)
    input_cell = InputCell()
    input_cell.store(Cell(h=20.0 * np.eye(3)))

    return f"""
<simulation verbosity='quiet'>
  <output prefix='sim'>
    <properties stride='100'> [ step, conserved ] </properties>
  </output>
  <total_steps> 100 </total_steps>
  <prng> <seed> 31415 </seed> </prng>
  <ffdirect name='harm'> <pes> harmonic </pes> <parameters> {{k1: 0.01}} </parameters> </ffdirect>
  <system>
    {input_beads.write("beads")}
    {input_cell.write("cell")}
    <initialize nbeads='{NBEADS}'>
      <velocities mode='thermal' units='kelvin'> 300 </velocities>
    </initialize>
    <forces> <force forcefield='harm'/> </forces>
    <ensemble> <temperature units='kelvin'> 300 </temperature> </ensemble>
    {motion}
    <normal_modes propagator='{propagator}'>
      {frequencies}
      {nm_extra}
    </normal_modes>
  </system>
</simulation>
"""


@pytest.fixture
def eco_sim(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return InteractiveSimulation(eco_xml())


def test_omegak_and_masses(eco_sim):
    """Checks the eco frequencies against the fit, and the dynamical masses."""

    nm = eco_sim.syslist[0].nm
    omegak = dstrip(nm.omegak)
    xmax = float(dstrip(nm.nm_freqs)[0] * nm.nbeads / nm.omegan)
    np.testing.assert_allclose(
        omegak, nm.omegan * nmtransform.eco_eva(nm.nbeads, xmax), rtol=1e-10
    )
    assert omegak[0] == 0.0
    # eco changes the springs but not the dynamical masses
    np.testing.assert_allclose(dstrip(nm.nm_factor), np.ones(NBEADS))


def test_spring_and_kinetic_td(eco_sim):
    """Checks the spring energy against an explicit normal-mode sum, and the
    consistency of the primitive kinetic energy estimator with it."""

    sys = eco_sim.syslist[0]
    nm = sys.nm
    qnm = dstrip(nm.qnm)
    m3 = dstrip(sys.beads.m3)
    wk2 = dstrip(nm.omegak) ** 2
    vspring = 0.5 * np.sum(wk2[:, np.newaxis] * m3 * qnm**2)
    np.testing.assert_allclose(dstrip(nm.vspring), vspring, rtol=1e-12)

    ktd = sys.properties["kinetic_td"][0]
    kt = Constants.kb * sys.ensemble.temp
    np.testing.assert_allclose(
        ktd, 1.5 * NATOMS * NBEADS * kt - vspring / NBEADS, rtol=1e-10
    )


def test_temperature_refit(eco_sim):
    """Checks that a temperature change refits the frequencies exactly
    through the depend chain (exercising the warm-start path)."""

    sys = eco_sim.syslist[0]
    nm = sys.nm
    gamma_old = dstrip(nm.omegak).copy() / nm.omegan
    sys.ensemble.temp = sys.ensemble.temp * 1.13
    omegak = dstrip(nm.omegak)
    xmax = float(dstrip(nm.nm_freqs)[0] * nm.nbeads / nm.omegan)
    np.testing.assert_allclose(
        omegak, nm.omegan * nmtransform.eco_eva(nm.nbeads, xmax), rtol=1e-8
    )
    # the dimensionless spectrum must change shape, not just rescale
    gamma_new = omegak / nm.omegan
    assert not np.allclose(gamma_new[1:], gamma_old[1:], rtol=1e-6)


@pytest.mark.parametrize("propagator", ["exact", "cayley", "bab"])
def test_eco_dynamics(tmp_path, monkeypatch, propagator):
    """Runs a few steps of eco PIMD with each free-ring-polymer propagator."""

    monkeypatch.chdir(tmp_path)
    sim = InteractiveSimulation(eco_xml(propagator=propagator))
    e0 = sim.properties("conserved")
    sim.run(4, write_outputs=False)
    assert np.isfinite(sim.properties("conserved"))
    assert abs(sim.properties("conserved") - e0) < 1e-2 * abs(e0)


def test_guard_bosons(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="bosons or open paths"):
        InteractiveSimulation(eco_xml(nm_extra="<bosons id='index'> [0] </bosons>"))


def test_guard_open_paths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="bosons or open paths"):
        InteractiveSimulation(eco_xml(nm_extra="<open_paths> [0] </open_paths>"))


def test_guard_nm_freqs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="requires one frequency"):
        InteractiveSimulation(
            eco_xml(
                frequencies="<frequencies style='eco' units='inversecm'> [2500, 3000] </frequencies>"
            )
        )


def test_guard_suzuki_chin(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="Suzuki-Chin"):
        InteractiveSimulation(eco_xml(dynamics="sc"))


def test_guard_nm_gle(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    a_matrix = str([1e-3] * NBEADS)
    thermostat = f"<thermostat mode='nm_gle'> <A shape='({NBEADS},1,1)'> {a_matrix} </A> </thermostat>"
    with pytest.raises(ValueError, match="nm_gle"):
        InteractiveSimulation(eco_xml(thermostat=thermostat))


def test_guard_instanton(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    motion = """<motion mode='instanton'>
      <instanton mode='rate'> <opt> nichols </opt> </instanton>
    </motion>"""
    with pytest.raises(ValueError, match="Instanton"):
        InteractiveSimulation(eco_xml(motion=motion))


@pytest.mark.parametrize(
    "prop",
    [
        "isotope_zetatd(alpha=1.1;atom=H)",
        "chin_weight",
        "ti_weight",
        "kinetic_prsc",
    ],
)
def test_guard_estimators(eco_sim, prop):
    """Estimators that assume Trotter springs must refuse to run under eco."""

    with pytest.raises(ValueError, match="Trotter"):
        eco_sim.syslist[0].properties[prop]


def test_trotter_estimators_unaffected(tmp_path, monkeypatch):
    """The same estimators keep working with Trotter springs."""

    monkeypatch.chdir(tmp_path)
    sim = InteractiveSimulation(
        eco_xml(frequencies="<frequencies style='rpmd'> [] </frequencies>")
    )
    value = sim.syslist[0].properties["isotope_zetatd(alpha=1.1;atom=H)"][0]
    assert np.all(np.isfinite(value))


def test_input_roundtrip(eco_sim):
    """Checks that mode and frequency survive an input write/parse cycle."""

    nm = eco_sim.syslist[0].nm
    inm = InputNormalModes()
    inm.store(nm)
    reparsed = InputNormalModes()
    reparsed.parse(xml_parse_string(inm.write("normal_modes")).fields[0][1])
    nm2 = reparsed.fetch()
    assert nm2.mode == "eco"
    np.testing.assert_allclose(nm2.nm_freqs, dstrip(nm.nm_freqs))
