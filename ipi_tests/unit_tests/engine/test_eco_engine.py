"""Engine-level tests for economised (eco) path integrals: input parsing,
binding, the omegak depend chain, spring/estimator consistency, propagation
and the guards against incompatible options."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import pytest
import numpy as np

from ipi.inputs.normalmodes import InputNormalModes
from ipi.utils import nmtransform
from ipi.utils.array_backend import to_numpy
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
    # spreads the beads so that the internal modes are excited
    q = np.array([0.0, 0.0, 0.0, 2.0, 0.0, 0.0] * NBEADS).reshape(
        NBEADS, 3 * NATOMS
    ) + 0.1 * rng.normal(size=(NBEADS, 3 * NATOMS))

    return f"""
<simulation verbosity='quiet'>
  <output prefix='sim'>
    <properties stride='100'> [ step, conserved ] </properties>
  </output>
  <total_steps> 100 </total_steps>
  <prng> <seed> 31415 </seed> </prng>
  <ffdirect name='harm'> <pes> harmonic </pes> <parameters> {{k1: 0.01}} </parameters> </ffdirect>
  <system>
    <beads natoms='{NATOMS}' nbeads='{NBEADS}'>
      <q shape='({NBEADS},{3 * NATOMS})'> {q.ravel().tolist()} </q>
      <m shape='({NATOMS})'> {[1837.36] * NATOMS} </m>
      <names shape='({NATOMS})'> [{", ".join(["H"] * NATOMS)}] </names>
    </beads>
    <cell shape='(3,3)'> [20, 0, 0, 0, 20, 0, 0, 0, 20] </cell>
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
def eco_sim(simulation):
    with simulation(eco_xml()) as sim:
        yield sim


def test_omegak_and_masses(eco_sim):
    """Checks the eco frequencies against the fit, and the dynamical masses."""

    nm = eco_sim.syslist[0].nm
    omegak = to_numpy(nm.omegak)
    xmax = float(nm.nm_freqs[0] * nm.nbeads / nm.omegan)
    np.testing.assert_allclose(
        omegak, nm.omegan * nmtransform.eco_eva(nm.nbeads, xmax), rtol=1e-10
    )
    assert omegak[0] == 0.0
    # eco changes the springs but not the dynamical masses
    np.testing.assert_allclose(to_numpy(nm.nm_factor), np.ones(NBEADS))


def test_spring_and_kinetic_td(eco_sim):
    """Checks the spring energy against an explicit normal-mode sum, and the
    consistency of the primitive kinetic energy estimator with it."""

    sys = eco_sim.syslist[0]
    nm = sys.nm
    qnm = to_numpy(nm.qnm)
    m3 = to_numpy(sys.beads.m3)
    wk2 = to_numpy(nm.omegak) ** 2
    vspring = 0.5 * np.sum(wk2[:, np.newaxis] * m3 * qnm**2)
    np.testing.assert_allclose(float(nm.vspring), vspring, rtol=1e-12)

    ktd = float(sys.properties["kinetic_td"][0])
    kt = Constants.kb * sys.ensemble.temp
    np.testing.assert_allclose(
        ktd, 1.5 * NATOMS * NBEADS * kt - vspring / NBEADS, rtol=1e-10
    )


def test_temperature_refit(eco_sim):
    """Checks that a temperature change refits the frequencies exactly
    through the depend chain (exercising the warm-start path)."""

    sys = eco_sim.syslist[0]
    nm = sys.nm
    gamma_old = np.array(to_numpy(nm.omegak)) / nm.omegan
    sys.ensemble.temp = sys.ensemble.temp * 1.13
    omegak = to_numpy(nm.omegak)
    xmax = float(nm.nm_freqs[0] * nm.nbeads / nm.omegan)
    np.testing.assert_allclose(
        omegak, nm.omegan * nmtransform.eco_eva(nm.nbeads, xmax), rtol=1e-8
    )
    # the dimensionless spectrum must change shape, not just rescale
    gamma_new = omegak / nm.omegan
    assert not np.allclose(gamma_new[1:], gamma_old[1:], rtol=1e-6)


@pytest.mark.parametrize("propagator", ["exact", "cayley", "bab"])
def test_eco_dynamics(simulation, propagator):
    """Runs a few steps of eco PIMD with each free-ring-polymer propagator."""

    with simulation(eco_xml(propagator=propagator)) as sim:
        e0 = sim.properties("conserved")
        sim.run(4, write_outputs=False)
        assert np.isfinite(sim.properties("conserved"))
        assert abs(sim.properties("conserved") - e0) < 1e-2 * abs(e0)


def test_guard_bosons(simulation):
    with pytest.raises(ValueError, match="bosons or open paths"):
        with simulation(eco_xml(nm_extra="<bosons id='index'> [0] </bosons>")):
            pass


def test_guard_open_paths(simulation):
    with pytest.raises(ValueError, match="bosons or open paths"):
        with simulation(eco_xml(nm_extra="<open_paths> [0] </open_paths>")):
            pass


def test_guard_nm_freqs(simulation):
    with pytest.raises(ValueError, match="requires one frequency"):
        with simulation(
            eco_xml(
                frequencies="<frequencies style='eco' units='inversecm'> [2500, 3000] </frequencies>"
            )
        ):
            pass


def test_guard_suzuki_chin(simulation):
    with pytest.raises(ValueError, match="Suzuki-Chin"):
        with simulation(eco_xml(dynamics="sc")):
            pass


def test_guard_nm_gle(simulation):
    a_matrix = str([1e-3] * NBEADS)
    thermostat = f"<thermostat mode='nm_gle'> <A shape='({NBEADS},1,1)'> {a_matrix} </A> </thermostat>"
    with pytest.raises(ValueError, match="nm_gle"):
        with simulation(eco_xml(thermostat=thermostat)):
            pass


def test_guard_instanton(simulation):
    motion = """<motion mode='instanton'>
      <instanton mode='rate'> <opt> nichols </opt> </instanton>
    </motion>"""
    with pytest.raises(ValueError, match="Instanton"):
        with simulation(eco_xml(motion=motion)):
            pass


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


def test_trotter_estimators_unaffected(simulation):
    """The same estimators keep working with Trotter springs."""

    xml = eco_xml(frequencies="<frequencies style='rpmd'> [] </frequencies>")
    with simulation(xml) as sim:
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
    np.testing.assert_allclose(nm2.nm_freqs, nm.nm_freqs)
