"""Unit tests for the <ffmpi> input parsing (no MPI runtime needed)."""

from ipi.utils.io.inputs.io_xml import xml_parse_string
from ipi.inputs.forcefields import InputFFMPI
from ipi.inputs.simulation import InputSimulation


def _ffmpi_node(xml):
    return xml_parse_string(xml).fields[0][1]


def test_inputffmpi_parse_and_fetch():
    node = _ffmpi_node(
        "<ffmpi name='ffX' mode='mpi'><batch_size>4</batch_size></ffmpi>"
    )
    inp = InputFFMPI()
    inp.parse(node)
    ff = inp.fetch()
    assert ff.name == "ffX"
    assert ff.mode == "mpi"
    assert ff.interface.batch_size == 4
    # the driver ranks are polled every `latency` seconds
    assert ff.poll_interval == ff.latency


def test_inputffmpi_store_roundtrip():
    node = _ffmpi_node(
        "<ffmpi name='ffY' mode='mpi'><batch_size>2</batch_size></ffmpi>"
    )
    built = InputFFMPI()
    built.parse(node)
    ff = built.fetch()

    stored = InputFFMPI()
    stored.store(ff)
    assert stored.name.fetch() == "ffY"
    assert stored.mode.fetch() == "mpi"
    assert stored.batch_size.fetch() == 2


def test_inputffmpi_default_batch_size():
    node = _ffmpi_node("<ffmpi name='ffZ' mode='mpi'></ffmpi>")
    inp = InputFFMPI()
    inp.parse(node)
    assert inp.fetch().interface.batch_size == 1


def test_simulation_dispatches_ffmpi():
    """An <ffmpi> block inside <simulation> is routed to InputFFMPI (covers the
    ffmpi entry in the forcefield dispatch table in inputs/simulation.py)."""
    xml = """<simulation>
      <ffmpi name='ffX' mode='mpi'><batch_size>2</batch_size></ffmpi>
      <system>
        <initialize nbeads='1'>
          <cell mode='abc'>[10,10,10]</cell>
          <file mode='xyz'> dummy </file>
        </initialize>
        <forces><force forcefield='ffX'/></forces>
        <ensemble><temperature units='kelvin'>300</temperature></ensemble>
        <motion mode='dynamics'>
          <dynamics mode='nvt'>
            <timestep units='femtosecond'>0.5</timestep>
            <thermostat mode='langevin'><tau units='femtosecond'>25</tau></thermostat>
          </dynamics>
        </motion>
      </system>
    </simulation>"""
    isim = InputSimulation()
    isim.parse(xml_parse_string(xml).fields[0][1])
    ffmpi_blocks = [obj for tag, obj in isim.extra if tag == "ffmpi"]
    assert len(ffmpi_blocks) == 1
    assert isinstance(ffmpi_blocks[0], InputFFMPI)
    assert ffmpi_blocks[0].batch_size.fetch() == 2
