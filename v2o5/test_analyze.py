"""Unit tests for analyze.py parsers (run: python3 test_analyze.py)."""
from pathlib import Path
from analyze import total_energy_ry, eads_ev

SAMPLE_OK = """
     iteration #  5
!    total energy              =     -100.50000000 Ry
     convergence has been achieved
     bfgs converged in 3 scf cycles
     JOB DONE.
"""
SAMPLE_UNCONV = "!    total energy = -1.0 Ry\n"   # no JOB DONE / no bfgs converged


def test_total_energy_reads_last_bang():
    p = Path("/tmp/_v2o5_t.out"); p.write_text(SAMPLE_OK)
    assert abs(total_energy_ry(p) - (-100.5)) < 1e-9


def test_total_energy_none_if_unconverged():
    p = Path("/tmp/_v2o5_t2.out"); p.write_text(SAMPLE_UNCONV)
    assert total_energy_ry(p) is None


def test_total_energy_none_if_missing():
    assert total_energy_ry(Path("/tmp/_v2o5_does_not_exist.out")) is None


def test_eads_ev():
    # (E_complex - E_slab - E_gas) Ry -> eV
    assert abs(eads_ev(-110.0, -100.0, -9.5) - (-0.5 * 13.605693)) < 1e-6


if __name__ == "__main__":
    test_total_energy_reads_last_bang()
    test_total_energy_none_if_unconverged()
    test_total_energy_none_if_missing()
    test_eads_ev()
    print("analyze tests PASS")
