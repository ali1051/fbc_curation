"""Testing result."""
from pathlib import Path
from typing import Dict, Tuple

import libsbml
import pandas as pd
import pytest

from fbc_curation import EXAMPLE_DIR
from fbc_curation.compare import FrogComparison
from fbc_curation.curator.cobrapy_curator import CuratorCobrapy
from fbc_curation.frog import Creator, CuratorConstants, FrogReport


model_path: Path = EXAMPLE_DIR / "models" / "e_coli_core.xml"
frog_id: str = "1234"
curators = [
    Creator(
        familyName="König",
        givenName="Matthias",
        email=None,
        organization=None,
        site=None,
        orcid=None
    )
]
curator = CuratorCobrapy(model_path=model_path, frog_id=frog_id, curators=curators)
report: FrogReport = curator.run()
doc: libsbml.SBMLDocument = libsbml.readSBMLFromFile(str(model_path))
model: libsbml.Model = doc.getModel()


@pytest.fixture(scope="module")
def setup_report() -> Tuple[FrogReport, libsbml.Model]:
    model_path_fixture: Path = EXAMPLE_DIR / "models" / "e_coli_core.xml"
    frog_id_fixture: str = "1234"
    curators_fixture = [
        Creator(
            familyName="König",
            givenName="Matthias",
            email=None,
            organization=None,
            site=None,
            orcid=None
        )
    ]
    curator_fixture = CuratorCobrapy(model_path=model_path_fixture, frog_id=frog_id_fixture, curators=curators_fixture)
    report_fixture: FrogReport = curator_fixture.run()
    doc_fixture: libsbml.SBMLDocument = libsbml.readSBMLFromFile(str(model_path_fixture))
    model_fixture: libsbml.Model = doc_fixture.getModel()
    return report_fixture, model_fixture


def test_objective_df(setup_report: Tuple[FrogReport, libsbml.Model]) -> None:
    """Check objective."""
    report, _ = setup_report
    dfs: Dict[str, pd.DataFrame] = report.to_dfs()
    df = dfs[CuratorConstants.OBJECTIVE_KEY]

    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    obj_value = df["value"].values[0]
    assert obj_value > 0

    status_codes = df.status.unique()
    assert len(status_codes) <= 2
    assert "optimal" in status_codes


def test_fva_df(setup_report: Tuple[FrogReport, libsbml.Model]) -> None:
    """Check FVA DataFrame."""
    report, model = setup_report
    dfs: Dict[str, pd.DataFrame] = report.to_dfs()
    df = dfs[CuratorConstants.FVA_KEY]

    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    assert len(df) == model.getNumReactions()

    status_codes = df.status.unique()
    assert len(status_codes) <= 2
    assert "optimal" in status_codes


def test_gene_deletion_df(setup_report: Tuple[FrogReport, libsbml.Model]) -> None:
    """Check gene deletion."""
    report, model = setup_report
    dfs: Dict[str, pd.DataFrame] = report.to_dfs()
    df = dfs[CuratorConstants.GENEDELETIONS_KEY]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    fbc_model: libsbml.FbcModelPlugin = model.getPlugin("fbc")
    assert len(df) == fbc_model.getNumGeneProducts()

    status_codes = df.status.unique()
    assert len(status_codes) <= 2
    assert "optimal" in status_codes


def test_reaction_deletion_df(setup_report: Tuple[FrogReport, libsbml.Model], tmp_path: Path) -> None:
    """Check reaction deletion."""
    report, model = setup_report
    dfs: Dict[str, pd.DataFrame] = report.to_dfs()
    df: pd.DataFrame = dfs[CuratorConstants.REACTIONDELETIONS_KEY]
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    assert len(df) == model.getNumReactions()

    status_codes = df.status.unique()
    assert len(status_codes) <= 2
    assert "optimal" in status_codes


def test_report_write_read_tsv_equal(setup_report: Tuple[FrogReport, libsbml.Model], tmp_path: Path) -> None:
    """Test equality of report after writing/reading TSV."""
    report, _ = setup_report
    report.to_tsv(tmp_path)
    report2 = FrogReport.from_tsv(tmp_path)

    assert FrogComparison.compare_reports({"report": report, "report2": report2})


def test_report_write_read_json_equal(setup_report: Tuple[FrogReport, libsbml.Model], tmp_path: Path) -> None:
    """Test equality of report after writing/reading JSON."""
    report, _ = setup_report
    report.to_json(path=tmp_path / "frog.json")
    report2 = FrogReport.from_json(path=tmp_path / "frog.json")

    assert FrogComparison.compare_reports({"report": report, "report2": report2})
