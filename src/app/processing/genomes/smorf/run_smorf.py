# pylint: disable=too-many-arguments
"""
smORFinder execution utilities.

This module provides a thin wrapper around the external `smorf`
command-line tool, responsible for executing smORF prediction
on individual genome FASTA files.

The functions in this module intentionally remain lightweight
and focused on tool execution, leaving output parsing and
post-processing to downstream modules.
"""

from __future__ import annotations

# Standard Library Imports
import subprocess
from pathlib import Path
from typing import List


# Public API
def run_smorf_single_genome(
    *,
    fasta_path: Path,
    output_dir: Path,
    dsn1_cutoff: float = 0.9999,
    dsn2_cutoff: float = 0.9999,
    phmm_cutoff: float = 1e-6,
    smorf_conda_env: str | None = None,
    overwrite: bool = False,
) -> bool:
    """
    Execute smORFinder prediction on a single genome FASTA.

    This function invokes the external `smorf single` command
    and verifies successful execution by checking for the
    expected output artifact produced by the tool.

    Parameters
    ----------
    fasta_path : pathlib.Path
        Path to the input genome FASTA file.
    output_dir : pathlib.Path
        Directory where smORFinder output will be written.
    dsn1_cutoff : float, optional
        Probability cutoff for the DSN1 model.
    dsn2_cutoff : float, optional
        Probability cutoff for the DSN2 model.
    phmm_cutoff : float, optional
        Significance cutoff for profile HMM search.
    overwrite : bool, optional
        Whether to force overwrite of existing output directory.

    Returns
    -------
    bool
        True if smORFinder executed successfully and produced
        the expected output artifact, False otherwise.

    Raises
    ------
    RuntimeError
        If the smORFinder executable cannot be executed.
    """

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = _build_smorf_command(
        fasta_path=fasta_path,
        output_dir=output_dir,
        dsn1_cutoff=dsn1_cutoff,
        dsn2_cutoff=dsn2_cutoff,
        phmm_cutoff=phmm_cutoff,
        smorf_conda_env=smorf_conda_env,
        overwrite=overwrite,
    )

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

    except OSError as exc:
        raise RuntimeError(
            "Failed to execute smORFinder. "
            "Ensure the `smorf` executable is available in PATH."
        ) from exc

    if result.returncode != 0:
        raise RuntimeError(
            "smORFinder execution failed\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDERR:\n{result.stderr}"
        )

    # Verify expected artifact
    if not _has_expected_output(output_dir):
        return False

    return True


# Internal Helpers
def _build_smorf_command(
    *,
    fasta_path: Path,
    output_dir: Path,
    dsn1_cutoff: float,
    dsn2_cutoff: float,
    phmm_cutoff: float,
    smorf_conda_env: str | None,
    overwrite: bool,
) -> List[str]:
    """
    Construct the smORFinder command.

    Parameters
    ----------
    fasta_path : pathlib.Path
        Input genome FASTA file.
    output_dir : pathlib.Path
        Output directory for smORFinder results.
    dsn1_cutoff : float
        DSN1 probability cutoff.
    dsn2_cutoff : float
        DSN2 probability cutoff.
    phmm_cutoff : float
        Profile HMM significance cutoff.
    smorf_conda_env : str | None
        Conda environment name used to execute smORFinder.
        If None, the `smorf` executable will be invoked directly
        from the current runtime environment.
    overwrite : bool
        Whether to force overwrite existing outputs.

    Returns
    -------
    list[str]
        Command argument list suitable for subprocess execution.
    """
    cmd: List[str] = []

    if smorf_conda_env:
        cmd.extend(["conda", "run", "-n", smorf_conda_env])

    cmd.extend(
        [
            "smorf",
            "single",
            str(fasta_path),
            "-o",
            str(output_dir),
            "-idsn1",
            str(dsn1_cutoff),
            "-idsn2",
            str(dsn2_cutoff),
            "-iphmm",
            str(phmm_cutoff),
        ]
    )

    if overwrite:
        cmd.append("--force")

    return cmd


def _has_expected_output(output_dir: Path) -> bool:
    """
    Determine whether smORFinder produced expected output artifacts.

    The presence of `model_predictions.tsv` is used as the
    primary indicator of successful execution.

    Parameters
    ----------
    output_dir : pathlib.Path
        Directory where smORFinder results are written.

    Returns
    -------
    bool
        True if the expected artifact exists, False otherwise.
    """

    expected_file = output_dir / "tmp" / "model_predictions.tsv"
    return expected_file.is_file()
