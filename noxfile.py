"""Nox sessions."""
import os
import shutil
from pathlib import Path

import nox  # type: ignore
from nox_poetry import Session  # type: ignore
from nox_poetry import session  # type: ignore


package = "haptools"
python_versions = ["3.7", "3.8", "3.9", "3.10", "3.11", "3.12", "3.13"]
locked_python_version = "3.9"  # keep in sync with dev-env.yml
nox.needs_version = ">= 2022.11.21"
nox.options.sessions = (
    "docs",
    "lint",
    "tests",
)


@session(python=locked_python_version)
def docs(session: Session) -> None:
    """Build the documentation."""
    args = session.posargs or ["docs", "docs/_build"]
    if not session.posargs and "FORCE_COLOR" in os.environ:
        args.insert(0, "--color")

    build_dir = Path("docs", "_build")
    if build_dir.exists():
        shutil.rmtree(build_dir)

    session.run("sphinx-build", *args)


@session(python=locked_python_version)
def lint(session: Session) -> None:
    """Lint our code."""
    session.install("black")
    session.run("black", "--verbose", "--check", ".")


def install_handle_python_numpy(session):
    """
    handle incompatibilities with python and numpy versions
    see https://github.com/cjolowicz/nox-poetry/issues/1116
    """
    if session._session.python in ["3.11", "3.12", "3.13"]:
        session._session.install(".")
    else:
        session.install(".")


# detect whether conda/mamba is installed
if os.getenv("CONDA_EXE"):
    conda_cmd = "conda"
    if (Path(os.getenv("CONDA_EXE")).parent / "mamba").exists():
        conda_cmd = "mamba"
    conda_args = ["-c", "conda-forge"]

    @session(venv_backend=conda_cmd, venv_params=conda_args, python=python_versions)
    def tests(session: Session) -> None:
        """Run the test suite."""
        session.conda_install(
            "coverage[toml]",
            "pytest",
            "numpy>=1.20.0",
            channel="conda-forge",
        )
        install_handle_python_numpy(session)
        try:
            session.run(
                "coverage", "run", "--parallel", "-m", "pytest", *session.posargs
            )
        finally:
            if session.interactive:
                session.notify("coverage", posargs=[])

else:

    @session(python=python_versions)
    def tests(session: Session) -> None:
        """Run the test suite."""
        session.install("coverage[toml]", "pytest")
        install_handle_python_numpy(session)
        try:
            session.run(
                "coverage", "run", "--parallel", "-m", "pytest", *session.posargs
            )
        finally:
            if session.interactive:
                session.notify("coverage", posargs=[])


@session(python=locked_python_version)
def coverage(session: Session) -> None:
    """Produce the coverage report."""
    session.install("coverage[toml]")
    args = session.posargs or ["report"]

    if not session.posargs and any(Path().glob(".coverage.*")):
        session.run("coverage", "combine")

    session.run("coverage", *args)
