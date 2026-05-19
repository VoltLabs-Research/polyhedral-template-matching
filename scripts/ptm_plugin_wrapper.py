#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PLUGIN_ROOT = Path(os.environ.get('PLUGIN_PROJECT_DIR', SCRIPT_DIR.parent)).resolve()
PLUGINS_ROOT = PLUGIN_ROOT.parent if PLUGIN_ROOT.parent.name == 'plugins' else None
EMBEDDED_LOADER = (PLUGIN_ROOT / 'lib/ld-linux-x86-64.so.2').resolve()
EMBEDDED_LIBRARY_DIR = (PLUGIN_ROOT / 'lib').resolve()


class WrapperError(RuntimeError):
    pass


def resolve_existing_path(label: str, candidates: list[Path | None]) -> Path:
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate.resolve()

    pretty_candidates = '\n'.join(f'  - {candidate}' for candidate in candidates if candidate)
    raise WrapperError(f'No pude resolver {label}. Paths probados:\n{pretty_candidates}')


def resolve_binary_from_env_or_candidates(env_var: str, label: str, candidates: list[Path], which_name: str | None = None) -> Path:
    env_value = os.environ.get(env_var, '').strip()
    env_candidate = Path(env_value) if env_value else None
    which_candidate = Path(shutil.which(which_name)).resolve() if which_name and shutil.which(which_name) else None
    return resolve_existing_path(label, [env_candidate, *candidates, which_candidate])


def resolve_embedded_runtime_command(command: list[str]) -> list[str]:
    if not command:
        return command

    if not EMBEDDED_LOADER.exists() or not EMBEDDED_LIBRARY_DIR.exists():
        return command

    binary_path = Path(command[0]).resolve()
    try:
        binary_path.relative_to(PLUGIN_ROOT)
    except ValueError:
        return command

    if binary_path == EMBEDDED_LOADER:
        return command

    return [
        str(EMBEDDED_LOADER),
        '--library-path', str(EMBEDDED_LIBRARY_DIR),
        str(binary_path),
        *command[1:]
    ]


def run(command: list[str]) -> None:
    command = resolve_embedded_runtime_command(command)
    print(f"[ptm-plugin] {' '.join(shlex.quote(part) for part in command)}", flush=True)
    completed = subprocess.run(command)
    if completed.returncode != 0:
        raise WrapperError(f'El comando falló con exit code {completed.returncode}: {command[0]}')


def ensure_parent_directory(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def require_file(path: Path) -> Path:
    if not path.exists():
        raise WrapperError(f'Falta el archivo requerido: {path}')
    return path


def normalize_structure_token(value: str) -> str:
    return value.strip().replace('-', '_').replace(' ', '_').upper()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='ptm_plugin_wrapper.py',
        description='Wrapper portable para Polyhedral Template Matching.'
    )
    parser.add_argument('input_dump')
    parser.add_argument('output_base')
    parser.add_argument('--crystalStructure', default='FCC')
    parser.add_argument('--rmsd', type=float, default=0.1)
    parser.add_argument('--dissolveSmallClusters', action='store_true')

    args, unknown = parser.parse_known_args()
    if unknown:
        raise WrapperError(f'Argumentos no soportados por PTM wrapper: {unknown}')
    return args


def resolve_runtime_paths() -> dict[str, Path]:
    repo_candidates = []
    if PLUGINS_ROOT is not None:
        repo_candidates = [
            PLUGINS_ROOT / 'PolyhedralTemplateMatching/build-local/build/Release/build/Release/polyhedral-template-matching',
            PLUGINS_ROOT / 'PolyhedralTemplateMatching/build/build/Release/polyhedral-template-matching',
            PLUGINS_ROOT / 'PolyhedralTemplateMatching/build-manual/build/Release/polyhedral-template-matching'
        ]

    return {
        'ptm': resolve_binary_from_env_or_candidates(
            'VOLT_PTM_BINARY',
            'binario PTM',
            [PLUGIN_ROOT / 'bin/polyhedral-template-matching', *repo_candidates],
            'polyhedral-template-matching'
        )
    }


def build_command(args: argparse.Namespace, runtime_paths: dict[str, Path]) -> list[str]:
    command = [
        str(runtime_paths['ptm']),
        args.input_dump,
        args.output_base,
        '--crystalStructure', normalize_structure_token(args.crystalStructure),
        '--rmsd', str(args.rmsd)
    ]
    if args.dissolveSmallClusters:
        command.append('--dissolveSmallClusters')
    return command


def main() -> int:
    args = parse_args()
    runtime_paths = resolve_runtime_paths()
    output_base = Path(args.output_base)
    ensure_parent_directory(output_base)
    run(build_command(args, runtime_paths))
    require_file(Path(f'{args.output_base}_ptm_analysis.msgpack'))
    require_file(Path(f'{args.output_base}_atoms.msgpack'))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except WrapperError as error:
        print(f'[ptm-plugin] error: {error}', file=sys.stderr)
        raise SystemExit(1)


# --- Volt persistent plugin entrypoint (auto-inyectado) ---
_VOLT_RUNTIME_FLAGS_WITH_VALUE = {
    "--selectedTimesteps",
    "--selected-timesteps",
}


def _volt_filter_reserved_args(args):
    filtered = []
    iterator = iter(range(len(args)))
    for i in iterator:
        token = str(args[i])
        if token in _VOLT_RUNTIME_FLAGS_WITH_VALUE:
            next(iterator, None)
            continue
        filtered.append(token)
    return filtered


def _volt_chmod_executables(plugin_root):
    """Ensure binaries inside the extracted plugin are executable.
    The daemon's unzipper extraction does not preserve Unix mode bits."""
    import os as _os
    bin_dir = plugin_root / "bin"
    targets = []
    if bin_dir.is_dir():
        for child in bin_dir.iterdir():
            if child.is_file():
                targets.append(child)
    lib_loader = plugin_root / "lib" / "ld-linux-x86-64.so.2"
    if lib_loader.exists():
        targets.append(lib_loader)
    for target in targets:
        try:
            mode = target.stat().st_mode
            if not (mode & 0o111):
                _os.chmod(target, mode | 0o755)
        except OSError:
            pass


def process(frame, config):  # noqa: D401 - declarative shim
    """Volt persistent plugin entry. Subprocess + prints MUST go to stderr,
    because stdout is reserved for the binary IPC protocol with the stub."""
    del frame
    import os as _volt_os
    import subprocess as _volt_subprocess
    if not isinstance(config, dict):
        raise RuntimeError("config debe ser un dict")
    args = config.get("args")
    if not isinstance(args, list):
        raise RuntimeError("config['args'] debe ser una lista")
    args = _volt_filter_reserved_args(args)

    try:
        _volt_chmod_executables(PLUGIN_ROOT)
    except Exception:
        pass

    original_argv = sys.argv
    sys.argv = [original_argv[0] if original_argv else __file__, *args]
    original_stdout_fd = _volt_os.dup(1)
    _volt_os.dup2(2, 1)
    original_run = _volt_subprocess.run

    def _stderr_run(*run_args, **run_kwargs):
        run_kwargs.setdefault("stdout", _volt_subprocess.DEVNULL)
        run_kwargs.setdefault("stderr", _volt_subprocess.DEVNULL)
        return original_run(*run_args, **run_kwargs)

    _volt_subprocess.run = _stderr_run
    try:
        rc = main()
    except WrapperError as error:
        raise RuntimeError(str(error))
    finally:
        _volt_subprocess.run = original_run
        _volt_os.dup2(original_stdout_fd, 1)
        _volt_os.close(original_stdout_fd)
        sys.argv = original_argv
    if rc not in (0, None):
        raise RuntimeError(f"plugin exit code {rc}")
    return {"ok": True}
