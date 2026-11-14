############################################
# Copyright (c) 2026 Microsoft Corporation
#
# Dependency-free tests for the Lean checking helper.
############################################
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT = _ROOT / "scripts" / "check_lean.sh"
_BASH = shutil.which("bash")
_TOOLCHAIN = (_ROOT / "lean" / "lean-toolchain").read_text().strip()
_FAKE_ELAN = """#!/bin/bash
printf '%s\\0' "$@" >> "$ELAN_TEST_LOG"
printf '\\0' >> "$ELAN_TEST_LOG"
printf '%s\\n' "$PWD" >> "$ELAN_TEST_CWD"
printf '%s\\n' "${LEAN_PATH:-}" >> "$ELAN_TEST_CWD.paths"
if [[ "${4:-}" == "build" ]]; then
    exit "${ELAN_TEST_BUILD_STATUS:-0}"
fi
exit "${ELAN_TEST_LEAN_STATUS:-0}"
"""


@unittest.skipUnless(os.name == "posix" and _BASH, "requires Bash on POSIX")
class TestCheckLean(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name).resolve()
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.elan = self.bin / "elan"
        self.elan.write_text(_FAKE_ELAN)
        self.elan.chmod(0o755)
        self.log = self.root / "calls"
        self.cwd_log = self.root / "directories"
        self.env = dict(os.environ)
        self.env.update({
            "PATH": str(self.bin) + os.pathsep + os.defpath,
            "HOME": str(self.root / "home"),
            "ELAN_TEST_LOG": str(self.log),
            "ELAN_TEST_CWD": str(self.cwd_log),
        })
        self.proof = self.root / "example proof.txt"
        self.proof.write_text("example : True := True.intro\n")

    def run_helper(self, *arguments):
        return subprocess.run(
            [_BASH, str(_SCRIPT), *arguments],
            cwd=self.root, env=self.env, text=True, capture_output=True, check=False,
        )

    def calls(self):
        if not self.log.exists():
            return []
        return [record.decode().split("\0")
                for record in self.log.read_bytes().split(b"\0\0") if record]

    def test_default_builds_the_pinned_project(self):
        result = self.run_helper()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls(), [["run", _TOOLCHAIN, "lake", "build"]])
        self.assertEqual(self.cwd_log.read_text().splitlines(), [str(_ROOT / "lean")])
        self.assertIn("Lean checks passed.", result.stdout)

    def test_relative_paths_spaces_and_strict_flags(self):
        result = self.run_helper(self.proof.name)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls(), [
            ["run", _TOOLCHAIN, "lake", "build"],
            ["run", _TOOLCHAIN, "lake", "env", "lean",
             "--trust=0", "-DwarningAsError=true", "--memory=8192", str(self.proof)],
        ])
        self.assertEqual(self.cwd_log.read_text().splitlines(), [str(_ROOT / "lean")] * 2)

    def test_library_build_has_the_same_flags(self):
        lakefile = (_ROOT / "lean" / "lakefile.toml").read_text()
        self.assertIn('moreLeanArgs = ["--trust=0", "-DwarningAsError=true", "--memory=8192"]',
                      lakefile.splitlines())

    def test_absolute_paths_and_multiple_files(self):
        second = self.root / "second.lean"
        second.write_text("example : True := True.intro\n")
        result = self.run_helper("--", str(self.proof), str(second))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([call[-1] for call in self.calls()[1:]], [str(self.proof), str(second)])

    def test_missing_file_is_rejected_before_building(self):
        result = self.run_helper("missing.lean")
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot read file", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_directory_is_not_a_source_file(self):
        result = self.run_helper(str(self.root))
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot read file", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_build_failure_stops_checking(self):
        self.env["ELAN_TEST_BUILD_STATUS"] = "7"
        result = self.run_helper(str(self.proof))
        self.assertEqual(result.returncode, 7)
        self.assertEqual(len(self.calls()), 1)
        self.assertNotIn("Lean checks passed.", result.stdout)

    def test_lean_failure_is_preserved_and_stops_later_files(self):
        self.env["ELAN_TEST_LEAN_STATUS"] = "9"
        result = self.run_helper(str(self.proof), str(self.proof))
        self.assertEqual(result.returncode, 9)
        self.assertEqual(len(self.calls()), 2)
        self.assertNotIn("Lean checks passed.", result.stdout)

    def chain(self, count=3):
        files = []
        for index in range(count):
            files.append(self.root / ("Module_%d.lean" % index))
            files[-1].write_text("example : True := True.intro\n")
        return files

    def test_chain_compiles_every_file_but_the_last_for_the_later_ones(self):
        files = self.chain()
        self.env["LEAN_PATH"] = "/somewhere/else"
        result = self.run_helper("--chain", files[0].name, str(files[1]), files[2].name)
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual(calls[0], ["run", _TOOLCHAIN, "lake", "build"])
        strict = ["run", _TOOLCHAIN, "lake", "env", "lean", "--trust=0", "-DwarningAsError=true", "--memory=8192"]
        objects = []
        for call, path in zip(calls[1:3], files[:2]):
            # The same flags as a single file; the module name is the file's name under its directory.
            self.assertEqual(call[:len(strict)], strict)
            self.assertEqual(call[len(strict)], "--root=" + str(self.root))
            self.assertEqual(call[len(strict) + 1], "-o")
            objects.append(Path(call[len(strict) + 2]))
            self.assertEqual(call[len(strict) + 3:], [str(path)])
        self.assertEqual([path.name for path in objects], ["Module_0.olean", "Module_1.olean"])
        # The last file is checked as a single file is: nothing imports it.
        self.assertEqual(calls[3], strict + [str(files[2])])
        self.assertEqual(len(calls), 4)
        # One fresh directory holds the chain's .olean files, is on the search path of its Lean
        # runs (before what the caller had there), and is gone when the run ends.
        directory = objects[0].parent
        self.assertEqual(objects[1].parent, directory)
        self.assertNotIn(directory, (self.root, _ROOT / "lean"))
        self.assertFalse(directory.exists())
        paths = Path(str(self.cwd_log) + ".paths").read_text().splitlines()
        self.assertEqual(paths[1:], ["%s:/somewhere/else" % directory] * 3)
        self.assertIn("Lean checks passed.", result.stdout)

    def test_chain_of_one_file_is_a_single_file(self):
        result = self.run_helper("--chain", "--", str(self.proof))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls()[1][-4:], ["--trust=0", "-DwarningAsError=true", "--memory=8192", str(self.proof)])
        self.assertEqual(self.run_helper("--chain").returncode, 0)

    def test_chain_rejects_an_imported_file_whose_name_is_not_a_module_name(self):
        files = self.chain(2)
        for name in ("example proof.txt", "1st.lean", "a.b.lean", "a-b.lean"):
            with self.subTest(name=name):
                (self.root / name).write_text("example : True := True.intro\n")
                result = self.run_helper("--chain", name, str(files[1]))
                self.assertEqual(result.returncode, 2)
                self.assertIn("must be named IDENTIFIER.lean", result.stderr)
                self.assertEqual(self.calls(), [])
        # The last file of a chain is imported by nothing, so any source file will do.
        result = self.run_helper("--chain", str(files[0]), str(self.proof))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.calls()), 3)

    def test_chain_stops_at_the_first_module_lean_rejects(self):
        files = self.chain()
        self.env["ELAN_TEST_LEAN_STATUS"] = "9"
        result = self.run_helper("--chain", *map(str, files))
        self.assertEqual(result.returncode, 9)
        calls = self.calls()
        self.assertEqual(len(calls), 2)
        self.assertFalse(Path(calls[1][-2]).parent.exists())  # the .olean directory is removed on failure too
        self.assertNotIn("Lean checks passed.", result.stdout)

    def test_missing_elan_has_an_actionable_error(self):
        self.env["PATH"] = os.defpath
        result = self.run_helper()
        self.assertEqual(result.returncode, 127)
        self.assertIn("elan was not found", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_standard_elan_home_is_a_fallback(self):
        fallback = Path(self.env["HOME"]) / ".elan" / "bin" / "elan"
        fallback.parent.mkdir(parents=True)
        shutil.copyfile(self.elan, fallback)
        fallback.chmod(0o755)
        self.env["PATH"] = os.defpath
        result = self.run_helper(str(self.proof))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.calls()), 2)

    def test_help_does_not_need_elan(self):
        self.env["PATH"] = os.defpath
        result = self.run_helper("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("Usage:", result.stdout)
        self.assertEqual(self.calls(), [])


if __name__ == "__main__":
    unittest.main()
