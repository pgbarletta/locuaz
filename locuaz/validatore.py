# type: ignore
# seems that cerberus doesn't provide proper stubs for mypy
import glob
import os
from pathlib import Path
from warnings import warn

from cerberus import Validator


class Validatore(Validator):
    def _validate_contains_any_of(self, constraints, field, value):
        """_validate_contains_any_of

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if constraints:
            set_fields = set(value.keys())
            if len(set_fields.intersection(constraints)) == 0:
                self._error(field, f"Must contain any of: {constraints}")

            if ("current_branches" in set_fields) and (
                "previous_branches" not in set_fields
            ):
                warn(
                    "Warning: `current_branches` is set, but `previous_branches` isn't. "
                    "Won't be able to prune the current branches. Make sure there are enough branches."
                )

    def _validate_step_bigger_than(self, other, field, value):
        """_validate_step_bigger_than

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            n = self.document[other]
            if n > value[1]:
                self._error(
                    field,
                    f"step ({value[1]}) is lower than `{other}`({n}). "
                    "This would dedicate threads to more than 1 run.",
                )

    def _validate_sorted(self, flag, field, value):
        """_validate_sorted

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if sorted(value) != value:
                self._error(field, "should be an incrementally sorted list.")

    def _validate_same_length(self, other, field, value):
        """_validate_same_length

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            try:
                if len(value) != len(self.document[other]):
                    self._error(field, f" should have the same length as {other}")
            except KeyError:
                # If we've gotten to this point, other must be optional, so this is ok.
                pass

    def _validate_same_as_length_of(self, other, field, value):
        """_validate_lower_than_length_of

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            try:
                if value != len(self.document[other]):
                    self._error(field, f" and the length of {other} must be the same.")
            except KeyError:
                pass

    def _validate_lower_than_length_of(self, other, field, value):
        """_validate_lower_than_length_of

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            try:
                if value > len(self.document[other]):
                    self._error(field, f"cannot be higher than the length of {other}.")
            except KeyError:
                pass

    def _validate_higher_than_length_of(self, other, field, value):
        """_validate_higher_than_length_of

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other and other in self.document:
            try:
                if value < len(self.document[other]):
                    self._error(field, f"cannot be lower than the length of {other}.")
            except KeyError:
                pass

    def _validate_higher_than(self, other, field, value):
        """_validate_higher_than

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other and other in self.document:
            try:
                if value <= self.document[other]:
                    self._error(field, f"cannot be lower than {other}.")
            except KeyError:
                pass

    def _validate_warn_when_above(self, threshold, field, value):
        """_validate_warn_when_above

        The rule's arguments are validated against this schema:
        {'type': 'integer'}
        """
        if threshold and value > threshold:
            warn(
                f"Warning: {field} set to {value}. Make sure you have enough resources."
            )

    def _validate_warn_thread_availability(self, others, field, value):
        """_validate_warn_thread_availability

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if others:
            assert len(others) == 2, f"The schema should have 2 elements in "
            f"{field}'s `warn_thread_availability`"
            mpi = self.document[others[0]]
            omp = self.document[others[1]]
            necessary = value * mpi * omp
            try:
                available_procs = int(os.getenv("SLURM_CPUS_ON_NODE"))
            except TypeError:
                # TODO: also try to check TORQUE PBS environmental variable.
                available_procs = len(os.sched_getaffinity(0))
            if available_procs < necessary:
                warn(
                    f"Warning, {value} gpus, {mpi} MPI processors and {omp} "
                    f"OMP threads requested. {necessary} threads are necessary, "
                    f"but only {available_procs} are available.\n "
                    "Continue only if you know what you're doing.",
                )

    def _validate_unique(self, flag, field, value):
        """_validate_unique

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if len(value) != len(set(value)):
                self._error(field, "should be a list of sorted unique values.")

    def _validate_is_file(self, flag, field, value):
        """_validate_is_file

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if not Path(value).is_file():
                self._error(field, "should be an existing file.")

    def _validate_is_directory(self, flag, field, value):
        """_validate_is_directory

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if not Path(value).is_dir():
                self._error(field, "should be an existing directory.")

    def _validate_is_branch_dir(self, flag, field, value):
        """_validate_is_branch_dir

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            try:
                nbr, *chains_resnames = Path(value).name.split("-")
            except ValueError:
                self._error(field, f"{value} is not a valid branch folder.")
                return

            if not nbr.isnumeric():
                self._error(
                    field,
                    f"{value} is not a valid branch folder. {nbr} is not a "
                    "valid epoch number.",
                )
            for chain_resname in chains_resnames:
                try:
                    chainID, resname = chain_resname.split("_")
                except ValueError:
                    self._error(field, f"{value} is not a valid branch folder.")
                    return

                if not len(chainID) == 1:
                    self._error(
                        field,
                        f"{value} is not a valid branch folder. {chainID} "
                        "is not a valid chainID.",
                    )
                if not resname.isalpha() or not resname.isupper():
                    self._error(
                        field,
                        f"{value} is not a valid branch folder. {resname} "
                        "is not a valid resname sequence.",
                    )

    def _validate_contains_branch_dirs(self, flag, field, value):
        """_validate_contains_branch_dirs

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        for filename in glob.glob(str(Path(value, "*"))):
            if Path(filename).is_dir():
                self._validate_is_branch_dir(flag, field, filename)

    def _validate_warn_overrides(self, other, field, value):
        """_validate_warn_overrides

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other and self.document.get(other, False):
            warn(
                f"Warning: both `{field}` and {other} are set, the former will override the latter. "
            )

    def _validate_scoring_end(self, flag, field, value):
        """_validate_scoring_end

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            start = self.document.get("start", False)
            if value != -1 and value < start:
                self._error(field, f"{value} is less than the starting frame: {start} ")

    def _validate_unique_values(self, flag, field, value):
        """
        _validate_unique_values
        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if len(set(value)) != len(value):
                self._error(field, f"should have unique values.")

    def _validate_crosscheck_radius(self, flag, field, value):
        """_validate_crosscheck_radius

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            r = self.document.get("reconstruct_radius")
            if value == "dlpr" and not r:
                self._error(field, f"{value} needs a `reconstruct_radius`")
            elif value != "dlpr" and r:
                warn(f"Warning: `{value}` does not make use of `reconstruct_radius`.")

    def _validate_forbidden_if_true_mandatory_if_false(self, others, field, value):
        """_validate_forbidden_if_true_mandatory_if_false

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if others:
            forbidden_opts = set(others)
            specified_opts = set(self.document.keys())
            if value:
                if len(forbidden_opts.intersection(specified_opts)) != 0:
                    self._error(
                        field,
                        f" is set. The following options cannot be also set: {others}.",
                    )
            elif not forbidden_opts.issubset(specified_opts):
                self._error(
                    field, f" is not set. The following options must be set: {others}."
                )

    def _validate_enforce_true(self, relationship, field, value):
        """_validate_enforce_true

        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        for key, other_opt in relationship.items():
            if key == value and not self.document[other_opt]:
                self._error(field, f" is set to {key}. {other_opt} must be true.")

    def _validate_warn_dependency_mutator(self, relationship, field, value):
        """_validate_warn_dependency_mutator

        The rule's arguments are validated against this schema:
        {'type': 'dict'}
        """
        for other_opt, other_opt_val in relationship.items():
            allowed_nonstandard_residues_is_set = len(value) >= 0
            not_dlpr = self.document[other_opt] != other_opt_val
            if allowed_nonstandard_residues_is_set and not_dlpr:
                warn(
                    f"`allowed_nonstandard_residues` is set, you should "
                    "probably be using `dlpr` as a mutator. evoef2 can't see "
                    "ligands and dlp may create clashes."
                )
