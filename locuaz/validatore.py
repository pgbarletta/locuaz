import os
from pathlib import Path

from cerberus import Validator


class Validatore(Validator):
    def _validate_contains_any_of(self, constraints, field, value):
        """_validate_contains_any_of

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if constraints:
            if len(set(value.keys()).intersection(constraints)) == 0:
                self._error(field, f"Must contain any of: {constraints}")

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
        """_validate_same_length _summary_

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            if len(value) != len(self.document[other]):
                self._error(field, f" should have the same length as {other}")

    def _validate_lower_than_length(self, other, field, value):
        """lower_than_length

        The rule's arguments are validated against this schema:
        {'type': 'string'}
        """
        if other:
            if value > len(self.document[other]):
                self._error(field, f"cannot be higher than the length of {other}.")

    def _validate_warn_when_above(self, threshold, field, value):
        """_validate_warn_when_above

        The rule's arguments are validated against this schema:
        {'type': 'integer'}
        """
        if threshold and value > threshold:
            print(
                f"Warning: {field} set to {value}. Make sure you have enough resources.",
                flush=True,
            )

    def _validate_warn_thread_availability(self, others, field, value):
        """_validate_warn_thread_availability

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if others:
            assert len(others) == 2, f"The schema should have 2 elemnts in "
            f"{field}'s `warn_thread_availability`"
            mpi = self.document[others[0]]
            omp = self.document[others[1]]
            necessary = value * mpi * omp
            try:
                available_procs = int(os.getenv("SLURM_CPUS_ON_NODE"))
            except TypeError as e:
                # TODO: also try to check TORQUE PBS environmental variable.
                available_procs = len(os.sched_getaffinity(0))
            if available_procs < necessary:
                print(
                    f"Warning, {value} gpus, {mpi} MPI processors and {omp} "
                    f"OMP threads requested. {necessary} threads are necessary, "
                    f"but only {available_procs} are available.\n "
                    "Continue only if you know what you're doing.",
                    flush=True,
                )

    def _validate_unique(self, flag, field, value):
        """_validate_unique

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if len(value) != len(set(value)):
                self._error(field, "should be a list of sorted unique values.")

    def _validate_crosscheck(self, others, field, value):
        """_validate_crosscheck

        The rule's arguments are validated against this schema:
        {'type': 'list'}
        """
        if others:
            branches = self.document[others[0]]
            if ("SPM" in value) and (branches > 19):
                self._error(
                    field,
                    f"{others[0]} set to {branches} but this binder generator outputs "
                    "mutations for a single position, hence, 19 is the maximum number "
                    "of possible mutations.",
                )

    def _validate_is_directory(self, flag, field, value):
        """_validate_is_directory

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            if not Path(value).is_dir():
                self._error(field, "should be an existing directory.")

    def _validate_is_iteration_dir(self, flag, field, value):
        """_validate_is_iteration_dir

        The rule's arguments are validated against this schema:
        {'type': 'boolean'}
        """
        if flag:
            nbr, *chains_resnames = Path(value).name.split("-")
            if not nbr.isnumeric():
                self._error(
                    field,
                    f"{value} is not a valid iteration folder. {nbr} is not a "
                    "valid epoch number.",
                )
            for chain_resname in chains_resnames:
                chainID, resname = chain_resname.split("_")
                if not len(chainID) == 1:
                    self._error(
                        field,
                        f"{value} is not a valid iteration folder. {chainID} "
                        "is not a valid chainID.",
                    )
                if not resname.isalpha() or not resname.isupper():
                    self._error(
                        field,
                        f"{value} is not a valid iteration folder. {resname} "
                        "is not a valid resname sequence.",
                    )
