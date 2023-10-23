from typing import Iterator, TypeVar, Generic, Set, Iterable
from random import choice

T = TypeVar("T")


class InfinitePool(Generic[T]):
    """Self-replacing set. It takes a set to initialize itself, and then it will
    only be empty when the last element has been picked. Future accesses will
    reset it to the original set."""
    pool: Set[T]
    reservoir: Set[T]

    def __init__(self, input_set: Iterable[T]) -> None:
        self.reservoir = set(input_set)
        self.pool = set(self.reservoir)

    def __reset__(self) -> None:
        if len(self.pool) == 0:
            self.pool = set(self.reservoir)

    def __contains__(self, item) -> bool:
        return item in self.pool

    def __iter__(self) -> Iterator:
        return iter(self.pool)

    def __len__(self) -> T:
        return len(self.pool)

    def __sub__(self, other: "InfinitePool") -> "InfinitePool":
        self.__reset__()
        return InfinitePool(self.pool - other.pool)

    def __isub__(self, other: "InfinitePool") -> "InfinitePool":
        self.__reset__()
        self.pool -= other.pool
        return self

    def __or__(self, other: "InfinitePool") -> "InfinitePool":
        self.__reset__()
        return InfinitePool(self.pool | other.pool)

    def __and__(self, other: "InfinitePool") -> "InfinitePool":
        self.__reset__()
        return InfinitePool(self.pool & other.pool)

    def __repr__(self) -> str:
        return self.pool.__repr__()

    def __str__(self) -> str:
        return self.pool.__str__()

    def force_reset(self) -> None:
        self.pool = set(self.reservoir)

    def discard(self, item) -> None:
        self.__reset__()
        return self.pool.discard(item)

    def pick(self) -> T:
        self.__reset__()
        return choice(tuple(self.pool))

    def pop(self) -> T:
        self.__reset__()
        item = self.pick()
        self.discard(item)
        return item

    def difference_update(self, other: Set) -> None:
        self.__reset__()
        self.pool.difference_update(other)


class BinPool(InfinitePool):
    """InfinitePool specifically tailored for ints. Used to represent the
    indices of the amino acid bins when selecting them for mutation creation."""

    def __init__(self, input_set: Iterable[int]) -> None:
        super().__init__(input_set)

    @classmethod
    def from_size(cls, size: int) -> "BinPool":
        return cls(set(range(0, size)))
