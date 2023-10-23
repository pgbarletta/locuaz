from typing import Final

from attrs import define, field

@define(frozen=True)
class Site:
    """
    Site represents a site to be mutated. It's immutable and besides carrying
    the chainID and the resSeq for the site, it also carries the indices that
    correspond to the chainID and the resSeq in the mutating chainIDs and
    mutating resSeq

    Parameters
    ----------
    idx_chain : Final[int]
    chainID : Final[str]
    idx_residue : Final[int]
    resSeq : Final[int]
    """
    idx_chain: Final[int] = field(converter=int, kw_only=True)
    chainID: Final[str] = field(converter=str, kw_only=True)
    idx_residue: Final[int] = field(converter=int, kw_only=True)
    resSeq: Final[int] = field(converter=int, kw_only=True)


@define(frozen=True)
class Mutation:
    """
    Mutation represents a Mutation. It's immutable and has the same info as Site
    plus the amino acid that's currently located at that position, and the new
    one to replace it. This means that, unlike Site, Mutation is tied to a
    particular Branch.

    Parameters
    ----------
    idx_chain : Final[int]
    chainID : Final[str]
    idx_residue : Final[int]
    resSeq : Final[int]
    """
    chainID: str = field(converter=str, kw_only=True)
    resSeq: int = field(converter=int, kw_only=True)
    old_aa: str = field(converter=str, kw_only=True)
    new_aa: str = field(converter=str, kw_only=True)
    chainID_idx: int = field(converter=int, kw_only=True)
    resSeq_idx: int = field(converter=int, kw_only=True)

    def get_mda_sel(self, backbone_only=True) -> str:
        return f"segid {self.chainID} and resnum {self.resSeq}" + (" and backbone" if backbone_only else "")

    def to_str(self) -> str:
        """
        Gives alternate string representation of the mutation
        Returns
        -------
        mutation string: str
        """
        return f"{self.chainID}:{self.old_aa}{self.resSeq}{self.new_aa}"

    @classmethod
    def from_site(cls, site: Site, *, old_aa: str, new_aa: str) -> "Mutation":
        return cls(
            chainID=site.chainID,
            resSeq=site.resSeq,
            old_aa=old_aa,
            new_aa=new_aa,
            chainID_idx=site.idx_chain,
            resSeq_idx=site.idx_residue)
