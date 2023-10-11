from typing import Final

from attrs import define, field


# from siteselector import Site

@define(frozen=True)
class Site:
    idx_chain: Final[int] = field(converter=int, kw_only=True)
    chainID: Final[str] = field(converter=str, kw_only=True)
    idx_residue: Final[int] = field(converter=int, kw_only=True)
    resSeq: Final[int] = field(converter=int, kw_only=True)


@define(frozen=True)
class Mutation:
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
