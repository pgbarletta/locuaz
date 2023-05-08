from attrs import define, field


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
