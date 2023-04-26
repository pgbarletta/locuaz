==========================================
Tutorial: running a simple optimization
==========================================

locuaz iteratively performs mutations on the selected residues of the binder, followed by minimization and
MD simulation. The trajectory is then scored with various scoring functions, to assess whether the mutation
results in higher or lower affinity. If the mutation results in better affinity, the mutation is accepted,
and the process is repeated. If the mutation does not significantly improve the affinity, then the mutants
are discarded and a new set of mutants are generated, based on the original complex(es).

