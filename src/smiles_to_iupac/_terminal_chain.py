"""Parent-chain naming for a chalcogenyl terminal ('R-sulfanyl', 'R-disulfanyl', ...) hung on an acyclic carbon chain.
The terminal is always a compound prefix (P-16.5.1.1), so it is parenthesized even without a locant; the locant `1`
is still omitted on methane and on a monosubstituted two-carbon chain (P-14.3.4.2). A `bare_name` terminal such as
'disulfanyl' carries no nested prefix and is cited unparenthesized in those no-locant positions."""

from ._common import group_substituents, longest_chains, substituent_locant_set_and_citation
from ._multiplicative_text import enclose
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, substituents_for_chain_forced_compound_terminals


def _name_from_grouped(chain_length, grouped, bare_name):
    if chain_length == 1 and grouped:
        (name,) = grouped
        if name == bare_name:
            return name + alkane_name(chain_length)
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and total_subs == 1:
        (name,) = grouped
        display_name = name if name == bare_name else (enclose(name) if grouped[name]["compound"] else name)
        return display_name + alkane_name(chain_length)
    return format_substituent_prefixes(grouped) + alkane_name(chain_length)


def _candidate_key(chain_length, substituents, bare_name):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_grouped(chain_length, grouped, bare_name)
    return (-total_count, locant_set, citation_locants, name), name


def best_terminal_chain(full_graph, carbon_graph, terminals, bare_name=None, mol=None):
    """(key, chain, name) of the winning longest chain and direction; `terminals` is {atom_idx -> prefix name}."""
    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])

    best = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            substituents = substituents_for_chain_forced_compound_terminals(full_graph, candidate, terminals, mol=mol)
            key, name = _candidate_key(chain_length, substituents, bare_name)
            if best is None or key < best[0]:
                best = (key, candidate, name)
    return best
