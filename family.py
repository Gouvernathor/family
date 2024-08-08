from __future__ import annotations

import enum
import itertools
from typing import Callable, ClassVar, Iterable, Mapping, Self

from util import StoringFactoryDict

class Relationship(enum.Flag):
    # technical
    NONE = 0
    SELF = 1

    # direct
    PARENT = enum.auto()
    CHILD = enum.auto()
    SPOUSE = enum.auto()

    # direct but special
    CLONE = enum.auto()
    CLONE_ORIGINAL = enum.auto()

    DIRECT_DISTANT_ANCESTOR = enum.auto()
    DIRECT_DISTANT_DESCENDANT = enum.auto()

    # indirect
    SIBLING = enum.auto()
    HALF_SIBLING = enum.auto()
    STEP_SIBLING = enum.auto()

    GRANDPARENT = enum.auto()
    GRANDCHILD = enum.auto()
    GREAT_GRANDPARENT = enum.auto()
    GREAT_GRANDCHILD = enum.auto()

    AUNT_UNCLE = enum.auto()
    NEPHEW_NIECE = enum.auto()
    FIRST_COUSIN = enum.auto()

    GRAND_AUNT_UNCLE = enum.auto()
    GRAND_NEPHEW_NIECE = enum.auto()

    PARENT_IN_LAW = enum.auto()
    CHILD_IN_LAW = enum.auto()
    # this can be full sibling or half sibling, but not step sibling
    SIBLING_IN_LAW = enum.auto()

    CO_CLONE = enum.auto()

    def __pos__(self) -> Self:
        return self

    def __neg__(self) -> Relationship:
        match self:
            case Relationship.NONE: return Relationship.NONE
            case Relationship.SELF: return Relationship.SELF
            case Relationship.PARENT: return Relationship.CHILD
            case Relationship.CHILD: return Relationship.PARENT
            case Relationship.SPOUSE: return Relationship.SPOUSE
            case Relationship.CLONE: return Relationship.CLONE_ORIGINAL
            case Relationship.CLONE_ORIGINAL: return Relationship.CLONE
            case Relationship.DIRECT_DISTANT_ANCESTOR: return Relationship.DIRECT_DISTANT_DESCENDANT
            case Relationship.DIRECT_DISTANT_DESCENDANT: return Relationship.DIRECT_DISTANT_ANCESTOR
            case Relationship.SIBLING: return Relationship.SIBLING
            case Relationship.HALF_SIBLING: return Relationship.HALF_SIBLING
            case Relationship.STEP_SIBLING: return Relationship.STEP_SIBLING
            case Relationship.GRANDPARENT: return Relationship.GRANDCHILD
            case Relationship.GRANDCHILD: return Relationship.GRANDPARENT
            case Relationship.GREAT_GRANDPARENT: return Relationship.GREAT_GRANDCHILD
            case Relationship.GREAT_GRANDCHILD: return Relationship.GREAT_GRANDPARENT
            case Relationship.AUNT_UNCLE: return Relationship.NEPHEW_NIECE
            case Relationship.NEPHEW_NIECE: return Relationship.AUNT_UNCLE
            case Relationship.FIRST_COUSIN: return Relationship.FIRST_COUSIN
            case Relationship.GRAND_AUNT_UNCLE: return Relationship.GRAND_NEPHEW_NIECE
            case Relationship.GRAND_NEPHEW_NIECE: return Relationship.GRAND_AUNT_UNCLE
            case Relationship.PARENT_IN_LAW: return Relationship.CHILD_IN_LAW
            case Relationship.CHILD_IN_LAW: return Relationship.PARENT_IN_LAW
            case Relationship.SIBLING_IN_LAW: return Relationship.SIBLING_IN_LAW
            case Relationship.CO_CLONE: return Relationship.CO_CLONE

        # if we arrive here, we have an union of relationships
        rv = Relationship.NONE
        for r in Relationship:
            if r in self:
                rv |= -r
        return rv

    _testers: ClassVar[dict[Relationship, Callable[[FamilyMember, FamilyMember], bool]]] = enum.nonmember({})
    def _register_tester(self, func, /) -> Callable[[FamilyMember, FamilyMember], bool]:
        """
        Should be used (on an appropriate method) at least once per relationship enum member,
        except for Relationship.NONE.
        The order matters : the first matching relation will be preferred between two persons.
        The callable returns whether the second one is [the tested relationship] of the first one.
        """
        Relationship._testers[self] = func
        return func

    def get_human_readable(self, /, *, article=False, preposition=False, unique=False) -> str:
        """
        `unique` affects the chosen article. If `article` is False, `unique` is ignored.
        """
        art = "a "
        prep = " of"

        match self:
            case Relationship.NONE: art, main, prep = "", "not related", " to"
            case Relationship.SELF: art, main, prep = "", "the same person", " as"

            case Relationship.PARENT: main = "parent"
            case Relationship.CHILD: main = "child"
            case Relationship.SPOUSE: main = "spouse"

            case Relationship.CLONE: main = "clone"
            case Relationship.CLONE_ORIGINAL: art, main = "some ", "source genetic material"

            case Relationship.DIRECT_DISTANT_ANCESTOR: main = "distant ancestor"
            case Relationship.DIRECT_DISTANT_DESCENDANT: main = "distant descendant"

            case Relationship.SIBLING: main = "sibling"
            case Relationship.HALF_SIBLING: main = "half sibling"
            case Relationship.STEP_SIBLING: main = "step sibling"

            case Relationship.GRANDPARENT: main = "grandparent"
            case Relationship.GRANDCHILD: main = "grandchild"
            case Relationship.GREAT_GRANDPARENT: main = "great grandparent"
            case Relationship.GREAT_GRANDCHILD: main = "great grandchild"

            case Relationship.AUNT_UNCLE: art, main = "an ", "aunt or uncle"
            case Relationship.NEPHEW_NIECE: main = "nephew or niece"
            case Relationship.FIRST_COUSIN: main = "first cousin"

            case Relationship.GRAND_AUNT_UNCLE: main = "grand aunt or uncle"
            case Relationship.GRAND_NEPHEW_NIECE: main = "grand nephew or niece"

            case Relationship.PARENT_IN_LAW: main = "parent in law"
            case Relationship.CHILD_IN_LAW: main = "child in law"
            case Relationship.SIBLING_IN_LAW: main = "sibling in law"

            case Relationship.CO_CLONE: main = "co-clone"

            case _: raise ValueError(f"Unknown relationship {self!r}")

        if unique and art:
            art = "the "
        if article:
            main = art + main
        if preposition:
            main = main + prep
        return main

    human_readable = property(get_human_readable)

class FamilyMember:
    """
    This is a Person class.
    """

    _parents: set[FamilyMember]
    _children: set[FamilyMember]
    _spouses: set[FamilyMember]

    _clones: set[FamilyMember]
    _clone_originals: set[FamilyMember]

    _distant_ancestors: set[FamilyMember]
    _distant_descendants: set[FamilyMember]

    def __init__(self):
        self._parents = set()
        self._children = set()
        self._spouses = set()

        self._clones = set()
        self._clone_originals = set()

        self._distant_ancestors = set()
        self._distant_descendants = set()

    @Relationship.SELF._register_tester
    def is_self(self, other: FamilyMember, /) -> bool:
        return self == other

    @property
    def parents(self) -> frozenset[FamilyMember]:
        return frozenset(self._parents)
    def add_parents(self, /, *parents: FamilyMember) -> None:
        self._parents.update(parents)
        for parent in parents:
            parent._children.add(self)
    def remove_parents(self, /, *parents: FamilyMember) -> None:
        self._parents.difference_update(parents)
        for parent in parents:
            parent._children.discard(self)
    @Relationship.PARENT._register_tester
    def is_parent(self, other: FamilyMember, /) -> bool:
        """
        These tester functions return whether *other* is the parent of *self*.
        Not the other way around.
        """
        return other in self._parents

    @property
    def children(self) -> frozenset[FamilyMember]:
        return frozenset(self._children)
    def add_children(self, /, *children: FamilyMember) -> None:
        self._children.update(children)
        for child in children:
            child._parents.add(self)
    def remove_children(self, /, *children: FamilyMember) -> None:
        self._children.difference_update(children)
        for child in children:
            child._parents.discard(self)
    @Relationship.CHILD._register_tester
    def is_child(self, other: FamilyMember, /) -> bool:
        return other in self._children

    @property
    def spouses(self) -> frozenset[FamilyMember]:
        return frozenset(self._spouses)
    def add_spouses(self, /, *spouses: FamilyMember) -> None:
        self._spouses.update(spouses)
        for spouse in spouses:
            spouse._spouses.add(self)
    def remove_spouses(self, /, *spouses: FamilyMember) -> None:
        self._spouses.difference_update(spouses)
        for spouse in spouses:
            spouse._spouses.discard(self)
    @Relationship.SPOUSE._register_tester
    def is_spouse(self, other: FamilyMember, /) -> bool:
        return other in self._spouses

    @property
    def clones(self) -> frozenset[FamilyMember]:
        return frozenset(self._clones)
    def add_clones(self, /, *clones: FamilyMember) -> None:
        self._clones.update(clones)
        for clone in clones:
            clone._clone_originals.add(self)
    def remove_clones(self, /, *clones: FamilyMember) -> None:
        self._clones.difference_update(clones)
        for clone in clones:
            clone._clone_originals.discard(self)
    @Relationship.CLONE._register_tester
    def is_clone(self, other: FamilyMember, /) -> bool:
        return other in self._clones

    @property
    def clone_originals(self) -> frozenset[FamilyMember]:
        return frozenset(self._clone_originals)
    def add_clone_originals(self, /, *clone_originals: FamilyMember) -> None:
        self._clone_originals.update(clone_originals)
        for clone_original in clone_originals:
            clone_original._clones.add(self)
    def remove_clone_originals(self, /, *clone_originals: FamilyMember) -> None:
        self._clone_originals.difference_update(clone_originals)
        for clone_original in clone_originals:
            clone_original._clones.discard(self)
    @Relationship.CLONE_ORIGINAL._register_tester
    def is_clone_original(self, other: FamilyMember, /) -> bool:
        return other in self._clone_originals

    @property
    def distant_ancestors(self) -> frozenset[FamilyMember]:
        return frozenset(self._distant_ancestors)
    def add_distant_ancestors(self, /, *distant_ancestors: FamilyMember) -> None:
        self._distant_ancestors.update(distant_ancestors)
        for distant_ancestor in distant_ancestors:
            distant_ancestor._distant_descendants.add(self)
    def remove_distant_ancestors(self, /, *distant_ancestors: FamilyMember) -> None:
        self._distant_ancestors.difference_update(distant_ancestors)
        for distant_ancestor in distant_ancestors:
            distant_ancestor._distant_descendants.discard(self)
    @Relationship.DIRECT_DISTANT_ANCESTOR._register_tester
    def is_direct_distant_ancestor(self, other: FamilyMember, /) -> bool:
        return other in self._distant_ancestors

    @property
    def distant_descendants(self) -> frozenset[FamilyMember]:
        return frozenset(self._distant_descendants)
    def add_distant_descendants(self, /, *distant_descendants: FamilyMember) -> None:
        self._distant_descendants.update(distant_descendants)
        for distant_descendant in distant_descendants:
            distant_descendant._distant_ancestors.add(self)
    def remove_distant_descendants(self, /, *distant_descendants: FamilyMember) -> None:
        self._distant_descendants.difference_update(distant_descendants)
        for distant_descendant in distant_descendants:
            distant_descendant._distant_ancestors.discard(self)
    @Relationship.DIRECT_DISTANT_DESCENDANT._register_tester
    def is_direct_distant_descendant(self, other: FamilyMember, /) -> bool:
        return other in self._distant_descendants

    # siblings of different kinds
    def _is_sibling_or_half(self, other: FamilyMember, half: bool|None = None) -> bool:
        """
        If half is None, returns True if it's a full sibling or a half sibling.
        """
        self_parents = self._parents
        other_parents = other._parents

        if len(self_parents) < 2:
            if half is False:
                # not enough parents to be considered full siblings
                return False
        elif self_parents == other_parents:
            # same 2+ parents, full siblings, not half siblings
            return half is not True

        if half is not False:
            # common parents, and not same 2+ parents, half siblings
            return bool(self_parents & other_parents)
        else:
            # 2+ parents but not the same, not full siblings
            return False

    @Relationship.SIBLING._register_tester
    def is_sibling(self, other: FamilyMember, /) -> bool:
        """
        A sibling is someone who has the exact same parents as self,
        if these parents are two or more.
        """
        return self._is_sibling_or_half(other, half=False)

    @Relationship.HALF_SIBLING._register_tester
    def is_half_sibling(self, other: FamilyMember, /) -> bool:
        return self._is_sibling_or_half(other, half=True)

    @Relationship.STEP_SIBLING._register_tester
    def is_step_sibling(self, other: FamilyMember, /) -> bool:
        if self._is_sibling_or_half(other, half=None):
            return False

        for parent in self._parents:
            if other._parents & parent._spouses:
                return True
        return False

    # grand parents and grand children
    @Relationship.GRANDPARENT._register_tester
    def is_grandparent(self, other: FamilyMember, /) -> bool:
        for parent in self._parents:
            if other in parent._parents:
                return True
        return False

    @Relationship.GRANDCHILD._register_tester
    def is_grandchild(self, other: FamilyMember, /) -> bool:
        return other.is_grandparent(self)

    @Relationship.GREAT_GRANDPARENT._register_tester
    def is_great_grandparent(self, other: FamilyMember, /) -> bool:
        for parent in self._parents:
            for grandparent in parent._parents:
                if other in grandparent._parents:
                    return True
        return False

    @Relationship.GREAT_GRANDCHILD._register_tester
    def is_great_grandchild(self, other: FamilyMember, /) -> bool:
        return other.is_great_grandparent(self)

    # aunts, uncles, nephews, nieces, cousins
    @Relationship.AUNT_UNCLE._register_tester
    def is_aunt_uncle(self, other: FamilyMember, /) -> bool:
        """
        Returns True if other is an aunt or uncle of self.
        That's the way these relationships work.
        This method should not return True with a normal parent of self,
        but it should return True with a parent that is a sibling of another parent of self.
        """
        for parent in self._parents:
            if parent == other:
                continue
            for grandparent in parent._parents:
                if other in grandparent._children:
                    return True
        return False

    @Relationship.NEPHEW_NIECE._register_tester
    def is_nephew_niece(self, other: FamilyMember, /) -> bool:
        return other.is_aunt_uncle(self)

    @Relationship.FIRST_COUSIN._register_tester
    def is_first_cousin(self, other: FamilyMember, /) -> bool:
        for parent in self._parents:
            for grandparent in parent._parents:
                for auntuncle in grandparent._children:
                    if auntuncle == parent:
                        continue
                    if other in auntuncle._children:
                        return True
        return False

    @Relationship.GRAND_AUNT_UNCLE._register_tester
    def is_grand_aunt_uncle(self, other: FamilyMember, /) -> bool:
        for parent in self._parents:
            if parent.is_aunt_uncle(other):
                return True
        return False

    @Relationship.GRAND_NEPHEW_NIECE._register_tester
    def is_grand_nephew_niece(self, other: FamilyMember, /) -> bool:
        return other.is_grand_aunt_uncle(self)

    # in-laws
    @Relationship.PARENT_IN_LAW._register_tester
    def is_parent_in_law(self, other: FamilyMember, /) -> bool:
        for spouse in self._spouses:
            if other in spouse._parents:
                return True
        return False

    @Relationship.CHILD_IN_LAW._register_tester
    def is_child_in_law(self, other: FamilyMember, /) -> bool:
        for child in self._children:
            if other in child._spouses:
                return True
        return False

    @Relationship.SIBLING_IN_LAW._register_tester
    def is_sibling_in_law(self, other: FamilyMember, /) -> bool:
        # spouse of sibling
        for parent in self._parents:
            for sibling in parent._children:
                if sibling == self:
                    continue
                if other in sibling._spouses:
                    return True

        # sibling of spouse
        for spouse in self._spouses:
            if spouse._is_sibling_or_half(other, half=None):
                return True

        return False

    # clones
    @Relationship.CO_CLONE._register_tester
    def is_co_clone(self, other: FamilyMember, /) -> bool:
        for original in self._clone_originals:
            if other in original._clones:
                return True
        return False


    def is_of_relationship(self, relationship: Relationship, other: FamilyMember, /) -> bool:
        """
        Relationship unions are not supported, because it's not clear whether to return
        if the members share *any* of the relationships, or *all* of them.
        """
        tester = Relationship._testers.get(relationship, None)
        if tester is not None:
            return tester(self, other)

        # special case
        if relationship == Relationship.NONE:
            return not any(self.is_of_relationship(r, other) for r in Relationship)

        raise ValueError(f"Unknown relationship {relationship!r}")

    def copy(self, /) -> Self:
        """
        Not the same as cloning.
        Cloning is single-directional, copying is not.
        The copied instance will test the same as the original for any relationship with any other instance.
        """
        new = getattr(super(), "__copy__", type(self))()
        new.add_parents(*self._parents)
        new.add_children(*self._children)
        new.add_spouses(*self._spouses)
        new.add_clones(*self._clones)
        new.add_clone_originals(*self._clone_originals)
        new.add_distant_ancestors(*self._distant_ancestors)
        new.add_distant_descendants(*self._distant_descendants)
        return new

    __copy__ = copy


    def get_nominal_relationship(self, other: FamilyMember, /) -> Relationship:
        """
        Returns the first matching relationship, in priority order,
        or Relationship.NONE if no relationship is found.
        """
        for relationship, tester in Relationship._testers.items():
            if tester(self, other):
                return relationship

        return Relationship.NONE

    def get_full_relationship(self, other: FamilyMember, /) -> Relationship:
        """
        Uses the union of enum flags.
        Much slower than get_nominal_relationship.
        """
        rv = Relationship.NONE
        for relationship, tester in Relationship._testers.items():
            if tester(self, other):
                rv |= relationship
        return rv


    def make_familytree(self, /, *, full: bool = False) -> Mapping[FamilyMember, Relationship]:
        """
        Returns all the full relationships between all persons reachable from this person.
        The person itself is included.
        """
        if full:
            meth = self.get_full_relationship
        else:
            meth = self.get_nominal_relationship

        memberset = set()
        # the maximum separation matching a listed relationship is 4, grand uncle/niece
        # with a complimentary padding of 1
        self.get_reachable_family_members(memberset, 5)

        return {member: relationship for member in memberset if (relationship:=meth(member))}

    def _iter_all_close_relatives(self):
        yield from self._parents
        yield from self._children
        yield from self._spouses
        yield from self._clones
        yield from self._clone_originals
        yield from self._distant_ancestors
        yield from self._distant_descendants

    def _iter_all_close_relatives_with_relationship(self):
        repeat = itertools.repeat
        yield from zip(self._parents, repeat(Relationship.PARENT))
        yield from zip(self._children, repeat(Relationship.CHILD))
        yield from zip(self._spouses, repeat(Relationship.SPOUSE))
        yield from zip(self._clones, repeat(Relationship.CLONE))
        yield from zip(self._clone_originals, repeat(Relationship.CLONE_ORIGINAL))
        yield from zip(self._distant_ancestors, repeat(Relationship.DIRECT_DISTANT_ANCESTOR))
        yield from zip(self._distant_descendants, repeat(Relationship.DIRECT_DISTANT_DESCENDANT))

    def get_reachable_family_members(self, rv: set, /, distance: int = float('inf')) -> None: # type: ignore
        """
        Mutates a set in-place.
        Fills it with all people reachable from self within distance hops of base relationships.
        """
        if self in rv:
            return

        rv.add(self)

        if not distance:
            return

        for person in set(self._iter_all_close_relatives())-rv:
            person.get_reachable_family_members(rv, distance-1)


    def get_all_reachable_relationship_chains(self, /) -> dict[FamilyMember, tuple[Relationship, ...]]:

        def key_function(rl: tuple[Relationship, ...], /) -> tuple[int, int]:
            return (len(rl), sum(r.value for r in rl))

        def neighbord_initial_finder(node: FamilyMember, /) -> list[tuple[FamilyMember, tuple[Relationship]]]:
            return [(other, (relationship,)) for other, relationship in node._iter_all_close_relatives_with_relationship()]

        floating_distance_cache: dict[tuple[Relationship, ...], tuple[int, int]] = StoringFactoryDict(key_function)
        new_neighbors: dict[FamilyMember, list[tuple[FamilyMember, tuple[Relationship]]]] = StoringFactoryDict(neighbord_initial_finder)
        chains: dict[FamilyMember, tuple[Relationship, ...]] = {self: ()}
        unvisited: set[FamilyMember] = {self}
        while unvisited:
            current = min(unvisited, key=(lambda member: floating_distance_cache[chains[member]]))
            unvisited.remove(current)
            current_chain = chains[current]
            newly_unvisited = set()
            for relative, relationship_chain in new_neighbors[current]:
                new_chain = current_chain + relationship_chain
                old_chain = chains.get(relative)
                if old_chain is None or floating_distance_cache[new_chain] < floating_distance_cache[old_chain]:
                    chains[relative] = new_chain
                    newly_unvisited.add(relative)
            del new_neighbors[current][:]
            for visited in chains.keys()-newly_unvisited:
                for newcomer in newly_unvisited:
                    relationship = visited.get_nominal_relationship(newcomer)
                    if relationship is not Relationship.NONE:
                        new_chain = chains[visited] + (relationship,)
                        if floating_distance_cache[new_chain] < floating_distance_cache[chains[newcomer]]:
                            chains[newcomer] = new_chain
                            unvisited.add(visited)
                            new_neighbors[visited].append((newcomer, (relationship,)))
            unvisited.update(newly_unvisited)

        del chains[self]
        return chains


    def get_readable_relationship(self, other: FamilyMember, relationship: Relationship|Iterable[Relationship]|None = None, /) -> str:
        if relationship is None:
            relationship = self.get_nominal_relationship(other)
        relationship = tuple(relationship)

        if relationship == (Relationship.NONE,):
            return f"{other} is not (closely) related to {self}."
        elif relationship == (Relationship.SELF,):
            return f"{other} is the same person as {self}."

        return f"{other} is " + " ".join(r.get_human_readable(article=True, preposition=True) for r in relationship[::-1]) + f" {self}."
