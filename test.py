from family import *

class Person(FamilyMember):
    def __init__(self, name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.name!r})"

a = Person("Anthony")
b = Person("Betty")
a.add_parents(b)
c = Person("Charlie")
d = Person("Duncan")
b.add_spouses(c)
c.add_parents(d)
e = Person("Elizea")
b.add_parents(e)
f = Person("Francis")
g = Person("Greg")
b.add_children(f)
f.add_spouses(g)
h = Person("Helga")
g.add_parents(h)
print("\n".join(a.get_readable_relationship(p, r) for p, r in a.get_all_reachable_relationship_chains().items()))
