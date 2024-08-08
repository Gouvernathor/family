from typing import Callable

class StoringDefaultedDictMixin[K, V](dict[K, V]):
    __slots__ = ()

    def __missing__(self, key, /) -> V:
        self[key] = value = super().__missing__(key) # type: ignore
        return value

class ConstantMissingDict[K, V](dict[K, V]):
    __slots__ = ("constant",)

    def __init__(self, constant: V = None, /, *args, **kwargs):
        self.constant = constant
        super().__init__(*args, **kwargs)

    def __missing__(self, _key, /) -> V:
        return self.constant

class FactoryDict[K, V](dict[K, V]):
    __slots__ = ("factory",)

    def __init__(self, factory: Callable[[K], V], /, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.factory = factory

    def __missing__(self, key, /) -> V:
        return self.factory(key)

class StoringFactoryDict[K, V](StoringDefaultedDictMixin[K, V], FactoryDict[K, V]):
    __slots__ = ()
