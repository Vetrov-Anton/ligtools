import inspect
import hashlib
from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, field_serializer


def _is_settable(member):
    return isinstance(member, property) and member.fset is not None


def _to_immutable(obj):
    if isinstance(obj, np.ndarray):
        obj = obj.tolist()
    if isinstance(obj, set):
        return frozenset(obj)
    elif isinstance(obj, list):
        return tuple(_to_immutable(x) for x in obj)
    elif isinstance(obj, dict):
        return tuple((k, _to_immutable(v)) for k, v in sorted(obj.items()))
    return obj


class Model(BaseModel):
    """Base class that all classes should subclass.
    """

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        validate_assignment=True,
    )

    @field_serializer('*', mode='wrap', when_used='json', check_fields=False)
    def _serialize_arrays(self, value, handler):
        if isinstance(value, np.ndarray):
            return value.tolist()
        return handler(value)

    @property
    def _clsname(self):
        return type(self).__name__

    def __init__(__pydantic_self__, **data: Any) -> None:  # lgtm[py/not-named-self]
        __pydantic_self__.__pre_init__(**data)  # lgtm[py/init-calls-subclass]
        super().__init__(**data)
        __pydantic_self__.__post_init__(**data)  # lgtm[py/init-calls-subclass]

    def __pre_init__(self, **kwargs):
        pass

    def __post_init__(self, **kwargs):
        pass

    def __setattr__(self, attr, value):
        try:
            super().__setattr__(attr, value)
        except ValueError as e:
            setters = inspect.getmembers(self.__class__, predicate=_is_settable)
            for propname, _ in setters:
                if propname == attr:
                    return object.__setattr__(self, propname, value)
            raise e

    def get_hash(self):
        mash = hashlib.sha1()
        mash.update(self.model_dump_json().encode("utf-8"))
        return mash.hexdigest()

    def __hash__(self):
        return hash(self.get_hash())

    def __eq__(self, other):
        try:
            other_hash = hash(other)
        except TypeError:
            other_hash = hash(_to_immutable(other))
        return hash(self) == other_hash
