# Python stubs generated by omniidl from ..\..\..\..\..\idl\COS\RDITestTypes.idl
# DO NOT EDIT THIS FILE!

import omniORB, _omnipy
from omniORB import CORBA, PortableServer
_0_CORBA = CORBA


_omnipy.checkVersion(4,2, __file__, 1)

try:
    property
except NameError:
    def property(*args):
        return None


#
# Start of module "RDITestTypes"
#
__name__ = "RDITestTypes"
_0_RDITestTypes = omniORB.openModule("RDITestTypes", r"..\..\..\..\..\idl\COS\RDITestTypes.idl")
_0_RDITestTypes__POA = omniORB.openModule("RDITestTypes__POA", r"..\..\..\..\..\idl\COS\RDITestTypes.idl")


# typedef ... StringArrayFive
class StringArrayFive:
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StringArrayFive:1.0"
    def __init__(self, *args, **kw):
        raise RuntimeError("Cannot construct objects of this type.")
_0_RDITestTypes.StringArrayFive = StringArrayFive
_0_RDITestTypes._d_StringArrayFive  = (omniORB.tcInternal.tv_array, (omniORB.tcInternal.tv_string,0), 5)
_0_RDITestTypes._ad_StringArrayFive = (omniORB.tcInternal.tv_alias, StringArrayFive._NP_RepositoryId, "StringArrayFive", (omniORB.tcInternal.tv_array, (omniORB.tcInternal.tv_string,0), 5))
_0_RDITestTypes._tc_StringArrayFive = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._ad_StringArrayFive)
omniORB.registerType(StringArrayFive._NP_RepositoryId, _0_RDITestTypes._ad_StringArrayFive, _0_RDITestTypes._tc_StringArrayFive)
del StringArrayFive

# typedef ... StringArrayTen
class StringArrayTen:
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StringArrayTen:1.0"
    def __init__(self, *args, **kw):
        raise RuntimeError("Cannot construct objects of this type.")
_0_RDITestTypes.StringArrayTen = StringArrayTen
_0_RDITestTypes._d_StringArrayTen  = (omniORB.tcInternal.tv_array, (omniORB.tcInternal.tv_string,0), 10)
_0_RDITestTypes._ad_StringArrayTen = (omniORB.tcInternal.tv_alias, StringArrayTen._NP_RepositoryId, "StringArrayTen", (omniORB.tcInternal.tv_array, (omniORB.tcInternal.tv_string,0), 10))
_0_RDITestTypes._tc_StringArrayTen = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._ad_StringArrayTen)
omniORB.registerType(StringArrayTen._NP_RepositoryId, _0_RDITestTypes._ad_StringArrayTen, _0_RDITestTypes._tc_StringArrayTen)
del StringArrayTen

# enum UnionSwitch
_0_RDITestTypes.a = omniORB.EnumItem("a", 0)
_0_RDITestTypes.b = omniORB.EnumItem("b", 1)
_0_RDITestTypes.c = omniORB.EnumItem("c", 2)
_0_RDITestTypes.d = omniORB.EnumItem("d", 3)
_0_RDITestTypes.e = omniORB.EnumItem("e", 4)
_0_RDITestTypes.UnionSwitch = omniORB.Enum("IDL:research.att.com/RDITestTypes/UnionSwitch:1.0", (_0_RDITestTypes.a, _0_RDITestTypes.b, _0_RDITestTypes.c, _0_RDITestTypes.d, _0_RDITestTypes.e,))

_0_RDITestTypes._d_UnionSwitch  = (omniORB.tcInternal.tv_enum, _0_RDITestTypes.UnionSwitch._NP_RepositoryId, "UnionSwitch", _0_RDITestTypes.UnionSwitch._items)
_0_RDITestTypes._tc_UnionSwitch = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_UnionSwitch)
omniORB.registerType(_0_RDITestTypes.UnionSwitch._NP_RepositoryId, _0_RDITestTypes._d_UnionSwitch, _0_RDITestTypes._tc_UnionSwitch)

# union UnionType
_0_RDITestTypes.UnionType = omniORB.newEmptyClass()
class UnionType (omniORB.Union):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/UnionType:1.0"

_0_RDITestTypes.UnionType = UnionType

UnionType._m_to_d = {"aLong": _0_RDITestTypes.a, "bString": _0_RDITestTypes.b, "cShort": _0_RDITestTypes.c, "dArray": _0_RDITestTypes.d}
UnionType._d_to_m = {_0_RDITestTypes.a: "aLong", _0_RDITestTypes.b: "bString", _0_RDITestTypes.c: "cShort", _0_RDITestTypes.d: "dArray"}
UnionType._def_m  = "defaultBoolean"
UnionType._def_d  = _0_RDITestTypes.e

_0_RDITestTypes._m_UnionType  = ((_0_RDITestTypes.a, "aLong", omniORB.tcInternal.tv_long), (_0_RDITestTypes.b, "bString", (omniORB.tcInternal.tv_string,0)), (_0_RDITestTypes.c, "cShort", omniORB.tcInternal.tv_short), (_0_RDITestTypes.d, "dArray", omniORB.typeMapping["IDL:research.att.com/RDITestTypes/StringArrayFive:1.0"]), (_0_RDITestTypes.e, "defaultBoolean", omniORB.tcInternal.tv_boolean),)
_0_RDITestTypes._d_UnionType  = (omniORB.tcInternal.tv_union, UnionType, UnionType._NP_RepositoryId, "UnionType", omniORB.typeMapping["IDL:research.att.com/RDITestTypes/UnionSwitch:1.0"], 4, _0_RDITestTypes._m_UnionType, _0_RDITestTypes._m_UnionType[4], {_0_RDITestTypes.a: _0_RDITestTypes._m_UnionType[0], _0_RDITestTypes.b: _0_RDITestTypes._m_UnionType[1], _0_RDITestTypes.c: _0_RDITestTypes._m_UnionType[2], _0_RDITestTypes.d: _0_RDITestTypes._m_UnionType[3]})
_0_RDITestTypes._tc_UnionType = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_UnionType)
omniORB.registerType(UnionType._NP_RepositoryId, _0_RDITestTypes._d_UnionType, _0_RDITestTypes._tc_UnionType)
del UnionType

# typedef ... StringSeq
class StringSeq:
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StringSeq:1.0"
    def __init__(self, *args, **kw):
        raise RuntimeError("Cannot construct objects of this type.")
_0_RDITestTypes.StringSeq = StringSeq
_0_RDITestTypes._d_StringSeq  = (omniORB.tcInternal.tv_sequence, (omniORB.tcInternal.tv_string,0), 0)
_0_RDITestTypes._ad_StringSeq = (omniORB.tcInternal.tv_alias, StringSeq._NP_RepositoryId, "StringSeq", (omniORB.tcInternal.tv_sequence, (omniORB.tcInternal.tv_string,0), 0))
_0_RDITestTypes._tc_StringSeq = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._ad_StringSeq)
omniORB.registerType(StringSeq._NP_RepositoryId, _0_RDITestTypes._ad_StringSeq, _0_RDITestTypes._tc_StringSeq)
del StringSeq

# typedef ... DoubleSeq
class DoubleSeq:
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/DoubleSeq:1.0"
    def __init__(self, *args, **kw):
        raise RuntimeError("Cannot construct objects of this type.")
_0_RDITestTypes.DoubleSeq = DoubleSeq
_0_RDITestTypes._d_DoubleSeq  = (omniORB.tcInternal.tv_sequence, omniORB.tcInternal.tv_double, 0)
_0_RDITestTypes._ad_DoubleSeq = (omniORB.tcInternal.tv_alias, DoubleSeq._NP_RepositoryId, "DoubleSeq", (omniORB.tcInternal.tv_sequence, omniORB.tcInternal.tv_double, 0))
_0_RDITestTypes._tc_DoubleSeq = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._ad_DoubleSeq)
omniORB.registerType(DoubleSeq._NP_RepositoryId, _0_RDITestTypes._ad_DoubleSeq, _0_RDITestTypes._tc_DoubleSeq)
del DoubleSeq

# union ExampleUnion1
_0_RDITestTypes.ExampleUnion1 = omniORB.newEmptyClass()
class ExampleUnion1 (omniORB.Union):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/ExampleUnion1:1.0"

_0_RDITestTypes.ExampleUnion1 = ExampleUnion1

ExampleUnion1._m_to_d = {"l": 1}
ExampleUnion1._d_to_m = {1: "l"}
ExampleUnion1._def_m  = "d"
ExampleUnion1._def_d  = 0

_0_RDITestTypes._m_ExampleUnion1  = ((1, "l", omniORB.tcInternal.tv_long), (0, "d", omniORB.tcInternal.tv_double),)
_0_RDITestTypes._d_ExampleUnion1  = (omniORB.tcInternal.tv_union, ExampleUnion1, ExampleUnion1._NP_RepositoryId, "ExampleUnion1", omniORB.tcInternal.tv_boolean, 1, _0_RDITestTypes._m_ExampleUnion1, _0_RDITestTypes._m_ExampleUnion1[1], {1: _0_RDITestTypes._m_ExampleUnion1[0]})
_0_RDITestTypes._tc_ExampleUnion1 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_ExampleUnion1)
omniORB.registerType(ExampleUnion1._NP_RepositoryId, _0_RDITestTypes._d_ExampleUnion1, _0_RDITestTypes._tc_ExampleUnion1)
del ExampleUnion1

# union ExampleUnion2
_0_RDITestTypes.ExampleUnion2 = omniORB.newEmptyClass()
class ExampleUnion2 (omniORB.Union):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/ExampleUnion2:1.0"

_0_RDITestTypes.ExampleUnion2 = ExampleUnion2

ExampleUnion2._m_to_d = {"l": 1, "d": 2}
ExampleUnion2._d_to_m = {1: "l", 2: "d"}
ExampleUnion2._def_m  = None
ExampleUnion2._def_d  = None

_0_RDITestTypes._m_ExampleUnion2  = ((1, "l", omniORB.tcInternal.tv_long), (2, "d", omniORB.tcInternal.tv_double),)
_0_RDITestTypes._d_ExampleUnion2  = (omniORB.tcInternal.tv_union, ExampleUnion2, ExampleUnion2._NP_RepositoryId, "ExampleUnion2", omniORB.tcInternal.tv_long, -1, _0_RDITestTypes._m_ExampleUnion2, None, {1: _0_RDITestTypes._m_ExampleUnion2[0], 2: _0_RDITestTypes._m_ExampleUnion2[1]})
_0_RDITestTypes._tc_ExampleUnion2 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_ExampleUnion2)
omniORB.registerType(ExampleUnion2._NP_RepositoryId, _0_RDITestTypes._d_ExampleUnion2, _0_RDITestTypes._tc_ExampleUnion2)
del ExampleUnion2

# union ExampleUnion3
_0_RDITestTypes.ExampleUnion3 = omniORB.newEmptyClass()
class ExampleUnion3 (omniORB.Union):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/ExampleUnion3:1.0"

_0_RDITestTypes.ExampleUnion3 = ExampleUnion3

ExampleUnion3._m_to_d = {"l": 1, "d": 0}
ExampleUnion3._d_to_m = {1: "l", 0: "d"}
ExampleUnion3._def_m  = None
ExampleUnion3._def_d  = None

_0_RDITestTypes._m_ExampleUnion3  = ((1, "l", omniORB.tcInternal.tv_long), (0, "d", omniORB.tcInternal.tv_double),)
_0_RDITestTypes._d_ExampleUnion3  = (omniORB.tcInternal.tv_union, ExampleUnion3, ExampleUnion3._NP_RepositoryId, "ExampleUnion3", omniORB.tcInternal.tv_boolean, -1, _0_RDITestTypes._m_ExampleUnion3, None, {1: _0_RDITestTypes._m_ExampleUnion3[0], 0: _0_RDITestTypes._m_ExampleUnion3[1]})
_0_RDITestTypes._tc_ExampleUnion3 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_ExampleUnion3)
omniORB.registerType(ExampleUnion3._NP_RepositoryId, _0_RDITestTypes._d_ExampleUnion3, _0_RDITestTypes._tc_ExampleUnion3)
del ExampleUnion3

# struct StructExample1
_0_RDITestTypes.StructExample1 = omniORB.newEmptyClass()
class StructExample1 (omniORB.StructBase):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StructExample1:1.0"

    def __init__(self, d):
        self.d = d

_0_RDITestTypes.StructExample1 = StructExample1
_0_RDITestTypes._d_StructExample1  = (omniORB.tcInternal.tv_struct, StructExample1, StructExample1._NP_RepositoryId, "StructExample1", "d", omniORB.tcInternal.tv_double)
_0_RDITestTypes._tc_StructExample1 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_StructExample1)
omniORB.registerType(StructExample1._NP_RepositoryId, _0_RDITestTypes._d_StructExample1, _0_RDITestTypes._tc_StructExample1)
del StructExample1

# struct StructExample2
_0_RDITestTypes.StructExample2 = omniORB.newEmptyClass()
class StructExample2 (omniORB.StructBase):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StructExample2:1.0"

    def __init__(self, event_name, d):
        self.event_name = event_name
        self.d = d

_0_RDITestTypes.StructExample2 = StructExample2
_0_RDITestTypes._d_StructExample2  = (omniORB.tcInternal.tv_struct, StructExample2, StructExample2._NP_RepositoryId, "StructExample2", "event_name", (omniORB.tcInternal.tv_string,0), "d", omniORB.tcInternal.tv_double)
_0_RDITestTypes._tc_StructExample2 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_StructExample2)
omniORB.registerType(StructExample2._NP_RepositoryId, _0_RDITestTypes._d_StructExample2, _0_RDITestTypes._tc_StructExample2)
del StructExample2

# struct StructExample3
_0_RDITestTypes.StructExample3 = omniORB.newEmptyClass()
class StructExample3 (omniORB.StructBase):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StructExample3:1.0"

    def __init__(self, domain_name, type_name, filterable_data, d):
        self.domain_name = domain_name
        self.type_name = type_name
        self.filterable_data = filterable_data
        self.d = d

_0_RDITestTypes.StructExample3 = StructExample3
_0_RDITestTypes._d_StructExample3  = (omniORB.tcInternal.tv_struct, StructExample3, StructExample3._NP_RepositoryId, "StructExample3", "domain_name", (omniORB.tcInternal.tv_string,0), "type_name", (omniORB.tcInternal.tv_string,0), "filterable_data", (omniORB.tcInternal.tv_string,0), "d", omniORB.tcInternal.tv_double)
_0_RDITestTypes._tc_StructExample3 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_StructExample3)
omniORB.registerType(StructExample3._NP_RepositoryId, _0_RDITestTypes._d_StructExample3, _0_RDITestTypes._tc_StructExample3)
del StructExample3

# struct StructExample4
_0_RDITestTypes.StructExample4 = omniORB.newEmptyClass()
class StructExample4 (omniORB.StructBase):
    _NP_RepositoryId = "IDL:research.att.com/RDITestTypes/StructExample4:1.0"

    def __init__(self, part1, part2, part3):
        self.part1 = part1
        self.part2 = part2
        self.part3 = part3

_0_RDITestTypes.StructExample4 = StructExample4
_0_RDITestTypes._d_StructExample4  = (omniORB.tcInternal.tv_struct, StructExample4, StructExample4._NP_RepositoryId, "StructExample4", "part1", omniORB.typeMapping["IDL:research.att.com/RDITestTypes/StructExample1:1.0"], "part2", omniORB.typeMapping["IDL:research.att.com/RDITestTypes/StructExample2:1.0"], "part3", omniORB.typeMapping["IDL:research.att.com/RDITestTypes/StructExample3:1.0"])
_0_RDITestTypes._tc_StructExample4 = omniORB.tcInternal.createTypeCode(_0_RDITestTypes._d_StructExample4)
omniORB.registerType(StructExample4._NP_RepositoryId, _0_RDITestTypes._d_StructExample4, _0_RDITestTypes._tc_StructExample4)
del StructExample4

#
# End of module "RDITestTypes"
#
__name__ = "RDITestTypes_idl"

_exported_modules = ( "RDITestTypes", )

# The end.
