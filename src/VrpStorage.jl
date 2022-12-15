mutable struct VrpNodeInfoUnit <: Coluna.ColunaBase.AbstractNewStorageUnit 
    rcsp_state::Vector{Ptr{Cvoid}}
end

Coluna.ColunaBase.new_storage_unit(::Type{VrpNodeInfoUnit}, _) = VrpNodeInfoUnit(111)

struct VrpNodeInfo <: Coluna.ColunaBase.AbstractNewRecord
    rcsp_state::Vector{Ptr{Cvoid}}
end

Coluna.ColunaBase.record_type(::Type{VrpNodeInfoUnit}) = VrpNodeInfo
Coluna.ColunaBase.storage_unit_type(::Type{VrpNodeInfo}) = VrpNodeInfoUnit

struct VrpNodeInfoKey <: Coluna.Algorithm.AbstractStorageUnitKey end

Coluna.Algorithm.key_from_storage_unit_type(::Type{VrpNodeInfoUnit}) = VrpNodeInfoKey()
Coluna.Algorithm.record_type_from_key(::VrpNodeInfoKey) = VrpNodeInfo

function Coluna.ColunaBase.new_record(
    ::Type{VrpNodeInfo}, id::Int, form::Coluna.MathProg.Formulation, unit::VrpNodeInfoUnit
)
    return VrpNodeInfo(unit.rcsp_state)
end

function Coluna.ColunaBase.restore_from_record!(
    ::Coluna.MathProg.Formulation, unit::VrpNodeInfoUnit, record::VrpNodeInfo
)
    unit.rcsp_state = record.rcsp_state
    return
end

function Coluna.Algorithm.get_branching_candidate_units_usage(
    ::Coluna.Algorithm.SingleVarBranchingCandidate, reform
)
    units_to_restore = Coluna.Algorithm.UnitsUsage()
    push!(
        units_to_restore.units_used,
        (Coluna.MathProg.getmaster(reform), Coluna.Algorithm.MasterBranchConstrsUnit)
    )
    push!(units_to_restore.units_used, (Coluna.MathProg.getmaster(reform), VrpNodeInfoUnit))
    return units_to_restore
end

Coluna.Algorithm.ismanager(::Coluna.Algorithm.BeforeCutGenAlgo) = false
Coluna.Algorithm.ismanager(::ImproveRelaxationAlgo) = false

function Coluna.Algorithm.get_units_usage(
    ::ImproveRelaxationAlgo, reform::Coluna.MathProg.Reformulation
) 
    units_usage = Tuple{
        Coluna.MathProg.AbstractModel,Coluna.ColunaBase.UnitType,Coluna.ColunaBase.UnitPermission
    }[]
    master = Coluna.MathProg.getmaster(reform)
    push!(units_usage, (master, VrpNodeInfoUnit, Coluna.ColunaBase.READ_AND_WRITE))
    return units_usage
end

function Coluna.Algorithm.get_child_algorithms(
    algo::Coluna.Algorithm.BeforeCutGenAlgo, reform::Coluna.MathProg.Reformulation
)
    child_algos = Tuple{Coluna.Algorithm.AbstractAlgorithm, Coluna.MathProg.AbstractModel}[]
    push!(child_algos, (algo.algorithm, reform))
    return child_algos
end
