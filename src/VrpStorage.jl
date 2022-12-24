struct RCSPState
    pricing::Ptr{Cvoid}
    state::Ptr{Cvoid}
end

mutable struct VrpNodeInfoUnit <: Coluna.ColunaBase.AbstractNewStorageUnit 
    rcsp_states::Vector{RCSPState}
    last_rcost_fix_gap::Float64
end

Coluna.ColunaBase.new_storage_unit(::Type{VrpNodeInfoUnit}, _) =
    VrpNodeInfoUnit(RCSPState[], Inf)

struct VrpNodeInfo <: Coluna.ColunaBase.AbstractNewRecord
    rcsp_states::Vector{RCSPState}
    last_rcost_fix_gap::Float64
end

Coluna.ColunaBase.record_type(::Type{VrpNodeInfoUnit}) = VrpNodeInfo
Coluna.ColunaBase.storage_unit_type(::Type{VrpNodeInfo}) = VrpNodeInfoUnit

struct VrpNodeInfoKey <: Coluna.Algorithm.AbstractStorageUnitKey end

Coluna.Algorithm.key_from_storage_unit_type(::Type{VrpNodeInfoUnit}) = VrpNodeInfoKey()
Coluna.Algorithm.record_type_from_key(::VrpNodeInfoKey) = VrpNodeInfo

function Coluna.ColunaBase.new_record(
    ::Type{VrpNodeInfo}, id::Int, form::Coluna.MathProg.Formulation, unit::VrpNodeInfoUnit
)
    # @info "In new_record $(unit.rcsp_states)"
    return VrpNodeInfo(unit.rcsp_states, unit.last_rcost_fix_gap)
end

function Coluna.ColunaBase.restore_from_record!(
    ::Coluna.MathProg.Formulation, unit::VrpNodeInfoUnit, record::VrpNodeInfo
)
    # @info "In restore_from_record! $(record.rcsp_states)"
    for rcsp_state in record.rcsp_states
        restore_rcsp_state(rcsp_state.pricing, rcsp_state.state)
    end
    unit.rcsp_states = record.rcsp_states
    unit.last_rcost_fix_gap = record.last_rcost_fix_gap
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
Coluna.Algorithm.ismanager(::RedCostFixAndEnumAlgorithm) = false

function Coluna.Algorithm.get_units_usage(
    ::RedCostFixAndEnumAlgorithm, reform::Coluna.MathProg.Reformulation
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
