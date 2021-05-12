using Catalyst, ModelingToolkit, Symbolics

# module SBML
# struct Model end
# struct Reaction end
# struct Compartment end
# struct Species end
# end

""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(model::Model; kwargs...)  # Todo: requires unique parameters (i.e. SBML must have been imported with localParameter promotion in libSBML)
    model = make_extensive(model)
    # @info "# REACTIONS" length(model.reactions)
    model = expand_reversible(model)
    rxs = mtk_reactions(model)
    @info "RXS" rxs
    t = Catalyst.DEFAULT_IV
    species = [create_var(k) for k in keys(model.species)]
    params = vcat([create_param(k) for k in keys(model.parameters)], [create_param(k) for k in keys(model.compartments)])
    ReactionSystem(rxs, t, species, params; kwargs...)
end

""" ReactionSystem constructor """
function ModelingToolkit.ReactionSystem(sbmlfile::String; kwargs...)
    model = readSBML(sbmlfile)
    ReactionSystem(model; kwargs...)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(model::Model; kwargs...)
    rs = ReactionSystem(model; kwargs...)
    model = make_extensive(model)  # PL: consider making `make_extensive!` to avoid duplicate calling in ReactionSystem and here
    u0map = get_u0(model)
    parammap = get_paramap(model)
    defaults = Dict(vcat(u0map, parammap))
    convert(ODESystem, rs, defaults=defaults)
end

""" ODESystem constructor """
function ModelingToolkit.ODESystem(sbmlfile::String; kwargs...)
    model = readSBML(sbmlfile)
    ODESystem(model; kwargs...)
end

""" ODEProblem constructor """
function ModelingToolkit.ODEProblem(model::Model, tspan;kwargs...)  # PL: Todo: add u0 and parameters argument
    odesys = ODESystem(model;kwargs...)
    ODEProblem(odesys, [], tspan)
end

""" ODEProblem constructor """
function ModelingToolkit.ODEProblem(sbmlfile::String, tspan;kwargs...)  # PL: Todo: add u0 and parameters argument
    odesys = ODESystem(sbmlfile;kwargs...)
    ODEProblem(odesys, [], tspan)
end

""" Convert intensive to extensive expressions """
function make_extensive(model)
    model = to_initial_amounts(model)
    # model = to_extensive_math(model)
    model  # Todo: For spevies with `hOSU=false` multiply all occurences in mathematical expressions by compartment size.
           # Also convert species initialConcentrations to initialAmounts
end

""" Convert initial_concentration to initial_amount """
function to_initial_amounts(model::Model)  # Test written
    @info "INITIAL AMOUNTS"
    model = deepcopy(model)
    for specie in values(model.species)
        if isequal(specie.initial_amount, nothing)
            compartment = model.compartments[specie.compartment]
            specie.initial_amount = (specie.initial_concentration[1] * compartment.size, "")
            specie.initial_concentration = nothing
        end
        end
    model
end

""" Convert intensive to extensive mathematical expression """
function to_extensive_math(model::Model)
    model = deepcopy(model)
    for reaction in model.reactions
        km = reaction.kinetic_math
        reaction.km = 1.  # PL: Todo: @Anand can you multiply species with `hOSU=true` with their compartment volume?
    end
    reaction
    end

""" Expand reversible reactions to two reactions """
function expand_reversible(model)
    count = 0
    for key in collect(keys(model.reactions))
        reaction = model.reactions[key]
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lb[1] < 0
            count += 1
            # #@info "ORIGINAL" model.reactions[key]
            reverse_reaction = deepcopy(reaction)
            reverse_reaction.lb = (max(0, -reaction.ub[1]), reaction.ub[2])
            reverse_reaction.ub = (-reaction.lb[1], reaction.lb[2])
            reverse_reaction.oc = reaction.oc * -1
            reaction.lb = (max(0, reaction.lb[1]), reaction.lb[2])
            reaction.ub = (max(0, reaction.ub[1]), reaction.ub[2])
            for (k, v) in reverse_reaction.stoichiometry
                reverse_reaction.stoichiometry[k] = v * -1
            end
            model.reactions[key * "_REVERSE"]  = reverse_reaction
            # #@info "REVERSE" model.reactions[key * "_REVERSE"]
        end
    end
    # @info "COUNT ADDED" count
    model  # Todo: convert all Reactions that are `reversible=true` to a forward and reverse reaction with `reversible=false`.
end

""" Get dictonary to change types in kineticLaw """
function _get_substitutions(model)
    subsdict = Dict()
    for k in keys(model.species)
        push!(subsdict, Pair(Num(Variable(Symbol(k))), create_var(k)))
    end
    for k in keys(model.parameters)
        push!(subsdict, Pair(Num(Variable(Symbol(k))), create_param(k)))
    end
    for k in keys(model.compartments)
        push!(subsdict, Pair(Num(Variable(Symbol(k))), create_param(k)))
    end
    subsdict
end

function to_symbol_string(k::MathVal)
    return string(k.val)
end

function to_symbol_string(k::MathIdent)
    return "Num(Variable(Symbol($k.id)))"
end

function to_symbol_string(k::MathApply)
    return "(" * to_symbol_string(k.args[1]) * k.fn * to_symbol_string(k.args[2]) * ")"
end

function make_symbolic(k::Math)
    return eval(Meta.parse(to_symbol_string(k)))
end

""" Convert SBML.Reaction to MTK.Reaction """
function mtk_reactions(model::Model)
    rxs = []
    count = 0
    for reaction in values(model.reactions)
        # #@info "REACTION" reaction
        reactants = Num[]
        rstoich = Num[]
        products = Num[]
        pstoich = Num[]
        for (k, v) in reaction.stoichiometry
            if v < 0
                push!(reactants, create_var(k))
                push!(rstoich, -v)
            elseif v > 0
                push!(products, create_var(k))
                push!(pstoich, v)
            else
                @error("Stoichiometry of $k must be non-zero")
            end
            end
        # if (length(reactants) == 0) reactants = nothing; rstoich = nothing end
        # if (length(products) == 0) products = nothing; pstoich = nothing end
        subsdict = _get_substitutions(model)
        # PL: Todo: @Anand: can you convert kinetic_mathto Symbolic expression. Perhaps it would actually better if kinetic Math would be a Symbolics.jl expression rather than of type `Math`? But Mirek wants `Math`, I think.
        symbolic_math = make_symbolic(reaction.kinetic_math)
        kl = substitute(symbolic_math, subsdict)  # PL: Todo: might need conversion of kinetic_math to Symbolic MTK expression
        if (isnothing(reactants) && isnothing(products))
            # @info "EMPTY" reaction
            count += 1
            # continue
        end
        push!(rxs, ModelingToolkit.Reaction(kl, reactants, products, rstoich, pstoich;only_use_rate=true))
    end
    @info "COUNT EMPTY" count
    rxs
    end
    

""" Extract u0map from Model """
function get_u0(model)
    u0map = []
    for (k, v) in model.species
        # println(v)
        push!(u0map, Pair(create_var(k), v.initial_amount[1]))
    end
    u0map
end

""" Extract paramap from Model """
function get_paramap(model)
    paramap = Pair{Num,Float64}[]
    for (k, v) in model.parameters
        push!(paramap, Pair(create_param(k), v))
    end
    for (k, v) in model.compartments
        push!(paramap, Pair(create_param(k), v.size))
    end
    paramap
end

create_var(x) = Num(Variable(Symbol(x)))
# # create_var(x, iv) = Num(Sym{FnType{Tuple{Real}}}(Symbol(x))(Variable(Symbol(iv)))).val
# # create_var(x, iv) = Num(Variable{Symbolics.FnType{Tuple{Any},Real}}(Symbol(x)))(Variable(Symbol(iv)))
function create_param(x)
    p = Sym{Real}(Symbol(x))
    ModelingToolkit.toparam(p)
end
